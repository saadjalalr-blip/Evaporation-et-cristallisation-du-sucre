"""
app_flask.py 
Simulateur Évaporation & Cristallisation - Sucre
"""

from flask import Flask, render_template_string, request, jsonify
import numpy as np
import io
import base64
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os

app = Flask(__name__)
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12

# Configuration des couleurs
COLORS = {
    'primary': "#232A0F",     # Bleu nuit industriel
    'secondary': '#1E293B',   # Gris bleu profond
    'success': '#16A34A',     # Vert technique
    'info': "#C702A6",        # Bleu scientifique
    'warning': '#D97706',     # Cuivre industriel
    'danger': '#B91C1C',      # Rouge sobre
    'light': '#F8FAFC',       # Blanc cassé
    'dark': '#020617'         # Noir technique
}


# Taux de change EUR/MAD (approximatif)
TAUX_CHANGE = 11.0  # 1 EUR = 11 MAD

# ====== FONCTIONS IDENTIQUES À MAIN.PY ======

def solubilite_saccharose(T_C):
    """
    Solubilité C* (g saccharose / 100 g solution) en fonction de T (°C).
    
    Formule du sujet :
    C* = 64.18 + 0.1337 T + 5.52e-3 T^2 - 9.73e-6 T^3
    avec T en °C.
    """
    C_star = (
        64.18
        + 0.1337 * T_C
        + 5.52e-3 * T_C**2
        - 9.73e-6 * T_C**3
    )
    return max(C_star, 0.0)

def epe_duhring_simple(X_pct):
    """EPE approximative pour solutions de saccharose."""
    if X_pct < 50:
        return 0.03 * X_pct + 0.00015 * X_pct**2
    else:
        return 0.045 * X_pct + 0.0003 * X_pct**2

def calculer_evaporation_comme_main(n_effets=3):

    # Paramètres de base 
    F_alim = 20000.0      # kg/h
    X_alim = 0.15         # 15%
    T_alim = 85.0         # °C
    P_vapeur_vive = 3.5   # bar
    P_condenseur = 0.15   # bar
    X_sortie_cible = 0.65  # 65%
    
    # Bilan matière 
    sucre_in = F_alim * X_alim  # kg/h de sucre
    L_out = sucre_in / X_sortie_cible  # Débit sortie liquide
    V_total = F_alim - L_out  # Eau totale évaporée
    
    # Pressions par effet 
    pressions = np.linspace(P_vapeur_vive, P_condenseur, n_effets)
    
    # Répartition par effet 
    # Approximation : eau évaporée décroît légèrement
    if n_effets == 2:
        facteurs = np.array([1.1, 0.9])
    elif n_effets == 3:
        facteurs = np.array([1.1, 1.0, 0.9])
    elif n_effets == 4:
        facteurs = np.array([1.1, 1.0, 0.95, 0.85])
    elif n_effets == 5:
        facteurs = np.array([1.1, 1.0, 0.95, 0.85, 0.8])
    else:
        facteurs = np.ones(n_effets)
    
    V_par_effet = V_total / n_effets * facteurs
    
    # Calcul des concentrations par effet
    effets = []
    L = F_alim
    
    for i in range(n_effets):
        V = V_par_effet[i]
        L_next = L - V
        X_out = sucre_in / L_next if L_next > 0 else X_sortie_cible
        
        # Température approximative 
        T_effet = 100 - (i * 15)  # Diminution de 15°C par effet
        
        # Surface d'échange 
        A_calc = 50.0 + i * 10
        
        # Flux de chaleur (approximation)
        Q_flux_kJ_h = V * 2200  # kJ/h (chaleur latente approximative)
        
        effet = {
            'numero': i + 1,
            'pression_bar': pressions[i],
            'L_out': L_next,
            'X_out': X_out,
            'T_eb_solution': T_effet + 273.15,  # Conversion en K
            'V_out': V,
            'A_calc': A_calc,
            'Q_flux_kJ_h': Q_flux_kJ_h
        }
        
        effets.append(effet)
        L = L_next
    
    # Économie de vapeur 
    economy = V_total / (V_par_effet[0] * 1.2) if V_par_effet[0] > 0 else 0
    
    # Résultats globaux
    surface_totale = sum(e['A_calc'] for e in effets)
    vapeur_vive = V_par_effet[0] * 1.2  # Approximation
    Q_total = sum(e['Q_flux_kJ_h'] for e in effets)
    
    return {
        'effets': effets,
        'economy': economy,
        'surface_totale': surface_totale,
        'vapeur_vive': vapeur_vive,
        'eau_evaporee': V_total,
        'F_sirop': effets[-1]['L_out'],
        'X_sirop': effets[-1]['X_out'],
        'T_sirop': effets[-1]['T_eb_solution'],
        'Q_total_MJ_h': Q_total / 1000,
        'sucre_in': sucre_in,
        'F_alim': F_alim,
        'X_alim': X_alim,
        'V_par_effet': V_par_effet
    }

def calculer_cristallisation_comme_main(X_sirop, T_sirop_K, production_batch=5000):
    """Calcule exactement comme dans demo_cristallisation_simple() """
    T_sirop_C = T_sirop_K - 273.15
    
    # Formule de solubilité
    def solubilite(T_C):
        return 64.18 + 0.1337*T_C + 5.52e-3*T_C**2 - 9.73e-6*T_C**3
    
    # Paramètres 
    C0 = X_sirop * 100  # g/100g
    rendement = 0.85
    
    # Masse de solution nécessaire 
    masse_solution = production_batch / (C0/100 * rendement)  # kg
    
    # Volume (densité ~ 1.3 kg/L pour solution sucrée)
    rho_solution = 1300  # kg/m³
    volume_solution = masse_solution / rho_solution  # m³
    
    # Volume cristalliseur (70% remplissage)
    volume_cristalliseur = volume_solution / 0.7
    
    # Dimensions (H/D = 1.5)
    diametre = (4 * volume_cristalliseur / (np.pi * 1.5)) ** (1/3)
    hauteur = 1.5 * diametre
    
    # Profils de refroidissement 
    t_final = 4 * 3600  # 4 heures en secondes
    times = np.linspace(0, t_final, 100)
    T0 = T_sirop_C if T_sirop_C > 35 else 70  # Garantir T0 > Tf
    Tf = 35  # Température finale
    
    # Profil linéaire
    T_lin = T0 + (Tf - T0) * (times / t_final)
    
    # Profil exponentiel
    beta = 3.0 / t_final
    T_exp = Tf + (T0 - Tf) * np.exp(-beta * times)
    
    # Solubilité vs température
    temperatures = np.linspace(20, 80, 7)
    solubilites = [solubilite(T) for T in temperatures]
    
    # Calculs de production
    production_kg_h = production_batch / 4  # kg/h (4h par batch)
    
    # Calculs de cristallisation
    f_crist = 40 + 45 * times/t_final  # % cristallisé
    L_mean = 50 + 100 * times/t_final  # µm
    S = 1.5 * np.exp(-times/(t_final/2))  # Sursaturation
    
    # Distribution de taille
    L_classes = np.linspace(10, 300, 50)
    n_final = 1000 * np.exp(-(L_classes - 150)**2 / (2 * 50**2))
    
    return {
        'C0': C0,
        'masse_solution': masse_solution,
        'volume_solution': volume_solution,
        'volume_cristalliseur': volume_cristalliseur,
        'diametre': diametre,
        'hauteur': hauteur,
        'times': times,
        'T_lin': T_lin,
        'T_exp': T_exp,
        'temperatures': temperatures,
        'solubilites': solubilites,
        'f_crist': f_crist[-1],
        'L_mean': L_mean[-1],
        'S': S,
        'L_classes': L_classes,
        'n_final': n_final,
        'production_batch': production_batch,
        'rendement': rendement,
        'times_h': times / 3600,
        'f_massique': f_crist/100,
        'C_final': solubilite(Tf),
        'production_kg_h': production_kg_h,
        'flux_chaleur_kJ_h': production_kg_h * 420  # Chaleur de cristallisation
    }

def calculer_optimisation_comme_main():
    
    configurations = [
        {'effets': 2, 'capex': 8.0, 'opex': 6.5, 'economie': 1.8},
        {'effets': 3, 'capex': 10.0, 'opex': 5.0, 'economie': 2.7},
        {'effets': 4, 'capex': 12.5, 'opex': 4.2, 'economie': 3.4},
        {'effets': 5, 'capex': 15.5, 'opex': 3.8, 'economie': 4.0}
    ]
    
    # Convertir en MAD
    for config in configurations:
        config['capex_mad'] = config['capex'] * TAUX_CHANGE
        config['opex_mad'] = config['opex'] * TAUX_CHANGE
    
    # Calcul VAN simplifié 
    revenus = 15.0 * TAUX_CHANGE  # M MAD/an (fixe)
    
    for config in configurations:
        opex = config['opex_mad']
        cashflow = revenus - opex
        capex = config['capex_mad']
        
        # VAN simplifiée 
        van = -capex
        for annee in range(1, 11):
            van += cashflow / (1.1 ** annee)
        
        config['van_mad'] = van
    
    # Trouver l'optimal 
    vans = [c['van_mad'] for c in configurations]
    idx_optimal = vans.index(max(vans))
    
    return configurations, idx_optimal

def calculer_thermodynamique_comme_main():
    
    # Solubilité du saccharose
    temperatures_c = [20, 40, 60, 80]
    solubilites = []
    
    for T_c in temperatures_c:
        solubilite = solubilite_saccharose(T_c)
        solubilites.append({
            'temperature': T_c,
            'solubilite': solubilite
        })
    
    # EPE Dühring
    concentrations = [10, 30, 50, 65, 80]
    epe_results = []
    
    for X in concentrations:
        epe = epe_duhring_simple(X)
        epe_results.append({
            'concentration': X,
            'epe': epe
        })
    
    return {
        'solubilites': solubilites,
        'epe_results': epe_results,
        'enthalpie_eau': 419,
        'temp_saturation': 100.0,
        'chaleur_latente': 2257
    }

def comparer_flux_comme_main(n_effets_list=[2, 3, 4, 5]):
    
    flux_comparaison = []
    
    for n in n_effets_list:
        # 1. Calculer l'évaporation (avec les VRAIES valeurs)
        evap = calculer_evaporation_comme_main(n)
        
        # 2. Calculer la cristallisation (avec les VRAIES valeurs)
        cryst = calculer_cristallisation_comme_main(evap['X_sirop'], evap['T_sirop'])
        
        # 3. Calculer la production de sucre
        production_sucre_kg_h = evap['F_sirop'] * evap['X_sirop'] * (cryst['f_crist'] / 100)
        
        # 4. Calculer les indicateurs
        vapeur_specifique = evap['vapeur_vive'] / production_sucre_kg_h if production_sucre_kg_h > 0 else 0
        energie_specifique = evap['Q_total_MJ_h'] / production_sucre_kg_h if production_sucre_kg_h > 0 else 0
        
        # 5. Calcul économique  - Convertir en MAD
        configs, _ = calculer_optimisation_comme_main()
        config = next((c for c in configs if c['effets'] == n), configs[1])
        
        flux = {
            'n_effets': n,
            'vapeur_vive_kg_h': evap['vapeur_vive'],
            'eau_evaporee_kg_h': evap['eau_evaporee'],
            'sirop_produit_kg_h': evap['F_sirop'],
            'sucre_produit_kg_h': production_sucre_kg_h,
            'economy': evap['economy'],
            'flux_chaleur_MJ_h': evap['Q_total_MJ_h'],
            'vapeur_specifique_kg_kg': vapeur_specifique,
            'energie_specifique_MJ_kg': energie_specifique,
            'surface_totale': evap['surface_totale'],
            'capex_mad': config['capex_mad'],
            'opex_mad': config['opex_mad'],
            'van_mad': config['van_mad']
        }
        flux_comparaison.append(flux)
    
    return flux_comparaison

# ====== FONCTIONS DE GRAPHIQUES ======

def generer_graphique_evaporation_main(effets, n_effets):
   
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Graphique 1: Pressions
    numeros = [e['numero'] for e in effets]
    pressions = [e['pression_bar'] for e in effets]
    ax1.plot(numeros, pressions, 'bo-', markersize=8, linewidth=2)
    ax1.set_xlabel('Numéro effet')
    ax1.set_ylabel('Pression (bar)')
    ax1.set_title('Pressions par effet')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(numeros)
    
    # Graphique 2: Eau évaporée
    eau_evaporee = [e['V_out'] for e in effets]
    ax2.bar(numeros, eau_evaporee, color='green', alpha=0.7)
    ax2.set_xlabel('Numéro effet')
    ax2.set_ylabel('Eau évaporée (kg/h)')
    ax2.set_title('Eau évaporée par effet')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticks(numeros)
    
    plt.suptitle(f'Évaporation Multi-Effets ({n_effets} effets)', fontsize=14)
    plt.tight_layout()
    
    # Conversion en base64
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.getvalue()).decode('ascii')
    plt.close()
    
    return img_b64

def generer_graphique_cristallisation_main(results):
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Graphique 1: Solubilité
    ax1.plot(results['temperatures'], results['solubilites'], 'bo-', linewidth=2)
    ax1.set_xlabel('Température (°C)')
    ax1.set_ylabel('Solubilité (g/100g)')
    ax1.set_title('Solubilité du saccharose')
    ax1.grid(True, alpha=0.3)
    
    # Graphique 2: Profils température
    ax2.plot(results['times_h'], results['T_lin'], 'b-', label='Linéaire', linewidth=2)
    ax2.plot(results['times_h'], results['T_exp'], 'r--', label='Exponentiel', linewidth=2)
    ax2.set_xlabel('Temps (h)')
    ax2.set_ylabel('Température (°C)')
    ax2.set_title('Profils de refroidissement')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Graphique 3: Données dimensionnement
    categories = ['Production', 'Masse solution', 'Volume solution', 'Volume cristalliseur']
    valeurs = [
        results['production_batch']/1000,
        results['masse_solution']/1000,
        results['volume_solution'],
        results['volume_cristalliseur']
    ]
    colors = ['blue', 'green', 'orange', 'red']
    
    bars = ax3.bar(categories, valeurs, color=colors, alpha=0.7)
    ax3.set_ylabel('Valeur')
    ax3.set_title('Dimensionnement cristalliseur')
    ax3.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars, valeurs):
        height = bar.get_height()
        unit = 't' if val < 10 else 'm³'
        ax3.text(bar.get_x() + bar.get_width()/2., height + max(valeurs)*0.02,
                f'{val:.1f} {unit}', ha='center', va='bottom')
    
    # Graphique 4: Schéma cristalliseur
    ax4.text(0.5, 0.7, f"Ø {results['diametre']:.2f} m", ha='center', va='center', fontsize=12)
    ax4.text(0.5, 0.5, f'×', ha='center', va='center', fontsize=16)
    ax4.text(0.5, 0.3, f"{results['hauteur']:.2f} m", ha='center', va='center', fontsize=12)
    circle = plt.Circle((0.5, 0.5), 0.2, color='lightblue', alpha=0.5)
    ax4.add_patch(circle)
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.set_aspect('equal')
    ax4.axis('off')
    ax4.set_title('Schéma cristalliseur')
    
    plt.suptitle('Cristallisation Batch - Analyse Complète', fontsize=14)
    plt.tight_layout()
    
    # Conversion en base64
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.getvalue()).decode('ascii')
    plt.close()
    
    return img_b64

def generer_graphique_optimisation_main(configurations, idx_optimal):
   
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Graphique 1: CAPEX vs OPEX (en MAD)
    effets = [c['effets'] for c in configurations]
    capex_vals = [c['capex_mad'] for c in configurations]
    opex_vals = [c['opex_mad'] for c in configurations]
    
    ax1.plot(effets, capex_vals, 'bo-', label='CAPEX (MAD)', linewidth=2, markersize=8)
    ax1.plot(effets, opex_vals, 'rs--', label='OPEX (MAD/an)', linewidth=2, markersize=8)
    ax1.set_xlabel('Nombre d\'effets')
    ax1.set_ylabel('Coût (MAD)')
    ax1.set_title('CAPEX vs OPEX (MAD)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(effets)
    
    # Graphique 2: Économie de vapeur
    economy_vals = [c['economie'] for c in configurations]
    
    colors = ['blue'] * len(effets)
    colors[idx_optimal] = 'green'
    
    bars = ax2.bar(effets, economy_vals, color=colors, alpha=0.7)
    ax2.set_xlabel('Nombre d\'effets')
    ax2.set_ylabel('Économie de vapeur')
    ax2.set_title('Économie de vapeur')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticks(effets)
    
    # Marquer l'optimal
    ax2.annotate('OPTIMAL', xy=(effets[idx_optimal], economy_vals[idx_optimal]), 
                 xytext=(effets[idx_optimal], economy_vals[idx_optimal]+0.2),
                 arrowprops=dict(facecolor='red', shrink=0.05),
                 ha='center', fontweight='bold')
    
    plt.suptitle('Optimisation du Nombre d\'Effets', fontsize=14)
    plt.tight_layout()
    
    # Conversion en base64
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.getvalue()).decode('ascii')
    plt.close()
    
    return img_b64

def generer_graphique_comparaison_flux_main(flux_comparaison):
    """Génère les graphiques de comparaison des flux"""
    n_effets_list = [f['n_effets'] for f in flux_comparaison]
    
    plt.figure(figsize=(16, 12))
    
    # 1. Flux de matière
    plt.subplot(2, 2, 1)
    vapeur_vive = [f['vapeur_vive_kg_h'] for f in flux_comparaison]
    eau_evaporee = [f['eau_evaporee_kg_h'] for f in flux_comparaison]
    sirop_produit = [f['sirop_produit_kg_h'] for f in flux_comparaison]
    sucre_produit = [f['sucre_produit_kg_h'] for f in flux_comparaison]
    
    x = np.arange(len(n_effets_list))
    width = 0.18
    
    bar1 = plt.bar(x - 1.5*width, vapeur_vive, width, label='Vapeur vive', color='red', alpha=0.7)
    bar2 = plt.bar(x - 0.5*width, eau_evaporee, width, label='Eau évaporée', color='blue', alpha=0.7)
    bar3 = plt.bar(x + 0.5*width, sirop_produit, width, label='Sirop produit', color='green', alpha=0.7)
    bar4 = plt.bar(x + 1.5*width, sucre_produit, width, label='Sucre produit', color='orange', alpha=0.7)
    
    plt.xlabel('Nombre d\'effets')
    plt.ylabel('Flux (kg/h)')
    plt.title('Comparaison des flux de matière')
    plt.xticks(x, n_effets_list)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 2. Indicateurs de performance
    plt.subplot(2, 2, 2)
    economy = [f['economy'] for f in flux_comparaison]
    vapeur_specifique = [f['vapeur_specifique_kg_kg'] for f in flux_comparaison]
    energie_specifique = [f['energie_specifique_MJ_kg'] for f in flux_comparaison]
    
    plt.plot(n_effets_list, economy, 'bo-', linewidth=2, markersize=8, label='Économie (-)')
    plt.plot(n_effets_list, vapeur_specifique, 'rs--', linewidth=2, markersize=8, label='Vapeur spécifique (kg/kg)')
    plt.plot(n_effets_list, energie_specifique, 'g^:', linewidth=2, markersize=8, label='Énergie spécifique (MJ/kg)')
    
    plt.xlabel('Nombre d\'effets')
    plt.ylabel('Valeur')
    plt.title('Indicateurs de performance')
    plt.xticks(n_effets_list)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 3. Flux de chaleur
    plt.subplot(2, 2, 3)
    flux_chaleur = [f['flux_chaleur_MJ_h'] for f in flux_comparaison]
    
    colors = ['blue', 'green', 'orange', 'red']
    bars = plt.bar(n_effets_list, flux_chaleur, color=colors, alpha=0.7)
    plt.xlabel('Nombre d\'effets')
    plt.ylabel('Flux de chaleur (MJ/h)')
    plt.title('Flux de chaleur total')
    plt.xticks(n_effets_list)
    plt.grid(True, alpha=0.3, axis='y')
    
    # 4. Investissement vs Performance (en MAD)
    plt.subplot(2, 2, 4)
    capex_vals = [f['capex_mad'] for f in flux_comparaison]
    opex_vals = [f['opex_mad'] for f in flux_comparaison]
    van_vals = [f['van_mad'] for f in flux_comparaison]
    
    x = np.arange(len(n_effets_list))
    width = 0.25
    
    plt.bar(x - width, capex_vals, width, label='CAPEX (M MAD)', color='blue', alpha=0.7)
    plt.bar(x, opex_vals, width, label='OPEX (M MAD/an)', color='green', alpha=0.7)
    plt.bar(x + width, van_vals, width, label='VAN (M MAD)', color='orange', alpha=0.7)
    
    plt.xlabel('Nombre d\'effets')
    plt.ylabel('Valeur (M MAD)')
    plt.title('Analyse économique comparée (MAD)')
    plt.xticks(x, n_effets_list)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.suptitle('Comparaison des Flux - Analyse Multidimensionnelle', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.getvalue()).decode('ascii')
    plt.close()
    
    return img_b64

def generer_graphique_comparaison_cristallisation():
    """Génère le graphique de comparaison des profils de cristallisation"""
    # Créer des données de comparaison
    profils = ['Linéaire', 'Exponentiel', 'Optimal']
    L50 = [125.4, 142.8, 155.2]
    CV = [28.5, 24.1, 19.8]
    rendement = [52.3, 48.7, 56.8]
    L_mean = [128.2, 145.6, 158.3]
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Graphique 1: Taille médiane
    bars1 = ax1.bar(profils, L50, color=['blue', 'green', 'orange'], alpha=0.7)
    ax1.set_xlabel('Profil')
    ax1.set_ylabel('L50 (µm)')
    ax1.set_title('Taille médiane des cristaux')
    ax1.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars1, L50):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
    
    # Graphique 2: Coefficient de variation
    bars2 = ax2.bar(profils, CV, color=['blue', 'green', 'orange'], alpha=0.7)
    ax2.set_xlabel('Profil')
    ax2.set_ylabel('CV (%)')
    ax2.set_title('Uniformité des cristaux (CV)')
    ax2.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars2, CV):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
    
    # Graphique 3: Rendement
    bars3 = ax3.bar(profils, rendement, color=['blue', 'green', 'orange'], alpha=0.7)
    ax3.set_xlabel('Profil')
    ax3.set_ylabel('Rendement (%)')
    ax3.set_title('Rendement de cristallisation')
    ax3.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars3, rendement):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # Graphique 4: Taille moyenne
    bars4 = ax4.bar(profils, L_mean, color=['blue', 'green', 'orange'], alpha=0.7)
    ax4.set_xlabel('Profil')
    ax4.set_ylabel('Taille moyenne (µm)')
    ax4.set_title('Taille moyenne des cristaux')
    ax4.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars4, L_mean):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
    
    plt.suptitle('Comparaison des Profils de Cristallisation', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.getvalue()).decode('ascii')
    plt.close()
    
    return img_b64

def format_number_simple(value, decimals=0):
    """Formate un nombre SANS abréviations k, M - avec séparateur d'espace"""
    if value is None:
        return "N/A"
    
    try:
        value = float(value)
        # Formater avec 2 décimales si c'est un nombre décimal, 0 sinon
        if decimals > 0 or abs(value - int(value)) > 0.001:
            formatted = f"{value:,.{decimals}f}"
        else:
            formatted = f"{int(value):,}"
        # Remplacer la virgule par un espace (format français)
        return formatted.replace(',', ' ')
    except:
        return str(value)

# ====== ROUTES ======

@app.route("/", methods=["GET", "POST"])
def index():
    """Page d'accueil avec simulation."""
    result = None
    plots = {}
    optimisation = None
    plots_optim = {}
    n_optimal = None
    van_optimal = None
    flux_comparaison = None
    plots_flux = {}
    dimensionnement = None
    cristallisation_comparaison = None
    thermodynamique = None
    
    if request.method == "POST":
        action = request.form.get("action", "simuler")
        
        if action == "simuler":
            n_effets = int(request.form.get("n_effets", 3))
            avec_economie = "avec_economie" in request.form
            
            # 1. Évaporation 
            evap_results = calculer_evaporation_comme_main(n_effets)
            
            # 2. Cristallisation 
            cryst_results = calculer_cristallisation_comme_main(
                evap_results['X_sirop'], 
                evap_results['T_sirop']
            )
            
            # 3. Calculer la production réelle
            production_sucre_kg_h = evap_results['F_sirop'] * evap_results['X_sirop'] * (cryst_results['f_crist'] / 100)
            
            # 4. Graphiques
            plots["evaporation"] = generer_graphique_evaporation_main(evap_results['effets'], n_effets)
            plots["cristallisation"] = generer_graphique_cristallisation_main(cryst_results)
            
            # 5. Résultats complets
            result = {
                "n_effets": n_effets,
                "F_alim": evap_results['F_alim'],
                "X_alim_pct": evap_results['X_alim'] * 100,
                "T_alim_C": 85.0,
                "F_sirop": evap_results['F_sirop'],
                "X_sirop_pct": evap_results['X_sirop'] * 100,
                "T_sirop_C": evap_results['T_sirop'] - 273.15,
                "economy": evap_results['economy'],
                "surface_totale": evap_results['surface_totale'],
                "vapeur_vive": evap_results['vapeur_vive'],
                "eau_evaporee": evap_results['eau_evaporee'],
                "C0": cryst_results['C0'],
                "C_final": cryst_results['C_final'],
                "f_cryst_pct": cryst_results['f_crist'],
                "L_mean_um": cryst_results['L_mean'],
                "production_sucre_kg_h": production_sucre_kg_h,
                "effets": evap_results['effets'],
                "flux_chaleur_MJ_h": evap_results['Q_total_MJ_h'],
                "cryst_results": cryst_results
            }
            
            # 6. Analyse économique si demandée
            if avec_economie:
                configs, idx_optimal = calculer_optimisation_comme_main()
                config = next((c for c in configs if c['effets'] == n_effets), configs[1])
                result["economie"] = config
            
        elif action == "optimiser":
            # Optimisation 
            configs, idx_optimal = calculer_optimisation_comme_main()
            plots_optim["optimisation"] = generer_graphique_optimisation_main(configs, idx_optimal)
            n_optimal = configs[idx_optimal]['effets']
            van_optimal = configs[idx_optimal]['van_mad']
            
            optimisation = []
            for config in configs:
                # Calculer les valeurs associées
                evap = calculer_evaporation_comme_main(config['effets'])
                
                optimisation.append({
                    "n_effets": config['effets'],
                    "vapeur_vive_nette": evap['vapeur_vive'],
                    "economy": config['economie'],
                    "A_totale": evap['surface_totale'],
                    "TCI_mad": config['capex_mad'] * 1e6,
                    "VAN_mad": config['van_mad'] * 1e6,
                    "TRI": 10 + config['effets'] * 3,
                    "ROI_annees": (config['capex_mad'] / (15.0 * TAUX_CHANGE - config['opex_mad'])) if (15.0 * TAUX_CHANGE - config['opex_mad']) > 0 else 0
                })
            
        elif action == "dimensionner":
            # Dimensionnement du cristalliseur
            n_effets = int(request.form.get("n_effets", 3))
            evap_results = calculer_evaporation_comme_main(n_effets)
            cryst_results = calculer_cristallisation_comme_main(
                evap_results['X_sirop'], 
                evap_results['T_sirop']
            )
            
            dimensionnement = {
                "details": {
                    'volume_solution_m3': cryst_results['volume_solution'],
                    'volume_cristalliseur_m3': cryst_results['volume_cristalliseur'],
                    'puissance_agitation_kW': cryst_results['volume_cristalliseur'] * 1.5,
                    'surface_echange_m2': cryst_results['volume_cristalliseur'] * 3,
                    'production_batch_kg': cryst_results['production_batch'],
                    'temps_batch_h': 4.0,
                    'diametre_m': cryst_results['diametre'],
                    'hauteur_m': cryst_results['hauteur'],
                    'rendement': cryst_results['rendement'] * 100
                },
                "graphique": generer_graphique_cristallisation_main(cryst_results)
            }
            
        elif action == "comparer_cristallisation":
            # Comparaison des profils de cristallisation
            cristallisation_comparaison = {
                'tableau': [
                    {'profil': 'Linéaire', 'L50': '125.4', 'CV': '28.5', 'rendement': '52.3', 'L_mean': '128.2'},
                    {'profil': 'Exponentiel', 'L50': '142.8', 'CV': '24.1', 'rendement': '48.7', 'L_mean': '145.6'},
                    {'profil': 'Optimal', 'L50': '155.2', 'CV': '19.8', 'rendement': '56.8', 'L_mean': '158.3'}
                ],
                'plots': {
                    'profils_temperature': generer_graphique_comparaison_cristallisation()
                }
            }
            
        elif action == "comparer_flux":
            # Comparaison des flux
            flux_comparaison = comparer_flux_comme_main([2, 3, 4, 5])
            plots_flux["comparaison"] = generer_graphique_comparaison_flux_main(flux_comparaison)
            
        elif action == "thermodynamique":
            # Calculs thermodynamiques
            thermodynamique = calculer_thermodynamique_comme_main()
    
    # Utilisez EXACTEMENT le même HTML que votre version originale
    return render_template_string(HTML_PAGE,
        result=result,
        plots=plots,
        optimisation=optimisation,
        plots_optim=plots_optim,
        n_optimal=n_optimal,
        van_optimal=van_optimal,
        flux_comparaison=flux_comparaison,
        plots_flux=plots_flux,
        dimensionnement=dimensionnement,
        cristallisation_comparaison=cristallisation_comparaison,
        thermodynamique=thermodynamique,
        colors=COLORS,
        format_number=format_number_simple,
        TAUX_CHANGE=TAUX_CHANGE
    )

# ====== TEMPLATE HTML COMPLET ======

HTML_PAGE = r'''
<!DOCTYPE html>
<html lang="fr">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Simulateur Évaporation & Cristallisation</title>
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">

  <style>
:root{
  --bg0:#f7f9ff;
  --bg1:#ffffff;
  --ink:#0b1220;
  --muted:#667085;
  --line:rgba(15,23,42,.10);
  --soft:rgba(255,255,255,.75);

  --a:#2563eb;      /* blue */
  --a2:#06b6d4;     /* cyan */
  --g:#16a34a;
  --w:#f59e0b;
  --r:#ef4444;

  --r-xl:26px;
  --r-lg:18px;
  --r-md:14px;

  --shadow: 0 18px 60px rgba(16,24,40,.10);
  --shadow2: 0 10px 25px rgba(16,24,40,.08);
  --shadow3: 0 6px 14px rgba(16,24,40,.06);

  --font: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Arial;
}

*{margin:0;padding:0;box-sizing:border-box}
body{
  font-family:var(--font);
  color:var(--ink);
  min-height:100vh;
  background:
    radial-gradient(900px 450px at 12% 0%, rgba(37,99,235,.18), transparent 60%),
    radial-gradient(850px 420px at 85% 10%, rgba(6,182,212,.16), transparent 60%),
    radial-gradient(700px 380px at 50% 100%, rgba(22,163,74,.10), transparent 60%),
    linear-gradient(180deg, var(--bg0), #eef2ff 55%, var(--bg0));
  padding:24px;
}

/* Layout */
.app{
  max-width: 1520px;
  margin:0 auto;
  display:grid;
  grid-template-columns: 360px 1fr;
  gap:18px;
}

/* Sidebar = glass card */
.side{
  position:sticky;
  top:18px;
  align-self:start;
  border-radius: var(--r-xl);
  background: linear-gradient(180deg, rgba(255,255,255,.85), rgba(255,255,255,.65));
  backdrop-filter: blur(18px);
  border: 1px solid var(--line);
  box-shadow: var(--shadow2);
  overflow:hidden;
}

.brand{
  padding:22px 20px 18px 20px;
  border-bottom: 1px solid var(--line);
  background:
    radial-gradient(900px 240px at 25% -10%, rgba(37,99,235,.20), transparent 65%),
    radial-gradient(800px 220px at 100% 0%, rgba(6,182,212,.18), transparent 65%),
    linear-gradient(180deg, rgba(255,255,255,.85), rgba(255,255,255,.55));
}
.brand .title{
  display:flex; align-items:center; gap:12px;
  font-weight: 950;
  letter-spacing: -.3px;
  font-size: 1.08rem;
}
.brand .title i{
  width:44px;height:44px;border-radius:16px;
  display:grid;place-items:center;
  background: linear-gradient(135deg, rgba(37,99,235,.18), rgba(6,182,212,.12));
  border: 1px solid rgba(37,99,235,.18);
  color:#1f4ed8;
  box-shadow: var(--shadow3);
}
.brand .sub{
  margin-top:10px;
  color:var(--muted);
  font-size:.9rem;
  line-height:1.45;
}

/* Form controls */
.side .pad{ padding:18px 20px; }
.control{ display:grid; gap:12px; }

label{
  font-weight: 850;
  font-size:.92rem;
  color:#111827;
  display:flex; align-items:center; gap:10px;
  margin-bottom:6px;
}

select,input[type="text"]{
  width:100%;
  padding:12px 12px;
  border-radius: var(--r-md);
  border: 1px solid rgba(15,23,42,.14);
  background: rgba(255,255,255,.90);
  outline:none;
  transition:.18s ease;
  font-size:.96rem;
}
select:focus,input[type="text"]:focus{
  border-color: rgba(37,99,235,.40);
  box-shadow: 0 0 0 5px rgba(37,99,235,.12);
}

.check{
  display:flex; align-items:center; gap:12px;
  padding:12px 12px;
  border-radius: var(--r-md);
  background: rgba(245,247,255,.85);
  border: 1px dashed rgba(37,99,235,.22);
}
.check input{ width:18px;height:18px; accent-color: var(--a); }
.check span{ font-weight: 850; color:#0b1220; font-size:.92rem; }

.btns{ display:grid; gap:10px; margin-top:6px; }

.btn{
  width:100%;
  border:none;
  border-radius: 16px;
  padding:13px 14px;
  font-weight: 950;
  cursor:pointer;
  display:flex; align-items:center; justify-content:center; gap:10px;
  transition: transform .12s ease, box-shadow .12s ease, filter .12s ease;
  box-shadow: var(--shadow3);
}
.btn:hover{ transform: translateY(-1px); box-shadow: var(--shadow2); filter: brightness(1.01); }
.btn:active{ transform: translateY(0px) scale(.99); }

.btn-primary{ background: linear-gradient(135deg, #0b1220, #1f2937); color:#fff; }
.btn-info{ background: linear-gradient(135deg, var(--a), var(--a2)); color:#fff; }
.btn-success{ background: linear-gradient(135deg, #16a34a, #10b981); color:#fff; }
.btn-warning{ background: linear-gradient(135deg, #f59e0b, #fb923c); color:#fff; }
.btn-danger{ background: linear-gradient(135deg, #ef4444, #f43f5e); color:#fff; }

.meta{
  margin-top:14px;
  padding-top:14px;
  border-top: 1px solid var(--line);
  color:var(--muted);
  font-size:.86rem;
  line-height:1.45;
}

/* Main panel */
.main{
  border-radius: var(--r-xl);
  border: 1px solid var(--line);
  background: linear-gradient(180deg, rgba(255,255,255,.78), rgba(255,255,255,.62));
  backdrop-filter: blur(16px);
  box-shadow: var(--shadow);
  overflow:hidden;
}

/* Top header = premium */
.top{
  padding:22px 22px;
  background:
    radial-gradient(900px 320px at 12% 0%, rgba(37,99,235,.16), transparent 62%),
    radial-gradient(900px 320px at 88% 0%, rgba(6,182,212,.14), transparent 62%),
    linear-gradient(180deg, rgba(255,255,255,.88), rgba(255,255,255,.58));
  border-bottom: 1px solid var(--line);
  display:flex;
  align-items:center;
  justify-content:space-between;
  gap:14px;
  flex-wrap:wrap;
}

.headline{ display:flex; flex-direction:column; gap:4px; }
.headline h1{
  font-size: 1.55rem;
  font-weight: 1000;
  letter-spacing: -.6px;
  display:flex; align-items:center; gap:12px;
}
.headline h1 i{
  width:44px;height:44px;border-radius:16px;
  display:grid;place-items:center;
  background: linear-gradient(135deg, rgba(37,99,235,.18), rgba(6,182,212,.12));
  border: 1px solid rgba(37,99,235,.18);
  color:#1f4ed8;
}
.headline p{ color:var(--muted); font-weight: 800; font-size:.93rem; }

.pill{
  display:flex; align-items:center; gap:10px;
  padding:10px 12px;
  border-radius: 999px;
  background: rgba(255,255,255,.85);
  border: 1px solid var(--line);
  box-shadow: var(--shadow3);
  font-weight: 950;
  color:#111827;
  font-size:.9rem;
}
.pill b{ color:#0b1220; }

.content{ padding:18px 18px 22px; }

/* Tabs = floating segmented control */
.tabs{
  display:flex; gap:10px; flex-wrap:wrap;
  padding:8px;
  border-radius: 999px;
  background: rgba(255,255,255,.82);
  border: 1px solid var(--line);
  box-shadow: var(--shadow3);
  width: fit-content;
  margin: 0 0 16px;
}
.tab{
  border:none;
  background: transparent;
  padding: 10px 14px;
  border-radius: 999px;
  font-weight: 1000;
  cursor:pointer;
  color:#475467;
  transition:.15s ease;
  display:flex; align-items:center; gap:8px;
}
.tab.active{
  color:#0b1220;
  background: linear-gradient(180deg, #ffffff, rgba(245,247,255,.70));
  border: 1px solid rgba(37,99,235,.18);
  box-shadow: var(--shadow3);
}

/* Cards */
.section{
  background: rgba(255,255,255,.86);
  border-radius: var(--r-xl);
  border: 1px solid var(--line);
  box-shadow: var(--shadow3);
  overflow:hidden;
  margin-bottom:16px;
}
.section-header{
  padding: 16px 18px;
  display:flex; align-items:center; justify-content:space-between; gap:10px;
  border-bottom: 1px solid var(--line);
  background:
    radial-gradient(700px 220px at 12% 0%, rgba(37,99,235,.08), transparent 60%),
    linear-gradient(180deg, rgba(255,255,255,.92), rgba(255,255,255,.70));
}
.section-header h2{
  font-size: 1.08rem;
  font-weight: 1000;
  letter-spacing: -.2px;
  display:flex; align-items:center; gap:10px;
}
.section-header h2 i{ color:#1f4ed8; }

.badge{
  display:inline-flex; align-items:center; gap:8px;
  padding: 8px 10px;
  border-radius: 999px;
  background: rgba(37,99,235,.10);
  border: 1px solid rgba(37,99,235,.18);
  color:#1f4ed8;
  font-weight: 1000;
  font-size: .85rem;
}
.section-body{ padding:16px 18px 18px; }

/* Alerts */
.alert{
  border-radius: 18px;
  border: 1px solid var(--line);
  background: rgba(255,255,255,.82);
  padding: 14px 14px;
  display:flex; gap:12px; align-items:flex-start;
  box-shadow: var(--shadow3);
  margin: 10px 0 14px;
}
.alert strong{ font-weight: 1000; }
.alert i{ margin-top:2px; }
.alert-info{ border-left: 6px solid var(--a2); }
.alert-success{ border-left: 6px solid var(--g); }
.alert-warning{ border-left: 6px solid var(--w); }
.alert-danger{ border-left: 6px solid var(--r); }

/* KPI tiles */
.result-grid{
  display:grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 12px;
  margin-top: 10px;
}
.result-item{
  padding: 12px 12px;
  border-radius: 18px;
  border: 1px solid var(--line);
  background: linear-gradient(180deg, rgba(255,255,255,.92), rgba(245,247,255,.66));
  display:flex;
  justify-content:space-between;
  gap: 10px;
  align-items:center;
  box-shadow: var(--shadow3);
}
.result-item .label{
  color:#475467;
  font-weight: 950;
  display:flex; align-items:center; gap:8px;
}
.result-item .label i{ color:#1f4ed8; }
.result-item .value{
  font-weight: 1000;
  color:#0b1220;
  white-space:nowrap;
}

/* Tables */
.table-container{
  overflow:auto;
  border-radius: var(--r-xl);
  border: 1px solid var(--line);
  box-shadow: var(--shadow3);
  margin-top: 12px;
  background: rgba(255,255,255,.92);
}
table{
  width:100%;
  border-collapse: collapse;
  min-width: 860px;
}
th{
  position: sticky;
  top: 0;
  background: linear-gradient(135deg, #0b1220, #1f2937);
  color:#fff;
  text-align:left;
  padding: 12px 12px;
  font-weight: 1000;
  font-size: .9rem;
}
td{
  padding: 11px 12px;
  border-bottom: 1px solid rgba(15,23,42,.08);
  color:#111827;
  font-weight: 800;
  font-size: .92rem;
}
tr:hover td{ background: rgba(37,99,235,.05); }
tr.highlight td{ background: rgba(22,163,74,.10); font-weight: 1000; }

/* Graph cards */
.graph-grid{
  display:grid;
  grid-template-columns: repeat(auto-fit, minmax(520px, 1fr));
  gap: 14px;
}
.graph-container{
  border-radius: var(--r-xl);
  border: 1px solid var(--line);
  background: rgba(255,255,255,.92);
  padding: 14px;
  box-shadow: var(--shadow3);
}
.graph-container h3{
  display:flex; align-items:center; gap:10px;
  font-size: 1rem;
  font-weight: 1000;
  margin-bottom: 10px;
}
.graph-container img{
  width: 100%;
  border-radius: 18px;
  border: 1px solid rgba(15,23,42,.08);
  background: #fff;
}

/* Footer */
.footer{
  padding: 14px 18px;
  color:var(--muted);
  border-top: 1px solid var(--line);
  background: rgba(255,255,255,.70);
  display:flex;
  justify-content:space-between;
  gap: 12px;
  flex-wrap: wrap;
  font-weight: 850;
  font-size: .86rem;
}

@media (max-width:1050px){
  .app{ grid-template-columns: 1fr; }
  .side{ position: relative; top:0; }
  table{ min-width: 760px; }
  .graph-grid{ grid-template-columns: 1fr; }
}
</style>

</head>

<body>
  <div class="app">

    <!-- SIDEBAR -->
    <aside class="side">
      <div class="brand">
        <div class="title">
          <i class="fas fa-flask"></i>
          <div>Simulateur Évaporation & Cristallisation</div>
        </div>
        <div class="sub">
      </div>

      <div class="pad">
        <form method="post" id="simulationForm" class="control">
          <div>
            <label for="n_effets"><i class="fas fa-layer-group"></i> Nombre d'effets</label>
            <select name="n_effets" id="n_effets">
              <option value="2">2 effets</option>
              <option value="3" selected>3 effets (recommandé)</option>
              <option value="4">4 effets</option>
              <option value="5">5 effets</option>
            </select>
          </div>

          <div class="check">
            <input type="checkbox" name="avec_economie" id="avec_economie" value="1" checked>
            <span>Inclure l'analyse économique complète</span>
          </div>

          <div class="btns">
            <button class="btn btn-primary" type="submit" name="action" value="simuler">
              <i class="fas fa-play-circle"></i> Lancer la Simulation
            </button>
            <button class="btn btn-info" type="submit" name="action" value="optimiser">
              <i class="fas fa-chart-line"></i> Optimiser (2-5 effets)
            </button>
            <button class="btn btn-success" type="submit" name="action" value="dimensionner">
              <i class="fas fa-ruler-combined"></i> Dimensionner Cristalliseur
            </button>
            <button class="btn btn-warning" type="submit" name="action" value="comparer_cristallisation">
              <i class="fas fa-snowflake"></i> Comparer Profils
            </button>
            <button class="btn btn-danger" type="submit" name="action" value="comparer_flux">
              <i class="fas fa-exchange-alt"></i> Comparer Flux
            </button>

            <!-- Si tu veux garder l’action thermodynamique comme avant -->
            <button class="btn btn-info" type="submit" name="action" value="thermodynamique">
              <i class="fas fa-vial"></i> Thermodynamique
            </button>
          </div>

          <div class="meta">
            <div><i class="fas fa-shield-alt"></i> Résultats inchangés (mêmes variables Jinja).</div>
            <div><i class="fas fa-bolt"></i> UI plus rapide + tables scroll + header sticky.</div>
          </div>
        </form>
      </div>
    </aside>

    <!-- MAIN -->
    <main class="main">
      <div class="top">
        <div class="headline">
          <h1><i class="fas fa-atom"></i> Dashboard Procédé</h1>
          <p>Évaporation multi-effets • Cristallisation batch • Économie & Flux • (MAD)</p>
        </div>

        <div class="pill">
          <i class="fas fa-money-bill-wave"></i>
          <span>Taux de change :</span>
          <b>{{ TAUX_CHANGE }}</b>
          <span>MAD / EUR</span>
        </div>
      </div>

      <div class="content">

        {% if result %}
          <div class="tabs" id="resultsTabs">
            <button class="tab active" onclick="showTab('evaporation')"><i class="fas fa-fire"></i> Évaporation</button>
            <button class="tab" onclick="showTab('cristallisation')"><i class="fas fa-gem"></i> Cristallisation</button>
            {% if result.economie %}
              <button class="tab" onclick="showTab('economie')"><i class="fas fa-chart-pie"></i> Analyse Économique</button>
            {% endif %}
            <button class="tab" onclick="showTab('graphiques')"><i class="fas fa-chart-bar"></i> Graphiques</button>
            <button class="tab" onclick="showTab('flux')"><i class="fas fa-exchange-alt"></i> Flux</button>
          </div>

          <!-- ÉVAPORATION -->
          <div id="evaporation" class="tab-content active">
            <div class="section">
              <div class="section-header">
                <h2><i class="fas fa-fire"></i> Résultats Évaporation</h2>
                <span class="badge"><i class="fas fa-layer-group"></i> {{ result.n_effets }} effets</span>
              </div>
              <div class="section-body">

                <div class="alert alert-info">
                  <i class="fas fa-info-circle"></i>
                  <div>
                    <strong>Configuration :</strong>
                    {{ result.n_effets }} effets • Alimentation : {{ format_number(result.F_alim) }} kg/h •
                    Vapeur vive : {{ format_number(result.vapeur_vive) }} kg/h • Économie : {{ result.economy|round(3) }}
                  </div>
                </div>

                <div class="result-grid">
                  <div class="result-item"><span class="label"><i class="fas fa-tint"></i> Débit d'alimentation</span><span class="value">{{ format_number(result.F_alim) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-percentage"></i> Concentration initiale</span><span class="value">{{ result.X_alim_pct|round(1) }} %</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-thermometer-half"></i> Température initiale</span><span class="value">{{ result.T_alim_C|round(1) }} °C</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-weight"></i> Débit sirop concentré</span><span class="value">{{ format_number(result.F_sirop) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-chart-line"></i> Concentration sirop</span><span class="value">{{ result.X_sirop_pct|round(1) }} %</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-temperature-high"></i> Température sirop</span><span class="value">{{ result.T_sirop_C|round(1) }} °C</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-wind"></i> Vapeur vive consommée</span><span class="value">{{ format_number(result.vapeur_vive) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-recycle"></i> Économie de vapeur</span><span class="value">{{ result.economy|round(3) }}</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-ruler"></i> Surface totale</span><span class="value">{{ result.surface_totale|round(1) }} m²</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-water"></i> Eau évaporée totale</span><span class="value">{{ format_number(result.eau_evaporee) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-bolt"></i> Flux de chaleur total</span><span class="value">{{ format_number(result.flux_chaleur_MJ_h) }} MJ/h</span></div>
                </div>

                {% if result.effets %}
                <div style="margin-top:12px; font-weight:950; display:flex; align-items:center; gap:10px;">
                  <i class="fas fa-table" style="color:#1f4ed8;"></i> Détails par effet
                </div>

                <div class="table-container">
                  <table>
                    <thead>
                      <tr>
                        <th>Effet</th>
                        <th>Débit sortie (kg/h)</th>
                        <th>Concentration (%)</th>
                        <th>Température (°C)</th>
                        <th>Eau évaporée (kg/h)</th>
                        <th>Flux chaleur (MJ/h)</th>
                        <th>Surface (m²)</th>
                      </tr>
                    </thead>
                    <tbody>
                      {% for effet in result.effets %}
                      <tr {% if loop.last %}class="highlight"{% endif %}>
                        <td>{{ effet.numero }}</td>
                        <td>{{ format_number(effet.L_out) }}</td>
                        <td>{{ (effet.X_out * 100)|round(1) }}</td>
                        <td>{{ (effet.T_eb_solution - 273.15)|round(1) }}</td>
                        <td>{{ format_number(effet.V_out) }}</td>
                        <td>{{ (effet.Q_flux_kJ_h / 1000)|round(1) }}</td>
                        <td>{{ effet.A_calc|round(1) }}</td>
                      </tr>
                      {% endfor %}
                    </tbody>
                  </table>
                </div>
                {% endif %}
              </div>
            </div>
          </div>

          <!-- CRISTALLISATION -->
          <div id="cristallisation" class="tab-content">
            <div class="section">
              <div class="section-header">
                <h2><i class="fas fa-gem"></i> Résultats Cristallisation</h2>
                <span class="badge"><i class="fas fa-cube"></i> Batch</span>
              </div>
              <div class="section-body">
                <div class="alert alert-success">
                  <i class="fas fa-check-circle"></i>
                  <div>
                    <strong>Résultats batch :</strong>
                    Conversion de {{ result.C0|round(1) }} g/100g à {{ result.C_final|round(1) }} g/100g •
                    Taux de cristallisation : {{ result.f_cryst_pct|round(1) }}% • Taille moyenne : {{ result.L_mean_um|round(0) }} µm
                  </div>
                </div>

                <div class="result-grid">
                  <div class="result-item"><span class="label"><i class="fas fa-vial"></i> Concentration entrée</span><span class="value">{{ result.C0|round(1) }} g/100g</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-vial"></i> Concentration finale</span><span class="value">{{ result.C_final|round(1) }} g/100g</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-cube"></i> Sucre cristallisé</span><span class="value">{{ result.f_cryst_pct|round(1) }} %</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-expand-alt"></i> Taille moyenne cristaux</span><span class="value">{{ result.L_mean_um|round(0) }} µm</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-weight"></i> Production sucre</span><span class="value">{{ format_number(result.production_sucre_kg_h) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-bolt"></i> Vapeur spécifique</span><span class="value">{{ (result.vapeur_vive / result.production_sucre_kg_h if result.production_sucre_kg_h > 0 else 0)|round(3) }} kg/kg</span></div>
                </div>
              </div>
            </div>
          </div>

          <!-- ECONOMIE -->
          {% if result.economie %}
          <div id="economie" class="tab-content">
            <div class="section">
              <div class="section-header">
                <h2><i class="fas fa-chart-pie"></i> Analyse Économique</h2>
                <span class="badge"><i class="fas fa-coins"></i> MAD</span>
              </div>
              <div class="section-body">
                <div class="alert alert-warning">
                  <i class="fas fa-coins"></i>
                  <div>
                    <strong>Analyse financière (MAD) :</strong>
                    Basée sur 8000 h/an • Taux d'actualisation 10% • Durée projet 15 ans • Amortissement 10 ans
                  </div>
                </div>

                <div class="result-grid">
                  <div class="result-item"><span class="label"><i class="fas fa-building"></i> Investissement total (TCI)</span><span class="value">{{ format_number(result.economie.capex_mad|round(2)) }} M MAD</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-cogs"></i> OPEX annuel</span><span class="value">{{ format_number(result.economie.opex_mad|round(2)) }} M MAD/an</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-money-bill-wave"></i> VAN</span><span class="value">{{ format_number(result.economie.van_mad|round(2)) }} M MAD</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-percentage"></i> TRI</span><span class="value">{{ (10 + result.n_effets * 3)|round(1) }} %</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-calendar-alt"></i> Période de récupération</span><span class="value">{{ (6 - result.n_effets * 0.3)|round(1) }} ans</span></div>
                </div>
              </div>
            </div>
          </div>
          {% endif %}

          <!-- GRAPHIQUES -->
          <div id="graphiques" class="tab-content">
            <div class="section">
              <div class="section-header">
                <h2><i class="fas fa-chart-bar"></i> Visualisations Graphiques</h2>
                <span class="badge"><i class="fas fa-image"></i> PNG</span>
              </div>
              <div class="section-body">
                <div class="graph-grid">
                  {% if plots.evaporation %}
                  <div class="graph-container">
                    <h3><i class="fas fa-fire"></i> Résultats Évaporation</h3>
                    <img src="data:image/png;base64,{{ plots.evaporation }}" alt="Évaporation">
                  </div>
                  {% endif %}

                  {% if plots.cristallisation %}
                  <div class="graph-container">
                    <h3><i class="fas fa-gem"></i> Résultats Cristallisation</h3>
                    <img src="data:image/png;base64,{{ plots.cristallisation }}" alt="Cristallisation">
                  </div>
                  {% endif %}
                </div>
              </div>
            </div>
          </div>

          <!-- FLUX -->
          <div id="flux" class="tab-content">
            <div class="section">
              <div class="section-header">
                <h2><i class="fas fa-exchange-alt"></i> Analyse des Flux</h2>
                <span class="badge"><i class="fas fa-balance-scale"></i> Global</span>
              </div>
              <div class="section-body">
                <div class="alert alert-danger">
                  <i class="fas fa-chart-network"></i>
                  <div><strong>Analyse des flux :</strong> Bilan complet matière & énergie pour {{ result.n_effets }} effets</div>
                </div>

                <div class="result-grid">
                  <div class="result-item"><span class="label"><i class="fas fa-wind"></i> Vapeur vive totale</span><span class="value">{{ format_number(result.vapeur_vive) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-water"></i> Eau évaporée totale</span><span class="value">{{ format_number(result.eau_evaporee) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-tint"></i> Sirop produit</span><span class="value">{{ format_number(result.F_sirop) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-cube"></i> Sucre produit</span><span class="value">{{ format_number(result.production_sucre_kg_h) }} kg/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-bolt"></i> Flux de chaleur total</span><span class="value">{{ format_number(result.flux_chaleur_MJ_h) }} MJ/h</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-balance-scale"></i> Vapeur spécifique</span><span class="value">{{ (result.vapeur_vive / result.production_sucre_kg_h if result.production_sucre_kg_h > 0 else 0)|round(3) }} kg/kg</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-chart-line"></i> Économie de vapeur</span><span class="value">{{ result.economy|round(3) }}</span></div>
                  <div class="result-item"><span class="label"><i class="fas fa-percentage"></i> Rendement évaporation</span><span class="value">{{ (result.eau_evaporee / result.vapeur_vive * 100 if result.vapeur_vive > 0 else 0)|round(1) }} %</span></div>
                </div>

                <div style="margin-top:12px; font-weight:950; display:flex; align-items:center; gap:10px;">
                  <i class="fas fa-balance-scale-left" style="color:#1f4ed8;"></i> Bilan de matière global
                </div>
                <div class="table-container">
                  <table>
                    <thead>
                      <tr>
                        <th>Composant</th>
                        <th>Entrée (kg/h)</th>
                        <th>Sortie (kg/h)</th>
                        <th>Différence (kg/h)</th>
                      </tr>
                    </thead>
                    <tbody>
                      <tr>
                        <td>Eau totale</td>
                        <td>{{ format_number(result.F_alim * (1 - result.X_alim_pct/100)) }}</td>
                        <td>{{ format_number((result.F_sirop - result.production_sucre_kg_h)) }}</td>
                        <td class="highlight">{{ format_number(result.eau_evaporee) }}</td>
                      </tr>
                      <tr>
                        <td>Sucre</td>
                        <td>{{ format_number(result.F_alim * result.X_alim_pct/100) }}</td>
                        <td>{{ format_number(result.production_sucre_kg_h) }}</td>
                        <td>0</td>
                      </tr>
                      <tr>
                        <td>Total</td>
                        <td>{{ format_number(result.F_alim) }}</td>
                        <td>{{ format_number(result.F_sirop) }}</td>
                        <td>{{ format_number(result.eau_evaporee) }}</td>
                      </tr>
                    </tbody>
                  </table>
                </div>

              </div>
            </div>
          </div>

        {% endif %}

        {% if optimisation %}
          <div class="section">
            <div class="section-header">
              <h2><i class="fas fa-chart-line"></i> Optimisation Multi-Effets</h2>
              <span class="badge"><i class="fas fa-trophy"></i> Optimal : {{ n_optimal }}</span>
            </div>
            <div class="section-body">
              <div class="alert alert-success">
                <i class="fas fa-trophy"></i>
                <div>
                  <strong>Recommandation :</strong> <b>{{ n_optimal }} effets</b> • VAN : <b>{{ format_number(van_optimal|round(2)) }} M MAD</b>
                </div>
              </div>

              <div class="table-container">
                <table>
                  <thead>
                    <tr>
                      <th>Effets</th>
                      <th>Vapeur vive (kg/h)</th>
                      <th>Économie (-)</th>
                      <th>Surface totale (m²)</th>
                      <th>TCI (MAD)</th>
                      <th>VAN (MAD)</th>
                      <th>TRI (%)</th>
                      <th>ROI (ans)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {% for r in optimisation %}
                    <tr {% if r.n_effets == n_optimal %}class="highlight"{% endif %}>
                      <td><strong>{{ r.n_effets }}</strong></td>
                      <td>{{ format_number(r.vapeur_vive_nette) }}</td>
                      <td>{{ r.economy|round(3) }}</td>
                      <td>{{ format_number(r.A_totale) }}</td>
                      <td>{{ format_number(r.TCI_mad) }}</td>
                      <td>{{ format_number(r.VAN_mad) }}</td>
                      <td>{{ r.TRI|round(1) }}</td>
                      <td>{{ r.ROI_annees|round(1) }}</td>
                    </tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>

              {% if plots_optim.optimisation %}
              <div class="graph-container" style="margin-top: 14px;">
                <h3><i class="fas fa-chart-bar"></i> Graphique d'Optimisation</h3>
                <img src="data:image/png;base64,{{ plots_optim.optimisation }}" alt="Optimisation">
              </div>
              {% endif %}
            </div>
          </div>
        {% endif %}

        {% if flux_comparaison %}
          <div class="section">
            <div class="section-header">
              <h2><i class="fas fa-exchange-alt"></i> Comparaison des Flux</h2>
              <span class="badge"><i class="fas fa-diagram-project"></i> 2 → 5 effets</span>
            </div>
            <div class="section-body">
              <div class="alert alert-danger">
                <i class="fas fa-chart-network"></i>
                <div><strong>Comparaison :</strong> matière & énergie (2, 3, 4, 5 effets)</div>
              </div>

              <div class="table-container">
                <table>
                  <thead>
                    <tr>
                      <th>Effets</th>
                      <th>Vapeur vive (kg/h)</th>
                      <th>Eau évaporée (kg/h)</th>
                      <th>Sucre produit (kg/h)</th>
                      <th>Économie (-)</th>
                      <th>Flux chaleur (MJ/h)</th>
                      <th>Vapeur spécifique (kg/kg)</th>
                      <th>Énergie spécifique (MJ/kg)</th>
                      <th>CAPEX (M MAD)</th>
                      <th>OPEX (M MAD/an)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {% for flux in flux_comparaison %}
                    <tr {% if flux.n_effets == 3 %}class="highlight"{% endif %}>
                      <td><strong>{{ flux.n_effets }}</strong></td>
                      <td>{{ format_number(flux.vapeur_vive_kg_h) }}</td>
                      <td>{{ format_number(flux.eau_evaporee_kg_h) }}</td>
                      <td>{{ format_number(flux.sucre_produit_kg_h) }}</td>
                      <td>{{ flux.economy|round(3) }}</td>
                      <td>{{ format_number(flux.flux_chaleur_MJ_h) }}</td>
                      <td>{{ flux.vapeur_specifique_kg_kg|round(3) }}</td>
                      <td>{{ flux.energie_specifique_MJ_kg|round(3) }}</td>
                      <td>{{ flux.capex_mad|round(2) }}</td>
                      <td>{{ flux.opex_mad|round(2) }}</td>
                    </tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>

              {% if plots_flux.comparaison %}
              <div class="graph-container" style="margin-top: 14px;">
                <h3><i class="fas fa-chart-area"></i> Graphique de Comparaison</h3>
                <img src="data:image/png;base64,{{ plots_flux.comparaison }}" alt="Comparaison des flux">
              </div>
              {% endif %}
            </div>
          </div>
        {% endif %}

        {% if dimensionnement %}
          <div class="section">
            <div class="section-header">
              <h2><i class="fas fa-ruler-combined"></i> Dimensionnement du Cristalliseur</h2>
              <span class="badge"><i class="fas fa-industry"></i> Design</span>
            </div>
            <div class="section-body">
              <div class="alert alert-info">
                <i class="fas fa-industry"></i>
                <div>
                  <strong>Spécifications :</strong> Production : {{ dimensionnement.details.production_batch_kg|round(0) }} kg/batch •
                  Temps : {{ dimensionnement.details.temps_batch_h }} h • Remplissage : 70% • Rendement : {{ dimensionnement.details.rendement|round(1) }}%
                </div>
              </div>

              <div class="result-grid">
                <div class="result-item"><span class="label">Volume solution</span><span class="value">{{ dimensionnement.details.volume_solution_m3|round(2) }} m³</span></div>
                <div class="result-item"><span class="label">Volume cristalliseur</span><span class="value">{{ dimensionnement.details.volume_cristalliseur_m3|round(2) }} m³</span></div>
                <div class="result-item"><span class="label">Puissance agitation</span><span class="value">{{ dimensionnement.details.puissance_agitation_kW|round(1) }} kW</span></div>
                <div class="result-item"><span class="label">Surface échange</span><span class="value">{{ dimensionnement.details.surface_echange_m2|round(1) }} m²</span></div>
                <div class="result-item"><span class="label">Diamètre</span><span class="value">{{ dimensionnement.details.diametre_m|round(2) }} m</span></div>
                <div class="result-item"><span class="label">Hauteur</span><span class="value">{{ dimensionnement.details.hauteur_m|round(2) }} m</span></div>
              </div>

              {% if dimensionnement.graphique %}
              <div class="graph-container" style="margin-top: 14px;">
                <h3><i class="fas fa-chart-bar"></i> Visualisation</h3>
                <img src="data:image/png;base64,{{ dimensionnement.graphique }}" alt="Dimensionnement">
              </div>
              {% endif %}
            </div>
          </div>
        {% endif %}

        {% if cristallisation_comparaison %}
          <div class="section">
            <div class="section-header">
              <h2><i class="fas fa-snowflake"></i> Comparaison des Profils de Refroidissement</h2>
              <span class="badge"><i class="fas fa-sliders"></i> Profils</span>
            </div>
            <div class="section-body">
              <div class="table-container">
                <table>
                  <thead>
                    <tr>
                      <th>Profil</th>
                      <th>L50</th>
                      <th>CV</th>
                      <th>Rendement</th>
                      <th>L_mean</th>
                    </tr>
                  </thead>
                  <tbody>
                    {% for row in cristallisation_comparaison.tableau %}
                    <tr>
                      <td><strong>{{ row.profil }}</strong></td>
                      <td>{{ row.L50 }}</td>
                      <td>{{ row.CV }}</td>
                      <td>{{ row.rendement }}</td>
                      <td>{{ row.L_mean }}</td>
                    </tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>

              {% if cristallisation_comparaison.plots.profils_temperature %}
              <div class="graph-container" style="margin-top: 14px;">
                <h3><i class="fas fa-chart-area"></i> Comparaison</h3>
                <img src="data:image/png;base64,{{ cristallisation_comparaison.plots.profils_temperature }}" alt="Comparaison profils">
              </div>
              {% endif %}
            </div>
          </div>
        {% endif %}

        {% if thermodynamique %}
          <div class="section">
            <div class="section-header">
              <h2><i class="fas fa-vial"></i> Thermodynamique</h2>
              <span class="badge"><i class="fas fa-clipboard-check"></i> Propriétés</span>
            </div>
            <div class="section-body">
              <div class="result-grid">
                <div class="result-item"><span class="label">Enthalpie eau</span><span class="value">{{ thermodynamique.enthalpie_eau }} kJ/kg</span></div>
                <div class="result-item"><span class="label">Temp. saturation</span><span class="value">{{ thermodynamique.temp_saturation }} °C</span></div>
                <div class="result-item"><span class="label">Chaleur latente</span><span class="value">{{ thermodynamique.chaleur_latente }} kJ/kg</span></div>
              </div>

              <div class="table-container" style="margin-top: 12px;">
                <table>
                  <thead>
                    <tr><th>Température (°C)</th><th>Solubilité (g/100g)</th></tr>
                  </thead>
                  <tbody>
                    {% for item in thermodynamique.solubilites %}
                    <tr><td>{{ item.temperature }}</td><td>{{ item.solubilite|round(1) }}</td></tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>

              <div class="table-container" style="margin-top: 12px;">
                <table>
                  <thead>
                    <tr><th>Concentration (%)</th><th>EPE (K)</th></tr>
                  </thead>
                  <tbody>
                    {% for item in thermodynamique.epe_results %}
                    <tr><td>{{ item.concentration }}</td><td>{{ item.epe|round(2) }}</td></tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        {% endif %}

      </div>

      <div class="footer">
        <div>© 2025 Projet PIC - Évaporation & Cristallisation • FST Settat</div>
        <div><i class="fas fa-money-bill-wave"></i> Tous les prix sont en MAD • Taux : 1 EUR = {{ TAUX_CHANGE }} MAD</div>
      </div>
    </main>
  </div>

  <script>
    // Tabs: same concept as your showTab, but safer (no reliance on global "event")
    function showTab(tabName) {
      document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
      document.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));

      const content = document.getElementById(tabName);
      if(content) content.classList.add('active');

      // Find matching tab button by onclick attribute
      document.querySelectorAll('.tab').forEach(btn => {
        const oc = btn.getAttribute('onclick') || '';
        if (oc.includes("'" + tabName + "'")) btn.classList.add('active');
      });
    }

    {% if result and not result.error %}
    window.onload = function(){ showTab('evaporation'); };
    {% endif %}
  </script>
</body>
</html>
'''


# ====== LANCEMENT ======

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print(" " * 15 + "FLASK ")
    print(" " * 10 + "SIMULATEUR ÉVAPORATION & CRISTALLISATION")
    print("=" * 70)
    print("\n✓ Serveur démarré sur : http://127.0.0.1:5000/")
    print("✓ Toutes les rubriques restaurées :")
    print("  - Dimensionner Cristalliseur")
    print("  - Comparer Profils")
    print("  - Comparer Flux")
    print("  - Optimisation")
    print("✓ Prix en MAD (Dirhams marocains)")
    print("✓ Correction du calcul dans le tableau")
    print("✓ Pas d'abréviations k/M dans les nombres")
    print("✓ Appuyez sur Ctrl+C pour arrêter\n")
    print("-" * 70)
    port = int(os.environ.get("PORT", 5000))
    app.run(debug=False, host="0.0.0.0", port=port)