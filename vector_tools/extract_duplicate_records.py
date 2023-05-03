import geopandas as gpd
import pandas as pd

# names
# inputs
troncon_hydro_name = '4_troncon_hydrographique_cours_d_eau_corr_conn_inv'
# outputs
duplicate_name = 'troncon_hydrographique_corr_duplicate'

# Sélectionner la colonne avec les valeurs à vérifier pour la duplication
colonne_a_verifier = 'cleabs'

# pgkg file exutoire
referentiel_hydro_gpkg = 'C:/Users/lmanie01/Documents/Projets/Mapdo/Data/fct/referentiel_hydrographique.gpkg'

# Charger le fichier shapefile
gdf = gpd.read_file(referentiel_hydro_gpkg, layer=troncon_hydro_name)

# Vérifier les valeurs dupliquées et les stocker dans un nouveau DataFrame
duplications = gdf[gdf.duplicated(colonne_a_verifier, keep=False)]

# Si des duplications ont été trouvées, créer une couche avec ces entités
if not duplications.empty:
    # Nom de la nouvelle couche geopackage
    nom_couche = duplicate_name+'.gpkg'
    
    # Créer une couche geopackage avec les entités dupliquées
    duplications.to_file(referentiel_hydro_gpkg, layer = duplicate_name, driver='GPKG')
    
    print(f'Couche créée avec succès : {duplicate_name}')
else:
    print('Aucune duplication trouvée.')
