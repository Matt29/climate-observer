# Fièvre océan

Page unique, sans dépendance, qui montre chaque jour l'écart à la normale des océans et des continents
(globe tournant + vues plates + écart local au clic) et un tracker El Niño (Niño 3.4 du jour, 12 mois,
Hovmöller équatorial, ONI depuis 1950).

Prototype validé en chat le 17/09/2026 sur les données du 17/08 au 15/09/2026.

## Structure

    pipeline/build_data.py   télécharge, calcule, écrit data/build.json (~5 Mo)
    build_page.py            injecte les données dans web/template.html -> dist/index.html
    web/template.html        toute l'UI (HTML/CSS/JS), zéro lib externe sauf Google Fonts
    web/coast.min.json       trait de côte Natural Earth 110m compacté (polylignes [lon,lat])
    .github/workflows/       cron quotidien + déploiement GitHub Pages (à adapter Vercel/R2)

## Lancer

    pip install -r requirements.txt
    python pipeline/build_data.py        # long la première fois (cache ~2,5 Go), ~1 min ensuite
    python build_page.py
    open dist/index.html

## Ce que la page embarque

- `META` : dates, grille 1°, séries globales (g60 = océan 60S-60N, gmean, sst60, frac_hot, lmean = terres)
- `B64` : cube int8 (t, 180, 360) anomalie mer ×10, -128 = terre/absent, base64
- `LB64` : idem pour la terre (CPC), -128 = océan ou pas de station
- `COAST` : polylignes
- `ENSO` : `hov` (dates, lons 120E-290E par 2°, rows, n34), `oni` [[saison, année, valeur]], `wk` [[date, N1+2, N3, N3.4, N4]]

Le JS : décodage, dilatation de 4 passes des NaN côtiers pour un rendu lisse, LUT couleur −6/+6,
rendu orthographique pixel par pixel (inverse orthographique, ombrage, terminateur, halo, étoiles,
trait de côte projeté), vues plates avec repli au 180e, lookup du point le plus proche avec donnée,
sparkline, onglets, graphes ENSO à axes temporels proportionnels (mélange 5 jours / quotidien).

## Sources

| Couche | Source | Latence | Base clim. |
|---|---|---|---|
| Mer | NOAA OISST v2.1 daily, fichiers NCEI `avhrr/YYYYMM/` | J-1 prélim., J-14 final | 1991-2020 (`sst` − ltm PSL, recalculée) |
| Terre | NOAA CPC Global Daily Temp Tmax/Tmin (PSL) | J-1 à J-2 | 1991-2020 (ltm PSL) |
| ENSO hebdo | CPC `wksst9120.for` | lundi | 1991-2020 |
| ONI | CPC `oni.ascii.txt` | ~10 du mois | glissante 30 ans |

## Pièges connus (déjà rencontrés)

1. **NCSS PSL corrompt les fichiers Tmin** (une ligne de latitude sur deux à 0), et PSL coupe les transferts
   tous les ~10 Mo, ce qu'une réponse NCSS générée à la volée ne sait pas reprendre (tmax bloqué à 12 Mo le 17/09).
   Le script télécharge donc `tmax|tmin.YYYY.nc` et `tmax|tmin.day.ltm.1991-2020.nc` en entier via `fileServer`.
   `get()` écrit dans `.part`, reprend avec `curl -C -` tant que le fichier grossit, et ne renomme qu'au succès
   de curl : un fichier en cache est toujours complet.
2. **Bases climatologiques différentes** : le champ `anom` des fichiers OISST est sur 1971-2000. Corrigé : l'anomalie
   mer est recalculée `sst − sst.day.mean.ltm.1991-2020.nc` (PSL, 0,25°, 1,4 Go en cache), même base que la terre,
   les indices hebdo et l'ONI. Le 29/02 prend la clim du 28/02 (clim de 365 j). Effet mesuré le 15/09/2026 : −0,09 °C sur
   60S-60N (−0,3 à −0,6 °C vers 30-70°N, ~0 dans l'océan Austral), −0,02 °C sur Niño 3.4. Premier téléchargement long (~140 reprises, cf. 1) ;
   un fichier tronqué donne `NetCDF: HDF error` à l'ouverture.
3. Le calendrier des fichiers ltm est `gregorian` année 1 : passer par `cftime.num2date`, pas par datetime.
4. CPC est un produit stations : trous (Antarctique, Groenland, Sahara, Amazonie) et quelques cellules
   aberrantes (jusqu'à −50 °C). Pour une couverture complète : ERA5 (CDS, J-5) ou analyse GFS 00Z.
5. Un artifact Claude publié ne peut pas fetcher : la page doit être régénérée et hébergée (Pages/Vercel/R2).

## Roadmap

- Prévision : ECMWF Open Data (T2m, 10 j) + plumes ENSO CPC/IRI.
- Servir le 0,25° en tuiles (COG/PMTiles) plutôt qu'embarquer le 1° ; petit Zarr pour le clic.
- Rendu vidéo quotidien du globe (matplotlib/cartopy ou Playwright + ffmpeg) pour TikTok.
- SOI, sous-surface TAO/TRITON, vagues de chaleur marines (catégories Hobday).
