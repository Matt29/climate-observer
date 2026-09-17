# Climate Observer

A daily, dependency-free web page showing how far today's ocean and land temperatures are from normal
(spinning globe, flat maps, local anomaly on click), with a 9-day forecast, an El Niño tracker
(today's Niño 3.4, last 12 months, equatorial Hovmöller, ONI since 1950) and ocean watch indicators
(sea ice, coral bleaching alerts).

Live page: https://oceandataconsulting.fr/climate-observer, rebuilt every morning by GitHub Actions.
Built by [OceanData Consulting](https://oceandataconsulting.fr).

## Layout

    pipeline/build_data.py   downloads, computes, writes data/build.json (~7 MB)
    pipeline/surveillance.py sea ice (NSIDC), corals (CRW), ENSO probabilities (CPC); `python pipeline/surveillance.py` runs a self-check
    build_page.py            injects the data into web/template.html -> dist/index.html
    web/template.html        the whole UI (HTML/CSS/JS), no external library except Google Fonts
    web/coast.min.json       compacted Natural Earth 110m coastline ([lon, lat] polylines)
    .github/workflows/       daily cron + GitHub Pages deployment (data/build.json is not versioned, CI regenerates it)

## Run locally

    uv venv -p 3.12 && uv pip install -r requirements.txt
    .venv/bin/python pipeline/build_data.py   # long the first time (~2.5 GB cache), a few minutes afterwards
    .venv/bin/python build_page.py
    open dist/index.html

The ocean forecast needs a free [Copernicus Marine](https://marine.copernicus.eu) account: either run
`copernicusmarine login` once, or set `COPERNICUSMARINE_SERVICE_USERNAME` / `COPERNICUSMARINE_SERVICE_PASSWORD`.
Without it the forecast is skipped and the page is still built.

## Embedded payload

- `META`: dates (observed then forecast), `nobs` (number of observed days), `fc_land`, 1° grid,
  1° series (`g60` = ocean 60S-60N, `gmean`, `frac_hot`, `lmean` = land)
- `B64`: int8 cube (t, 180, 360) of sea anomaly ×10, -128 = land/missing, base64; frames ≥ `nobs` are the forecast
- `LB64`: same for land (CPC), -128 = ocean or no station
- `COAST`: polylines
- `ENSO`: `hov` (dates, lons 120E-290E every 2°, rows, n34), `oni` [[season, year, value]], `wk` [[date, N1+2, N3, N3.4, N4]], `probs` [[season, La Niña, neutral, El Niño]]
- `SURV`: `ice` (NSIDC north/south: extent, normal, rank_low, nyears, same_date{year}), `ice_maps` (0.5° uint8 polar caps,
  `now` and `ref` per reference year), `coral` (alert1/alert2/warning shares). A missing block means that source failed.

The JS decodes the cubes, dilates coastal NaNs over 4 passes for smooth rendering, applies a −6/+6 colour LUT,
renders the orthographic globe pixel by pixel (inverse projection, shading, terminator, halo, stars, projected
coastline), draws flat maps with dateline wrap, looks up the nearest cell with data, and plots the ENSO charts
on proportional time axes (5-day / daily mix).

## Sources

| Layer | Source | Latency | Baseline |
|---|---|---|---|
| Sea | NOAA OISST v2.1 daily, NCEI files `avhrr/YYYYMM/` | D-1 preliminary, D-14 final | 1991-2020 (`sst` − PSL ltm, recomputed) |
| Land | NOAA CPC Global Daily Temperature Tmax/Tmin (PSL) | D-1 to D-2 | 1991-2020 (PSL ltm) |
| Weekly ENSO | CPC `wksst9120.for` | Mondays | 1991-2020 |
| ONI | CPC `oni.ascii.txt` | ~10th of the month | rolling 30 years |
| ENSO probabilities | CPC `enso/roni/probabilities/` (HTML table) | ~10th of the month | RONI |
| Sea forecast | Copernicus Marine GLO12 `thetao` 0.49 m (account required) | D0 | anomaly anchored on OISST 1991-2020 |
| Land forecast | ECMWF IFS open data `2t`, 00 UTC run of the last observed day | D0, ~4-day retention | anomaly anchored on CPC 1991-2020 |
| Sea ice | NSIDC Sea Ice Index v4 (figures) + OISST `ice` field (map) | D-1 | normal recomputed on 1991-2020 (NSIDC ships 1981-2010) |
| Corals | NOAA Coral Reef Watch, ERDDAP `NOAA_DHW` `CRW_BAA_7D_MAX` (0.5°) | D-1 | alert levels 0-4 |

Each dataset remains subject to its provider's terms of use.

## Known pitfalls

1. **PSL NCSS corrupts Tmin files** (every other latitude row set to 0), and PSL drops transfers every ~10 MB,
   which an on-the-fly NCSS response cannot resume. The pipeline therefore downloads whole `tmax|tmin.YYYY.nc` and
   `tmax|tmin.day.ltm.1991-2020.nc` files through `fileServer`. `get()` writes to `.part`, resumes with `curl -C -`
   as long as the file grows, and renames only when curl succeeds: a cached file is always complete.
2. **Mismatched climatologies**: the `anom` field of OISST files uses 1971-2000. The sea anomaly is recomputed as
   `sst − sst.day.mean.ltm.1991-2020.nc` (PSL, 0.25°, 1.4 GB), the same baseline as land, the weekly indices and ONI.
   Feb 29 uses the Feb 28 climatology (365-day climatology). Measured effect on 2026-09-15: −0.09 °C over 60S-60N
   (−0.3 to −0.6 °C around 30-70°N, ~0 in the Southern Ocean), −0.02 °C on Niño 3.4. A truncated file raises
   `NetCDF: HDF error` on open.
3. The ltm files use a `gregorian` calendar starting in year 1: go through `cftime.num2date`, not datetime.
4. CPC is a station product: gaps (Antarctica, Greenland, Sahara, Amazon) and a few outlier cells (down to −50 °C).
   For full coverage, use ERA5 (CDS, D-5) or the GFS 00Z analysis.
5. **Forecasts use the delta method**: `anomaly(d) = observed anomaly(D0) + [model(d) − model(D0)] − [clim(d) − clim(D0)]`.
   Subtracting the OISST/CPC climatology from the raw model field would create a step at the boundary (model bias).
   Checked on 2026-09-15: no jump (median day-to-day change 0.1 °C at sea, ~1 °C on land, same on both sides).
   IFS open data stops at 240 h, so the forecast is limited to the days covered by both sea and land (9 days).
   Series (`g60`, `frac_hot`, `lmean`) are computed at 1° over the whole cube for the same reason.
   Every optional source (forecasts, NSIDC, CRW, probabilities) fails gracefully: its block is omitted, the build goes on.
6. **Sea ice**: before 1987 NSIDC only has one value every 2 days, so reference years use the neighbouring day (±1).
   2020 (2nd-lowest Arctic year) would show "more ice than in 2020", hence 1982 as the default, and the map shows
   both lost and gained ice. **Corals**: `NOAA_DHW` has no reef mask, so the page shows the share of the tropical
   ocean under alert, not "of reefs".
7. **CI cache**: only the three 1991-2020 climatologies (~1.8 GB, immutable) are cached, under the fixed key
   `ltm-1991-2020-v1` (saved once, only when the job succeeds); the rest (~250 MB) is downloaded again every day.
   Changing the baseline means changing the key. Required secrets: `COPERNICUSMARINE_SERVICE_USERNAME` /
   `COPERNICUSMARINE_SERVICE_PASSWORD`. GitHub disables a scheduled workflow after 60 days without repository
   activity: re-enable it in the Actions tab if the page stops updating.

## Roadmap

- Serve the 0.25° field as tiles (COG/PMTiles) instead of embedding the 1° grid; small Zarr for click lookups.
- Daily video render of the globe.
- SOI, TAO/TRITON subsurface, marine heatwaves (Hobday categories, 90th percentile to compute).

## License

Code released under the [MIT License](LICENSE). The data belongs to its providers (NOAA, NSIDC, Copernicus Marine, ECMWF) and follows their own terms.
