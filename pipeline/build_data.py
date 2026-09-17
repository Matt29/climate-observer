#!/usr/bin/env python3
"""
Climate Observer, pipeline quotidien.

Produit data/build.json (tout ce que la page embarque) à partir de :
  - NOAA OISST v2.1 daily (mer, 0,25°) : 30 derniers jours + 12 mois à pas de 5 jours (ENSO)
  - NOAA CPC Global Daily Temperature (terre, 0,5°, Tmax/Tmin) + climatologie 1991-2020 (PSL)
  - NOAA CPC : indices Niño hebdo (wksst9120.for) et ONI (oni.ascii.txt)

Usage :
  python pipeline/build_data.py [--end YYYY-MM-DD] [--days 30] [--cache cache/]
puis
  python build_page.py

Notes importantes (voir README) :
  - L'anomalie OISST livrée par NOAA est sur base 1971-2000 : on la recalcule sur 1991-2020
    (sst quotidien - sst.day.mean.ltm.1991-2020.nc de PSL, ~1,4 Go en cache), comme la terre et les indices CPC.
  - Le NetCDF Subset Service de PSL renvoie des fichiers Tmin corrompus (une ligne sur deux à 0) et ne sait
    pas reprendre un transfert coupé. On télécharge donc les fichiers annuels complets pour tmax et tmin.
"""
import argparse, base64, datetime as dt, functools, json, os, re, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import xarray as xr
import cftime

import surveillance

OISST_BASE = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/"
PSL_FILES = "https://psl.noaa.gov/thredds/fileServer/Datasets/cpc_global_temp/"
GLO12 = "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m"
OISST_LTM = "https://psl.noaa.gov/thredds/fileServer/Datasets/noaa.oisst.v2.highres/sst.day.mean.ltm.1991-2020.nc"
CPC_WK = "https://www.cpc.ncep.noaa.gov/data/indices/wksst9120.for"
CPC_ONI = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt"

LAND = -128  # sentinel int8 = pas de donnée


def get(url, path, stall=3):
    """Téléchargement avec reprise (curl -C -) dans path.part, renommé seulement si curl a tout reçu.
    PSL coupe les transferts tous les ~10 Mo (fileServer et NCSS) : on reprend tant que ça avance,
    abandon après `stall` essais sans progrès. Un fichier non vide n'est pas un fichier complet."""
    if os.path.exists(path):
        return path
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    part, fails, size = path + ".part", 0, -1
    while fails < stall:
        rc = os.system(f'curl -s -f -C - -o "{part}" "{url}"') >> 8
        if rc == 0 and os.path.getsize(part) > 0:
            os.replace(part, path)
            return path
        if rc == 33 and os.path.exists(part):   # reprise refusée : on repart de zéro
            os.remove(part)
        new = os.path.getsize(part) if os.path.exists(part) else 0
        fails, size = (0, new) if new > size else (fails + 1, size)
    raise RuntimeError(f"échec téléchargement {url} (curl rc={rc}, {size} octets)")


# ---------------------------------------------------------------- OISST
def oisst_listing(month):
    """Liste des fichiers disponibles pour un mois (final ou _preliminary)."""
    html = urllib.request.urlopen(OISST_BASE + month + "/").read().decode()
    files = sorted(set(re.findall(r"oisst-avhrr-v02r01\.[0-9_a-z]*\.nc", html)))
    return {re.search(r"\.(\d{8})", f).group(1): f for f in files}


def oisst_files_between(d0, d1):
    """Dict date -> nom de fichier, tous les jours dispo entre d0 et d1 inclus."""
    out = {}
    m = dt.date(d0.year, d0.month, 1)
    while m <= d1:
        out.update({k: v for k, v in oisst_listing(m.strftime("%Y%m")).items()
                    if d0.strftime("%Y%m%d") <= k <= d1.strftime("%Y%m%d")})
        m = (m.replace(day=28) + dt.timedelta(days=4)).replace(day=1)
    return out


def fetch_oisst(dates_files, cache):
    os.makedirs(cache, exist_ok=True)
    def one(item):
        k, f = item
        return k, get(OISST_BASE + k[:6] + "/" + f, os.path.join(cache, f))
    with ThreadPoolExecutor(6) as ex:
        return dict(ex.map(one, sorted(dates_files.items())))


@functools.cache
def oisst_ltm(path):
    """Climatologie OISST 1991-2020 (PSL) + index (mois, jour) -> pas de temps.
    Calendrier gregorian année 1 : cftime, pas datetime (voir README)."""
    ds = xr.open_dataset(path, decode_times=False)
    days = cftime.num2date(ds.time.values, ds.time.attrs["units"], "gregorian")
    return ds.sst, {(d.month, d.day): k for k, d in enumerate(days)}


def open_anom(path, ltm_path):
    """Anomalie SST sur base 1991-2020 (le champ `anom` NOAA est sur 1971-2000)."""
    ds = xr.open_dataset(path)
    sst = ds.sst.isel(time=0, zlev=0).load()
    d = dt.date.fromisoformat(str(ds.time.values[0])[:10])
    clim, idx = oisst_ltm(ltm_path)
    k = idx.get((d.month, d.day), idx.get((2, 28)))   # 29/02 absent d'une clim 365 j
    c = clim.isel(time=k)
    assert np.allclose(c.lat, sst.lat) and np.allclose(c.lon, sst.lon), "grille ltm != OISST"
    return sst - c.values, sst, d.isoformat()


def to_180(da):
    return da.assign_coords(lon=((da.lon + 180) % 360) - 180).sortby("lon")


def encode_cube(arr, scale=10, clip=120):
    """float (t, lat, lon) -> int8 base64, NaN -> LAND."""
    q = np.clip(np.round(arr * scale), -clip, clip)
    q = np.where(np.isnan(arr), LAND, q).astype(np.int8)
    return base64.b64encode(q.tobytes()).decode()


def build_ocean(paths_by_date, ltm_path):
    """Cube 1° (xarray) des anomalies des 30 derniers jours."""
    anoms, dates = [], []
    for k in sorted(paths_by_date):
        a, _, d = open_anom(paths_by_date[k], ltm_path)
        anoms.append(to_180(a)); dates.append(d)
    return dates, xr.concat(anoms, "time").coarsen(lat=4, lon=4, boundary="exact").mean()


def ocean_forecast(obs1, t0, ltm_path, cache, days):
    """Anomalie mer prévue à 1° par la méthode des écarts : anomalie observée du jour t0
    + évolution prévue par GLO12 (Copernicus Marine) - évolution de la climatologie.
    Le biais du modèle par rapport à OISST s'annule : pas de marche entre observé et prévu."""
    path = os.path.join(cache, f"glo12_thetao_1deg_{t0:%Y%m%d}.nc")
    if not os.path.exists(path):
        import copernicusmarine   # identifiants : ~/.copernicusmarine ou COPERNICUSMARINE_SERVICE_USERNAME/PASSWORD
        os.makedirs(cache, exist_ok=True)
        ds = copernicusmarine.open_dataset(dataset_id=GLO12, variables=["thetao"], minimum_depth=0, maximum_depth=1,
                                           start_datetime=t0.isoformat(), end_datetime=(t0 + dt.timedelta(days=days)).isoformat())
        th = ds.thetao.isel(depth=0).coarsen(latitude=12, longitude=12, boundary="trim").mean()   # 1/12° -> 1°, blocs à bords entiers
        th = th.assign_coords(latitude=np.floor(th.latitude) + 0.5, longitude=np.floor(th.longitude) + 0.5)
        th.rename(latitude="lat", longitude="lon").load().to_netcdf(path + ".part")
        os.replace(path + ".part", path)
    th = xr.open_dataarray(path).reindex(lat=obs1.lat, lon=obs1.lon, method="nearest", tolerance=0.01)   # GLO12 s'arrête à 80°S
    tdates = [str(x)[:10] for x in th.time.values]
    assert tdates[0] == t0.isoformat(), f"GLO12 sans le jour d'ancrage {t0}"
    clim, idx = oisst_ltm(ltm_path)
    def clim1(d):
        c = clim.isel(time=idx.get((d.month, d.day), idx[(2, 28)])).load()
        return to_180(c).coarsen(lat=4, lon=4, boundary="exact").mean().values
    a0, m0, c0 = obs1.isel(time=-1).values, th.isel(time=0).values, clim1(t0)
    out = []
    for k in range(1, len(tdates)):
        d = dt.date.fromisoformat(tdates[k])
        out.append(a0 + (th.isel(time=k).values - m0) - (clim1(d) - c0))
    return np.array(out)


# ---------------------------------------------------------------- terre (CPC)
def build_land(d0, d1, cache, days):
    """Anomalie 0,5° des jours d0..d1 et climatologie (Tmax+Tmin)/2 des jours d1..d1+days (pour la prévision)."""
    year = d1.year
    n = (d1 - d0).days + 1

    def field(var):
        """Valeurs du jour et climatologie 1991-2020 alignée sur le jour de l'année.
        Fichiers complets : NCSS corrompt tmin et ne sait pas reprendre un transfert coupé (voir README)."""
        path = os.path.join(cache, f"{var}.{year}.nc")
        for _ in range(2):   # le fichier annuel grossit chaque jour : copie en cache trop courte -> on la reprend
            with xr.open_dataset(get(PSL_FILES + f"{var}.{year}.nc", path)) as ds:
                obs = ds[var].sel(time=slice(str(d0), str(d1))).load()
            if obs.sizes["time"] == n:
                break
            os.remove(path)
        assert obs.sizes["time"] == n, f"jours manquants dans {var}"
        ltm = xr.open_dataset(get(PSL_FILES + f"{var}.day.ltm.1991-2020.nc",
                                  os.path.join(cache, f"{var}.day.ltm.1991-2020.nc")), decode_times=False)
        days_ = cftime.num2date(ltm.time.values, ltm.time.attrs["units"], "gregorian")   # année 1 : cftime
        i0 = next(k for k, d in enumerate(days_) if (d.month, d.day) == (d0.month, d0.day))
        return obs, ltm[var].isel(time=[(i0 + k) % len(days_) for k in range(n + days)]).values

    tx, lxv = field("tmax")
    tn, ln = field("tmin")
    clim = (lxv + ln) / 2
    wrap = lambda arr: to_180(xr.DataArray(arr, dims=("time", "lat", "lon"),
                                           coords=dict(lat=tx.lat.values, lon=tx.lon.values))).sortby("lat")
    return wrap((tx.values + tn.values) / 2 - clim[:n]), wrap(clim[n - 1:])


def land_forecast(anom05, clim05, t0, cache, days):
    """Anomalie terre prévue à 1°, même méthode des écarts avec ECMWF IFS open data (2t, run 00 UTC du jour t0).
    Moyenne des échéances 00/06/12/18 UTC par jour, comparable à (Tmax+Tmin)/2 à la méthode près, ce qui s'annule dans l'écart."""
    path = os.path.join(cache, f"ifs_2t_{t0:%Y%m%d}00.grib2")
    if not os.path.exists(path):
        from ecmwf.opendata import Client
        os.makedirs(cache, exist_ok=True)
        Client(source="ecmwf").retrieve(date=t0.isoformat(), time=0, type="fc", param="2t",
                                        step=list(range(0, 24 * days, 6)), target=path + ".part")
        os.replace(path + ".part", path)
    t2 = xr.open_dataset(path, engine="cfgrib", indexpath="").t2m.sortby("latitude")
    # la grille 0,25° d'IFS contient exactement les centres 0,5° de CPC : sélection, pas d'interpolation (ni scipy)
    daily = t2.coarsen(step=4, boundary="trim").mean().sel(latitude=anom05.lat.values, longitude=anom05.lon.values, method="nearest").values
    a0 = anom05.isel(time=-1).values
    fc = np.array([a0 + (daily[k] - daily[0]) - (clim05.values[k] - clim05.values[0]) for k in range(1, len(daily))])
    return xr.DataArray(fc, dims=("time", "lat", "lon"), coords=dict(lat=anom05.lat, lon=anom05.lon)) \
             .coarsen(lat=2, lon=2, boundary="exact").mean().values


# ---------------------------------------------------------------- glace de mer (carte)
ICE_YEARS = (1982, 2000, 2012, 2020)


def ice_grid(path):
    """Concentration de glace OISST à 0,5°, uint8 0-100, 255 = terre. Lignes par latitude croissante, lon -180..180."""
    ds = xr.open_dataset(path)
    ice, sst = to_180(ds.ice.isel(time=0, zlev=0)), to_180(ds.sst.isel(time=0, zlev=0))
    c = ice.fillna(0).where(sst.notnull()).coarsen(lat=2, lon=2, boundary="exact").mean().values
    return np.where(np.isnan(c), 255, np.round(c * 100)).astype(np.uint8)


def ice_maps(today_path, t0, cache):
    """Calottes nord (>= 40°N) et sud (<= 50°S) du jour et de la même date les années repères."""
    cut = lambda g: dict(north=base64.b64encode(g[260:].tobytes()).decode(), south=base64.b64encode(g[:80].tobytes()).decode())
    ref = {}
    for y in ICE_YEARS:
        d = dt.date(y, t0.month, min(t0.day, 28) if t0.month == 2 else t0.day)
        ref[y] = cut(ice_grid(next(iter(fetch_oisst(oisst_files_between(d, d), os.path.join(cache, "oisst_ref")).values()))))
    return dict(step=0.5, nlon=720, lon0=-179.75, north_lat0=40.25, north_nlat=100, south_lat0=-89.75, south_nlat=80,
                now=cut(ice_grid(today_path)), ref=ref)


def series(ocean1, land1):
    """Moyennes quotidiennes à 1° sur tout le cube (observé + prévu), mêmes calculs des deux côtés de la frontière."""
    lat = ocean1.lat
    w = np.cos(np.deg2rad(lat)); band = slice(-60, 60)
    r = lambda da, nd: [None if np.isnan(x) else round(float(x), nd) for x in da.values]
    return dict(
        gmean=r(ocean1.weighted(w).mean(("lat", "lon")), 3),
        g60=r(ocean1.sel(lat=band).weighted(w.sel(lat=band)).mean(("lat", "lon")), 3),
        frac_hot=r((ocean1 > 1).where(ocean1.notnull()).weighted(w).mean(("lat", "lon")), 3),
        lmean=r(land1.weighted(w).mean(("lat", "lon")), 2),
    )


def soft(name, fn):
    """Sources optionnelles : un échec ne casse pas le build quotidien, le bloc est omis de la page."""
    try:
        return fn()
    except Exception as e:
        print(f"[dégradé] {name} : {type(e).__name__}: {e}", file=sys.stderr)
        return None


# ---------------------------------------------------------------- ENSO
def build_enso(paths_by_date, ltm_path):
    rows, n34, dates = [], [], []
    lons = None
    for k in sorted(paths_by_date):
        a, _, d = open_anom(paths_by_date[k], ltm_path)           # lon 0..360, pratique pour le Pacifique
        eq = a.sel(lat=slice(-5, 5)).mean("lat")
        pac = eq.sel(lon=slice(120, 290)).coarsen(lon=8, boundary="trim").mean()   # bins de 2°
        rows.append([None if np.isnan(x) else round(float(x), 2) for x in pac.values])
        n34.append(round(float(a.sel(lat=slice(-5, 5), lon=slice(190, 240)).mean()), 2))
        dates.append(d); lons = [round(float(x), 1) for x in pac.lon.values]
    oni = []
    for l in urllib.request.urlopen(CPC_ONI).read().decode().splitlines()[1:]:
        p = l.split()
        if len(p) == 4:
            oni.append([p[0], int(p[1]), float(p[3])])
    wk = []
    rx = re.compile(r"\s*(\d{2}[A-Z]{3}\d{4})\s+([\d.]+)\s*(-?[\d.]+)\s+([\d.]+)\s*(-?[\d.]+)\s+([\d.]+)\s*(-?[\d.]+)\s+([\d.]+)\s*(-?[\d.]+)")
    for l in urllib.request.urlopen(CPC_WK).read().decode().splitlines()[4:]:
        m = rx.match(l)
        if m:
            wk.append([m.group(1)] + [float(m.group(i)) for i in (3, 5, 7, 9)])
    return dict(hov=dict(dates=dates, lons=lons, rows=rows, n34=n34), oni=oni, wk=wk)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--end", help="dernier jour (défaut : dernier fichier OISST dispo)")
    ap.add_argument("--days", type=int, default=30)
    ap.add_argument("--cache", default="cache")
    ap.add_argument("--out", default="data/build.json")
    ap.add_argument("--fc-days", type=int, default=10)
    args = ap.parse_args()

    today = dt.date.today()
    avail = oisst_files_between(today - dt.timedelta(days=45), today)
    end = dt.date.fromisoformat(args.end) if args.end else dt.datetime.strptime(max(avail), "%Y%m%d").date()
    start = end - dt.timedelta(days=args.days - 1)
    print("fenêtre", start, "->", end)

    recent = {k: v for k, v in oisst_files_between(start, end).items()}
    recent_paths = fetch_oisst(recent, os.path.join(args.cache, "oisst"))

    # 12 mois à pas de 5 jours pour l'ENSO, puis les 30 jours quotidiens
    year_ago = end - dt.timedelta(days=365)
    hist = oisst_files_between(year_ago, start - dt.timedelta(days=1))
    keys = sorted(hist)[::5]
    hist_paths = fetch_oisst({k: hist[k] for k in keys}, os.path.join(args.cache, "oisst"))
    enso_paths = {**hist_paths, **recent_paths}

    ltm_path = get(OISST_LTM, os.path.join(args.cache, "oisst", "sst.day.mean.ltm.1991-2020.nc"))
    dates, ocean1 = build_ocean(recent_paths, ltm_path)
    anom05, clim05 = build_land(start, end, os.path.join(args.cache, "cpc"), args.fc_days)
    land1 = anom05.coarsen(lat=2, lon=2, boundary="exact").mean()
    enso = build_enso(enso_paths, ltm_path)

    # prévision : la mer fixe le nombre de jours, limité aux jours où la terre existe aussi (NaN si la terre échoue)
    ofc = soft("prévision mer GLO12", lambda: ocean_forecast(ocean1, end, ltm_path, os.path.join(args.cache, "glo12"), args.fc_days))
    lfc = soft("prévision terre ECMWF", lambda: land_forecast(anom05, clim05, end, os.path.join(args.cache, "ecmwf"), args.fc_days))
    nobs, nfc = len(dates), 0 if ofc is None else len(ofc) if lfc is None else min(len(ofc), len(lfc))
    if nfc:
        def pad(fc):   # IFS open data s'arrête à 240 h : la terre a un jour de moins que la mer
            out = np.full((nfc, 180, 360), np.nan)
            if fc is not None:
                out[:min(nfc, len(fc))] = fc[:nfc]
            return out
        stack = lambda obs, fc: xr.concat([obs, obs.isel(time=[0] * nfc).copy(data=pad(fc))], "time")
        ocean1, land1 = stack(ocean1, ofc), stack(land1, lfc)
        dates += [(end + dt.timedelta(days=k)).isoformat() for k in range(1, nfc + 1)]

    meta = dict(dates=dates, nobs=nobs, nlat=180, nlon=360, lat0=-89.5, lon0=-179.5, step=1.0,
                fc_land=lfc is not None, **series(ocean1, land1))
    enso["probs"] = soft("probabilités ENSO CPC", surveillance.enso_probs)
    surv = dict(ice=soft("NSIDC", surveillance.sea_ice), coral=soft("Coral Reef Watch", surveillance.coral),
                ice_maps=soft("cartes de glace", lambda: ice_maps(recent_paths[max(recent_paths)], end, args.cache)))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(dict(meta=meta, ocean_b64=encode_cube(ocean1.values), land_b64=encode_cube(land1.values, clip=127),
                   enso=enso, surv=surv, built_at=dt.datetime.now(dt.timezone.utc).isoformat()), open(args.out, "w"))
    t = nobs - 1
    print("ok ->", args.out, "| océan 60S-60N", meta["g60"][t], "| terre", meta["lmean"][t], "| Niño 3.4", enso["hov"]["n34"][-1],
          "| prévision", nfc, "j (terre", "ok" if lfc is not None else "absente", ")")


if __name__ == "__main__":
    sys.exit(main())
