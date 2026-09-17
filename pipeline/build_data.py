#!/usr/bin/env python3
"""
Fièvre océan, pipeline quotidien.

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

OISST_BASE = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/"
PSL_FILES = "https://psl.noaa.gov/thredds/fileServer/Datasets/cpc_global_temp/"
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
    """Cube 1° des 30 derniers jours + stats globales."""
    anoms, ssts, dates = [], [], []
    for k in sorted(paths_by_date):
        a, s, d = open_anom(paths_by_date[k], ltm_path)
        anoms.append(to_180(a)); ssts.append(to_180(s)); dates.append(d)
    anom = xr.concat(anoms, "time"); sst = xr.concat(ssts, "time")
    w = np.cos(np.deg2rad(anom.lat))
    band = anom.sel(lat=slice(-60, 60)); wb = w.sel(lat=slice(-60, 60))
    meta = dict(
        dates=dates, nlat=180, nlon=360, lat0=-89.5, lon0=-179.5, step=1.0,
        gmean=[round(float(x), 3) for x in anom.weighted(w).mean(("lat", "lon")).values],
        g60=[round(float(x), 3) for x in band.weighted(wb).mean(("lat", "lon")).values],
        sst60=[round(float(x), 2) for x in sst.sel(lat=slice(-60, 60)).weighted(wb).mean(("lat", "lon")).values],
        frac_hot=[round(float(x), 3) for x in ((anom > 1).where(anom.notnull())).weighted(w).mean(("lat", "lon")).values],
    )
    a1 = anom.coarsen(lat=4, lon=4, boundary="exact").mean()
    return meta, encode_cube(a1.values)


# ---------------------------------------------------------------- terre (CPC)
def build_land(d0, d1, cache):
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
        days = cftime.num2date(ltm.time.values, ltm.time.attrs["units"], "gregorian")   # année 1 : cftime
        i0 = next(k for k, d in enumerate(days) if (d.month, d.day) == (d0.month, d0.day))
        return obs, ltm[var].values[i0:i0 + n]

    tx, lxv = field("tmax")
    tn, ln = field("tmin")
    anom = (tx.values + tn.values) / 2 - (lxv + ln) / 2
    da = xr.DataArray(anom, dims=("time", "lat", "lon"),
                      coords=dict(time=tx.time.values, lat=tx.lat.values, lon=tx.lon.values))
    da = to_180(da).sortby("lat")
    w = np.cos(np.deg2rad(da.lat))
    lmean = da.weighted(w).mean(("lat", "lon")).values
    a1 = da.coarsen(lat=2, lon=2, boundary="exact").mean()
    return [None if np.isnan(x) else round(float(x), 2) for x in lmean], encode_cube(a1.values, clip=127)


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
    meta, ocean_b64 = build_ocean(recent_paths, ltm_path)
    lmean, land_b64 = build_land(start, end, os.path.join(args.cache, "cpc"))
    meta["lmean"] = lmean
    enso = build_enso(enso_paths, ltm_path)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(dict(meta=meta, ocean_b64=ocean_b64, land_b64=land_b64, enso=enso,
                   built_at=dt.datetime.utcnow().isoformat() + "Z"), open(args.out, "w"))
    print("ok ->", args.out, "| océan 60S-60N", meta["g60"][-1], "| terre", lmean[-1], "| Niño 3.4", enso["hov"]["n34"][-1])


if __name__ == "__main__":
    sys.exit(main())
