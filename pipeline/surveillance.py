"""
Surveillance des océans et prévision ENSO officielle, petites sources texte/NetCDF sans authentification :
  - NSIDC Sea Ice Index v4 (G02135) : étendue quotidienne Arctique et Antarctique depuis 1978
  - NOAA Coral Reef Watch via ERDDAP (NOAA_DHW) : niveaux d'alerte blanchissement, 7 jours max
  - NOAA CPC : probabilités ENSO officielles par saison (table HTML)
Chaque fonction lève en cas d'échec ; build_data.py les appelle en mode dégradé (bloc omis).
"""
import csv, datetime as dt, html, io, re, tempfile, urllib.request
import numpy as np
import xarray as xr

NSIDC = "https://noaadata.apps.nsidc.org/NOAA/G02135/{h}/daily/data/{H}_seaice_extent_daily_v4.0.csv"
CRW = ("https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW.nc?"
       "CRW_BAA_7D_MAX%5B(last)%5D%5B(35):10:(-35)%5D%5B(-180):10:(180)%5D,"
       "CRW_BAA_7D_MAX_mask%5B(last)%5D%5B(35):10:(-35)%5D%5B(-180):10:(180)%5D")
CPC_PROBS = "https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/probabilities/"


def fetch(url, timeout=300):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "climate-observer"}), timeout=timeout).read()


def sea_ice_from_csv(text):
    """Étendue du dernier jour, normale 1991-2020 du même jour, rang depuis 1979 et valeur de chaque année à cette date."""
    ext = {}
    for row in csv.reader(io.StringIO(text)):
        try:
            y, m, d, e, miss = int(row[0]), int(row[1]), int(row[2]), float(row[3]), float(row[4])
        except (ValueError, IndexError):
            continue   # 2 lignes d'en-tête
        if miss < 0.1 and not (m == 2 and d == 29):
            ext[dt.date(y, m, d)] = e
    last = max(ext)
    md = lambda day: (day.month, day.day)
    # 1991-2020 est complet (quotidien depuis 1987), pas de trou à combler
    clim = {}
    for day, e in ext.items():
        if 1991 <= day.year <= 2020:
            clim.setdefault(md(day), []).append(e)
    clim = {k: sum(v) / len(v) for k, v in clim.items()}
    same = []
    for y in range(1979, last.year + 1):   # avant 1987 une mesure tous les 2 jours : jour le plus proche à ±1
        for k in (0, 1, -1):
            day = dt.date(y, last.month, last.day) + dt.timedelta(days=k)
            if day in ext:
                same.append((ext[day], y)); break
    same.sort()
    rank = 1 + [y for _, y in same].index(last.year)
    return dict(date=last.isoformat(), extent=round(ext[last], 3), normal=round(clim[md(last)], 3),
                rank_low=rank, nyears=len(same), same_date={y: round(e, 3) for e, y in same})


def sea_ice():
    return {h: sea_ice_from_csv(fetch(NSIDC.format(h=h, H=h[0].upper())).decode()) for h in ("north", "south")}


def coral():
    """Part de l'océan tropical (35°S-35°N, eau libre, grille 0,5°) par niveau d'alerte CRW (0-4)."""
    with tempfile.NamedTemporaryFile(suffix=".nc") as f:
        f.write(fetch(CRW)); f.flush()
        ds = xr.open_dataset(f.name).load()
    baa = ds.CRW_BAA_7D_MAX.isel(time=0)
    water = ds.CRW_BAA_7D_MAX_mask.isel(time=0) == 0
    w = np.cos(np.deg2rad(ds.latitude)).broadcast_like(baa).where(water & baa.notnull())
    share = lambda sel: round(float(w.where(sel).sum() / w.sum()), 4)
    # ponytail: pas de masque récifs dans NOAA_DHW, on mesure l'océan tropical, pas « les récifs »
    return dict(date=str(ds.time.values[0])[:10], alert1=share(baa >= 3), alert2=share(baa >= 4), warning=share(baa >= 2))


def enso_probs_from_html(page):
    """Table CPC : saison, La Niña, Neutre, El Niño (%)."""
    txt = html.unescape(re.sub(r"<[^>]+>", " ", page))
    rows = re.findall(r"\b([A-Z]{3})\s+[A-Z][a-z]{2} [A-Z][a-z]{2} [A-Z][a-z]{2}\s+(\d{1,3})\s+(\d{1,3})\s+(\d{1,3})\b", txt)
    assert rows, "table CPC introuvable"
    return [[s, int(a), int(b), int(c)] for s, a, b, c in rows]


def enso_probs():
    return enso_probs_from_html(fetch(CPC_PROBS).decode())


if __name__ == "__main__":
    fake = "Year, Month, Day, Extent, Missing, Source\nYYYY, MM, DD, 10^6 sq km, 10^6 sq km, x\n" + "".join(
        f"{y}, 09, 15, {10 - (y - 1979) * 0.1:.3f}, 0.000, x\n" for y in range(1979, 2027))
    fake = fake.replace("1982, 09, 15", "1982, 09, 16")   # trou SMMR : le 16 remplace le 15
    s = sea_ice_from_csv(fake)
    assert s["date"] == "2026-09-15" and s["rank_low"] == 1 and s["nyears"] == 48 and s["same_date"][1982] == 9.7, s
    assert abs(s["normal"] - np.mean([10 - (y - 1979) * 0.1 for y in range(1991, 2021)])) < 1e-3
    p = enso_probs_from_html("<td>ASO</td><td>Aug Sep Oct</td><td>0</td><td>2</td><td>98</td> <td>SON</td><td>Sep Oct Nov</td><td>1</td><td>9</td><td>90</td>")
    assert p == [["ASO", 0, 2, 98], ["SON", 1, 9, 90]], p
    print("self-check ok")
