#!/usr/bin/env python3
"""
Build gw.html: gravitational-wave follow-up with 7DT and KMTNet, from public data.

For each event in tools/lists/gw.csv ("graceid,circulars") the script
  1. reads the superevent, its alert notices and the latest sky map from GraceDB,
  2. reads the GCN Circulars listed for it (titles, and any updated classification),
  3. reads the 7DT and KMTNet pointings and footprints from the Gravitational-wave Treasure Map,
  4. works out the 50% and 90% contours, the 90% area and the probability each instrument covered,
  5. fills gw_template.html and writes gw.html (the sky data are embedded in the page).

Everything it reads is public, so it runs anywhere with internet access:

    pip install astropy astropy-healpix numpy requests matplotlib
    python tools/make_gw.py

Downloaded sky maps are cached in tools/.gwcache/ (git-ignored).
"""
import csv, html, io, json, math, os, re, sys

import numpy as np
import requests

GRACEDB = "https://gracedb.ligo.org"
TM = "https://treasuremap.space"
GCN = "https://gcn.nasa.gov/circulars"
INSTRUMENTS = [(94, "7DT", "7-Dimensional Telescope"), (95, "KMTNet", "KMTNet")]
CACHE = "tools/.gwcache"
ORDER = 9             # HEALPix order for coverage (nside 512, about 0.11 deg)
HEAT_ORDER = 7        # HEALPix order of the probability cells drawn on the globe
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "seowony.github.io gw page builder"


def get(url, binary=False, cache=None):
    if cache and os.path.exists(cache):
        return open(cache, "rb").read()
    r = SESSION.get(url, timeout=120)
    r.raise_for_status()
    data = r.content
    if cache:
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        open(cache, "wb").write(data)
    return data if binary else data.decode()


def get_json(url):
    return json.loads(get(url))


# --------------------------------------------------------------------------- GraceDB
def pick_skymap(files):
    names = [f for f in files if "," not in f and f.endswith(".multiorder.fits")]
    offline = sorted((int(m.group(1)), f) for f in names if (m := re.match(r"Bilby\.offline(\d+)\.multiorder\.fits$", f)))
    if offline:
        return offline[-1][1]
    for f in ("Bilby.multiorder.fits", "bayestar.multiorder.fits"):
        if f in names:
            return f
    return names[-1] if names else None


def latest_notice(gid, files):
    for kind in ("update", "initial", "preliminary"):
        f = f"{gid}-{kind}.json"
        if f in files:
            try:
                return get_json(f"{GRACEDB}/api/superevents/{gid}/files/{f}")
            except Exception:
                continue
    return None


def read_skymap(path):
    """Flat nested probability at ORDER, plus distance mean and std (Mpc)."""
    from astropy.io import fits
    with fits.open(path) as h:
        d, hdr = h[1].data, h[1].header
        uniq = np.asarray(d["UNIQ"], dtype=np.int64)
        dens = np.asarray(d["PROBDENSITY"], dtype=float)
    order = (np.floor(np.log2(uniq) / 2) - 1).astype(int)
    ipix = uniq - 4 * 4 ** order.astype(np.int64)
    npix = 12 * 4 ** ORDER
    acc = np.zeros(npix)
    cnt = np.zeros(npix)
    for o in np.unique(order):
        m = order == o
        if o <= ORDER:
            k = 4 ** (ORDER - o)
            base = ipix[m] * k
            idx = (base[:, None] + np.arange(k)[None, :]).ravel()
            acc[idx] += np.repeat(dens[m], k)
            cnt[idx] += 1
        else:
            parent = ipix[m] >> (2 * (o - ORDER))
            np.add.at(acc, parent, dens[m] / 4 ** (o - ORDER))
            np.add.at(cnt, parent, 1 / 4 ** (o - ORDER))
    dens_map = np.where(cnt > 0, acc / np.maximum(cnt, 1e-30), 0)
    prob = dens_map * (4 * math.pi / npix)
    prob /= prob.sum()
    return prob, hdr.get("DISTMEAN"), hdr.get("DISTSTD")


def credible_levels(prob):
    order = np.argsort(prob)[::-1]
    cl = np.empty_like(prob)
    cl[order] = np.cumsum(prob[order])
    return cl


def contours(cl, levels=(0.5, 0.9), step=0.25):
    """Contour lines of the credible-level map as lists of [ra, dec] (deg)."""
    import contourpy
    from astropy_healpix import HEALPix
    import astropy.units as u
    hp = HEALPix(nside=2 ** ORDER, order="nested")
    ra = np.arange(0, 360 + step, step)
    dec = np.arange(-90, 90 + step, step)
    RA, DEC = np.meshgrid(ra, dec)
    pix = hp.lonlat_to_healpix(RA.ravel() * u.deg, DEC.ravel() * u.deg).reshape(RA.shape)
    gen = contourpy.contour_generator(ra, dec, cl[pix])
    out = {}
    for lev in levels:
        lines = []
        for seg in gen.lines(lev):
            if len(seg) < 4:
                continue
            seg = seg[:: max(1, len(seg) // 400)]
            lines.append([[round(float(x), 2), round(float(y), 2)] for x, y in seg])
        out[str(int(lev * 100))] = lines
    return out


def heat_cells(prob, cl, limit=0.9, max_cells=3000):
    """HEALPix cells inside the 90% region (order chosen to keep the page small): corners and shade 0-1."""
    from astropy_healpix import HEALPix
    for order in range(HEAT_ORDER, 3, -1):
        k = 4 ** (ORDER - order)
        p = prob.reshape(-1, k).sum(axis=1)
        c = cl.reshape(-1, k).min(axis=1)
        sel = np.where(c <= limit)[0]
        if sel.size <= max_cells:
            break
    if sel.size == 0:
        return []
    hp = HEALPix(nside=2 ** order, order="nested")
    lon, lat = hp.boundaries_lonlat(sel, step=1)
    ref = np.quantile(p[sel], 0.98)
    shade = np.sqrt(np.clip(p[sel] / ref, 0, 1))
    return [[[[round(float(lon[i][j].deg), 2), round(float(lat[i][j].deg), 2)] for j in range(4)],
             round(float(shade[i]), 2)] for i in range(sel.size)]


# --------------------------------------------------------------------------- Treasure Map
def parse_polygon(wkt):
    nums = [float(x) for x in re.findall(r"-?\d+\.?\d*", wkt)]
    pts = list(zip(nums[0::2], nums[1::2]))
    return pts[:-1] if pts and pts[0] == pts[-1] else pts


def tm_pointings(gid, inst):
    data = get_json(f"{TM}/api/v1/pointings?graceid={gid}&instrument={inst}")
    out = []
    for p in data:
        if p.get("status") != "completed":
            continue
        m = re.match(r"POINT\s*\(\s*(-?[\d.]+)\s+(-?[\d.]+)\s*\)", p["position"])
        if not m:
            continue
        ra, dec = float(m.group(1)) % 360, float(m.group(2))
        out.append(dict(ra=ra, dec=dec, pa=float(p.get("pos_angle") or 0), band=p.get("band") or "",
                        depth=p.get("depth"), time=(p.get("time") or "")[:19]))
    return out


def tm_footprint(inst):
    return [parse_polygon(f["footprint"]) for f in get_json(f"{TM}/api/v1/footprints?id={inst}")]


def coverage(prob, pointings, chips):
    """Total probability inside the union of all footprints."""
    from astropy_healpix import HEALPix
    from matplotlib.path import Path
    import astropy.units as u
    if not pointings:
        return 0.0
    hp = HEALPix(nside=2 ** ORDER, order="nested")
    radius = max(math.hypot(x, y) for c in chips for x, y in c) * 1.05
    paths = [Path(c) for c in chips]
    hit = np.zeros(prob.size, dtype=bool)
    done = set()
    for p in pointings:
        key = (round(p["ra"], 4), round(p["dec"], 4), round(p["pa"], 1))
        if key in done:
            continue
        done.add(key)
        idx = hp.cone_search_lonlat(p["ra"] * u.deg, p["dec"] * u.deg, radius * u.deg)
        if idx.size == 0:
            continue
        lon, lat = hp.healpix_to_lonlat(idx)
        ra, dec = np.radians(lon.deg), np.radians(lat.deg)
        ra0, dec0 = math.radians(p["ra"]), math.radians(p["dec"])
        cosc = math.sin(dec0) * np.sin(dec) + math.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
        x = np.degrees(np.cos(dec) * np.sin(ra - ra0) / cosc)      # east
        y = np.degrees((math.cos(dec0) * np.sin(dec) - math.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc)
        t = math.radians(p["pa"])
        xr, yr = x * math.cos(t) + y * math.sin(t), -x * math.sin(t) + y * math.cos(t)
        pts = np.column_stack([xr, yr])
        inside = np.zeros(idx.size, dtype=bool)
        for path in paths:
            inside |= path.contains_points(pts)
        hit[idx[inside & (cosc > 0)]] = True
    return float(prob[hit].sum())


# --------------------------------------------------------------------------- GCN
def circular(cid):
    try:
        d = get_json(f"{GCN}/{cid}.json")
        return dict(id=cid, title=d.get("subject", ""), body=d.get("body", ""))
    except Exception as e:
        print(f"  GCN {cid}: {e}")
        return dict(id=cid, title=f"GCN Circular {cid}", body="")


def updated_classification(circs):
    """The last 'updated classification' sentence in the circulars, if any."""
    found = None
    for c in circs:
        m = re.search(r"updated classification of the candidate, in order of descending probability, is ([^.]*)\.",
                      c["body"], flags=re.I)
        if m:
            found = (m.group(1).strip(), c["id"])
    return found


# --------------------------------------------------------------------------- page
def pct(x):
    return "<1%" if x < 0.01 else (">99%" if x > 0.99 else f"{x * 100:.0f}%")


def class_text(cls):
    names = {"BNS": "BNS", "NSBH": "NSBH", "BBH": "BBH", "Terrestrial": "terrestrial"}
    items = sorted(cls.items(), key=lambda kv: -kv[1])
    return ", ".join(f"{names.get(k, k)} {pct(v)}" for k, v in items if v >= 0.005)


def fmt_dates(times):
    days = sorted({t[:10] for t in times if t})
    if not days:
        return ""
    return days[0] if len(days) == 1 else f"{days[0]} to {days[-1]}"


def build_event(gid, circ_ids):
    from astropy.time import Time
    print(f"== {gid}")
    se = get_json(f"{GRACEDB}/api/superevents/{gid}/")
    files = get_json(f"{GRACEDB}/api/superevents/{gid}/files/")
    t0 = Time(se["t_0"], format="gps").utc.iso[:19]
    notice = latest_notice(gid, files) or {}
    ev = notice.get("event") or {}
    cls = ev.get("classification") or {}

    # latest sky map first, then the first low-latency map (what the observations were planned on)
    names = [pick_skymap(files)]
    if "bayestar.multiorder.fits" in files and names[0] != "bayestar.multiorder.fits":
        names.append("bayestar.multiorder.fits")
    maps = []
    for sky in names:
        path = os.path.join(CACHE, f"{gid}_{sky}")
        get(f"{GRACEDB}/api/superevents/{gid}/files/{sky}", binary=True, cache=path)
        prob, dmean, dstd = read_skymap(path)
        cl = credible_levels(prob)
        pixarea = 4 * math.pi * (180 / math.pi) ** 2 / prob.size
        area90 = float((cl <= 0.9).sum() * pixarea)
        maps.append(dict(name=sky.replace(".multiorder.fits", ""), prob=prob, cl=cl, area90=area90,
                         dist=[dmean, dstd]))
        print(f"  sky map {sky}: 90% area {area90:.0f} deg2, distance {dmean:.0f} +/- {dstd:.0f} Mpc")

    circs = [circular(c) for c in circ_ids]
    upd = updated_classification(circs)

    insts = []
    for iid, short, full in INSTRUMENTS:
        pts = tm_pointings(gid, iid)
        if not pts:
            continue
        chips = tm_footprint(iid)
        cov = [coverage(m["prob"], pts, chips) for m in maps]
        depths = [p["depth"] for p in pts if p["depth"]]
        print(f"  {short}: {len(pts)} pointings, " + ", ".join(f"{c * 100:.1f}% of {m['name']}" for c, m in zip(cov, maps)))
        insts.append(dict(id=iid, short=short, name=full, n=len(pts), coverage=cov,
                          bands=sorted({p["band"] for p in pts if p["band"]}),
                          dates=fmt_dates([p["time"] for p in pts]),
                          depth=float(np.median(depths)) if depths else None,
                          chips=[[[round(x, 4), round(y, 4)] for x, y in c] for c in chips],
                          pointings=[[round(p["ra"], 4), round(p["dec"], 4), round(p["pa"], 1), p["band"],
                                      round(p["depth"], 2) if p["depth"] else None, p["time"]] for p in pts]))
    return dict(id=gid, t0=t0, far=se.get("far"), classification=cls, updated=upd,
                maps=[dict(name=m["name"], area90=m["area90"], dist=m["dist"], contours=contours(m["cl"]),
                           cells=heat_cells(m["prob"], m["cl"])) for m in maps],
                instruments=insts,
                circulars=[dict(id=c["id"], title=c["title"]) for c in circs])


def card(e, i):
    esc = html.escape
    lines = []
    cls = class_text(e["classification"]) if e["classification"] else ""
    if e["updated"]:
        txt, cid = e["updated"]
        cls_html = (f"At first {esc(cls)}. After an updated significance estimate: {esc(txt)} "
                    f'(<a href="{GCN}/{cid}">GCN {cid}</a>).') if cls else f"{esc(txt)} (GCN {cid})."
    else:
        cls_html = f"{esc(cls)}." if cls else ""
    m0 = e["maps"][0]
    dist = f"{m0['dist'][0]:.0f} ± {m0['dist'][1]:.0f} Mpc" if m0["dist"][0] else ""
    labels = ["latest map"] + ["first alert map"] * (len(e["maps"]) - 1)
    ths = "".join(f"<th>Covered ({esc(m['name'])})</th>" for m in e["maps"])
    rows = []
    for inst in e["instruments"]:
        depth = f"{inst['depth']:.1f}" if inst["depth"] else ""
        rows.append(f'          <tr><td><span class="key k{inst["id"]}"></span>{esc(inst["short"])}</td>'
                    f'<td>{inst["n"]}</td><td>{esc(inst["dates"])}</td><td>{esc(", ".join(inst["bands"]))}</td>'
                    f'<td>{depth}</td>' + "".join(f'<td>{c * 100:.1f}%</td>' for c in inst["coverage"]) + '</tr>')
    circ = "\n".join(f'          <li><a href="{GCN}/{c["id"]}">GCN {c["id"]}</a> {esc(c["title"])}</li>'
                     for c in e["circulars"])
    tm_names = ", ".join(i["short"] for i in e["instruments"])
    return f'''
    <article class="gw-event" id="{e["id"]}">
      <h2>{e["id"]}</h2>
      <p class="gw-meta">{e["t0"]} UTC. {cls_html}</p>
      <dl class="gw-facts">
        <div><dt>Distance</dt><dd>{dist}</dd></div>
        <div><dt>90% area</dt><dd>{m0["area90"]:.0f} deg²</dd></div>
        <div><dt>Latest sky map</dt><dd>{esc(m0["name"])}</dd></div>
      </dl>
      <div class="gw-sky" data-event="{i}">
        <svg role="img" aria-label="Sky map of {e["id"]} with the {tm_names} pointings"></svg>
        <div class="gw-controls">
{"".join(f'          <label><input type="radio" name="map{i}" value="{k}"{" checked" if k == 0 else ""}> {esc(m["name"])} ({labels[k]}, {m["area90"]:.0f} deg²)</label>{chr(10)}' for k, m in enumerate(e["maps"])) if len(e["maps"]) > 1 else ""}          <label><input type="checkbox" data-layer="heat" checked> Probability</label>
          <label><input type="checkbox" data-layer="c90" checked> 90% / 50%</label>
{"".join(f'          <label><input type="checkbox" data-layer="i{inst["id"]}" checked><span class="key k{inst["id"]}"></span>{esc(inst["short"])}</label>{chr(10)}' for inst in e["instruments"])}          <button type="button" data-view="field">Field</button>
          <button type="button" data-view="sky">Whole sky</button>
          <span class="gw-readout" aria-live="polite"></span>
        </div>
      </div>
      <div class="gw-table-wrap">
        <table class="gw-table">
          <thead><tr><th>Telescope</th><th>Pointings</th><th>Dates (UTC)</th><th>Filters</th><th>Median depth (AB)</th>{ths}</tr></thead>
          <tbody>
{chr(10).join(rows)}
          </tbody>
        </table>
      </div>
      <p class="gw-links"><a href="{TM}/alerts?graceids={e["id"]}">Treasure Map</a>
        <a href="{GRACEDB}/superevents/{e["id"]}/view/">GraceDB</a></p>
      <details class="gw-circs">
        <summary>GCN Circulars</summary>
        <ul>
{circ}
        </ul>
      </details>
    </article>'''


def main():
    if not os.path.exists("gw_template.html"):
        sys.exit("Run this from the website root (gw_template.html not found).")
    with open("tools/lists/gw.csv") as f:
        rows = list(csv.DictReader(f))
    events = [build_event(r["graceid"].strip(), [c for c in r["circulars"].split(";") if c.strip()]) for r in rows]
    events.sort(key=lambda e: e["t0"], reverse=True)
    data = [dict(id=e["id"], maps=[{k: m[k] for k in ("name", "contours", "cells")} for m in e["maps"]],
                 instruments=e["instruments"]) for e in events]
    page = open("gw_template.html").read()
    page = page.replace("<!-- GW_EVENTS -->", "\n".join(card(e, i) for i, e in enumerate(events)))
    page = page.replace("/* GW_DATA */", "const GW = " + json.dumps(data, separators=(",", ":")) + ";")
    open("gw.html", "w").write(page)
    print(f"Wrote gw.html with {len(events)} events ({len(page) / 1e3:.0f} kB).")


if __name__ == "__main__":
    main()
