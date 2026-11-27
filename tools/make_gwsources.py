#!/usr/bin/env python3
"""
Build gwsources.html: real sources found in the GECKO KMTNet follow-up of GW alerts, with light curves.

Input: the folder written on the SNU server by
    python manage.py export_realsources S240422ed S250206dm --out realsources/
(sources.csv, <id>_{new,ref,sub}.png and lc/<id>.csv). The folder is not tracked by git.

    python tools/make_gwsources.py --input realsources

Writes assets/gwsources/<id>.jpg (new | ref | diff) and gwsources.html from gwsources_template.html. Light
curves are drawn as inline SVG, in days since the GW event (event times from GraceDB, cached in tools/.gwcache).
Only sources marked reported=yes in sources.csv are used unless --include-unreported is given.
"""
import argparse, csv, html, json, os, re, sys

import requests
from PIL import Image

OUT = "assets/gwsources"
CACHE = "tools/.gwcache/event_t0.json"
TYPES = {"YSO": "young stellar object", "CV": "cataclysmic variable", "II": "Type II supernova",
         "Ia": "Type Ia supernova", "IIn": "Type IIn supernova", "IIP": "Type IIP supernova", "Ib": "Type Ib supernova",
         "Ic": "Type Ic supernova", "QSO": "quasar", "AGN": "active galactic nucleus", "Mira": "Mira variable",
         "RRLyr": "RR Lyrae star", "LPV": "long-period variable", "EB": "eclipsing binary"}
COLOURS = {"B": "#2f62c9", "V": "#2e8b57", "R": "#c8402f", "I": "#7a1f17", "g": "#2e8b57", "r": "#c8402f",
           "i": "#7a1f17", "z": "#4b2a1f"}


def event_t0_mjd(event):
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    if event not in cache:
        from astropy.time import Time
        d = requests.get(f"https://gracedb.ligo.org/api/superevents/{event}/", timeout=60).json()
        cache[event] = Time(d["t_0"], format="gps").mjd
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        json.dump(cache, open(CACHE, "w"), indent=1)
    return cache[event]


def link_for(name):
    if name.startswith("Gaia"):
        return f"https://gsaweb.ast.cam.ac.uk/alerts/alert/{name}/"
    m = re.match(r"^(?:AT|SN)\s?(\d{4}[a-z]+)$", name)
    if m:
        return f"https://www.wis-tns.org/object/{m.group(1)}"
    return None


def pretty_name(name):
    m = re.match(r"^(AT|SN)\s?(\d{4}[a-z]+)$", name)
    return f"{m.group(1)} {m.group(2)}" if m else name


def lc_svg(rows, t0, col="mag_aper"):
    """SVG light curve from the `col` magnitudes (mag_aper or mag_auto; old files have mag/mag_err)."""
    if rows and col not in rows[0] and "mag" in rows[0]:
        col = "mag"
    pts = [(float(r["mjd"]) - t0, float(r[col]), float(r.get(col + "_err") or r.get("mag_err") or 0), r["filter"])
           for r in rows if r.get(col)]
    if not pts:
        return ""
    W, H, L, R, T, B = 460, 190, 44, 12, 10, 30
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    x0, x1 = min(xs), max(xs)
    if x1 - x0 < 0.5:
        x0, x1 = x0 - 0.25, x1 + 0.25
    pad = (x1 - x0) * 0.06
    x0, x1 = x0 - pad, x1 + pad
    y0, y1 = min(ys) - 0.3, max(ys) + 0.3                # magnitudes: bright at the top
    X = lambda x: L + (x - x0) / (x1 - x0) * (W - L - R)
    Y = lambda y: T + (y - y0) / (y1 - y0) * (H - T - B)
    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Light curve">',
           f'<line class="ax" x1="{L}" y1="{H - B}" x2="{W - R}" y2="{H - B}"/>',
           f'<line class="ax" x1="{L}" y1="{T}" x2="{L}" y2="{H - B}"/>']
    step = max(0.1, round((y1 - y0) / 4, 1))
    m = round(y0 + 0.05, 1)
    while m <= y1:
        out.append(f'<text x="{L - 6}" y="{Y(m) + 4:.1f}" text-anchor="end">{m:.1f}</text>')
        m = round(m + step, 1)
    nice = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200]
    tstep = next(s for s in nice if (x1 - x0) / s <= 6)
    t = (int(x0 / tstep) + (1 if x0 > 0 else 0)) * tstep
    while t <= x1:
        out.append(f'<text x="{X(t):.1f}" y="{H - B + 14}" text-anchor="middle">{t:g}</text>')
        t += tstep
    out.append(f'<text x="{(L + W - R) / 2:.0f}" y="{H - 2}" text-anchor="middle">days after the GW event</text>')
    out.append(f'<text x="11" y="{(T + H - B) / 2:.0f}" text-anchor="middle" '
               f'transform="rotate(-90 11 {(T + H - B) / 2:.0f})">mag</text>')
    filters = []
    for x, y, e, f in pts:
        c = COLOURS.get(f, "#555")
        if f not in filters:
            filters.append(f)
        if e:
            out.append(f'<line x1="{X(x):.1f}" y1="{Y(y - e):.1f}" x2="{X(x):.1f}" y2="{Y(y + e):.1f}" '
                       f'stroke="{c}" stroke-width="1.2"/>')
        out.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="3.2" fill="{c}"><title>{f} {y:.2f} ± {e:.2f}, '
                   f'{x:+.2f} d</title></circle>')
    for k, f in enumerate(filters):
        out.append(f'<circle cx="{W - R - 60 + k * 18}" cy="{T + 6}" r="3.5" fill="{COLOURS.get(f, "#555")}"/>'
                   f'<text x="{W - R - 54 + k * 18}" y="{T + 10}">{html.escape(f)}</text>')
    out.append("</svg>")
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="realsources")
    ap.add_argument("--size", type=int, default=240, help="pixels per stamp")
    ap.add_argument("--include-unreported", action="store_true")
    args = ap.parse_args()
    if not os.path.exists("gwsources_template.html"):
        sys.exit("Run this from the website root (gwsources_template.html not found).")
    table = os.path.join(args.input, "sources.csv")
    if not os.path.exists(table):
        sys.exit(f"{table} not found")
    os.makedirs(OUT, exist_ok=True)
    sources = [r for r in csv.DictReader(open(table))
               if args.include_unreported or r.get("reported", "yes") == "yes"]
    cards = []
    events = []
    for r in sources:
        stamps = []
        for k in ("new", "ref", "diff"):
            im = Image.open(os.path.join(args.input, r[k])).convert("RGB")
            stamps.append(im.resize((args.size, args.size), Image.LANCZOS))
        for k, im in zip(("new", "ref", "diff"), stamps):
            im.save(os.path.join(OUT, f"{r['id']}_{k}.jpg"), quality=86)
        lc = list(csv.DictReader(open(os.path.join(args.input, r["lc"])))) if r.get("lc") else []
        t0 = event_t0_mjd(r["event"])
        if r["event"] not in events:
            events.append(r["event"])
        name = pretty_name(r["name"])
        url = link_for(r["name"])
        title = f'<a href="{url}">{html.escape(name)}</a>' if url else html.escape(name)
        typ = TYPES.get(r.get("type", ""), r.get("type", ""))
        typ = "" if typ.lower() in ("", "candidate") else typ
        n = len(lc)
        sub = ", ".join(x for x in [typ[:1].upper() + typ[1:] if typ else "", f"{r['event']} follow-up",
                                    f"KMTNet {r['site'].replace('KMTN_', '')}" if r.get("site") else "",
                                    f"{n} detection{'s' if n != 1 else ''}"] if x)
        imgs = "".join(f'<div><img src="{OUT}/{r["id"]}_{k}.jpg" width="{args.size}" height="{args.size}" '
                       f'loading="lazy" alt="{lab} image of {html.escape(name)}"><span>{lab}</span></div>'
                       for k, lab in (("new", "New"), ("ref", "Reference"), ("diff", "Difference")))
        cards.append(f'      <article class="src-card" id="{r["id"]}">\n'
                     f'        <h2>{title}</h2>\n        <p class="src-sub">{html.escape(sub)}</p>\n'
                     f'        <div class="src-stamps">{imgs}</div>\n'
                     f'        <div class="src-lc" data-mag="aper">{lc_svg(lc, t0, "mag_aper")}</div>\n'
                     f'        <div class="src-lc" data-mag="auto">{lc_svg(lc, t0, "mag_auto")}</div>\n      </article>')
    if not cards:
        sys.exit("No sources: gwsources.html was not written.")
    page = open("gwsources_template.html").read()
    page = page.replace("<!-- SRC_CARDS -->", '    <div class="src-list">\n' + "\n".join(cards) + "\n    </div>")
    ev_txt = ", ".join(events[:-1]) + " and " + events[-1] if len(events) > 1 else events[0]
    page = page.replace("<!-- SRC_COUNT -->", f"{len(cards)} sources from the follow-up of {ev_txt}.")
    page = page.replace("<!-- SRC_NOTE -->",
                        "Magnitudes come from the KMTNet follow-up images as measured by the GECKO pipeline: either "
                        "the best-aperture magnitude or the AUTO (Kron-like) magnitude, chosen with the switch above. "
                        "Names link to the Gaia alerts and TNS pages. None of these sources is a counterpart of the "
                        "gravitational-wave event.")
    open("gwsources.html", "w").write(page)
    print(f"Wrote gwsources.html with {len(cards)} sources.")


if __name__ == "__main__":
    main()
