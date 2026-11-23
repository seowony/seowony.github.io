#!/usr/bin/env python3
"""
Write the object lists used by tools/make_atlas.py.

  tools/lists/ngc.csv       NGC galaxies with major axis >= --min-size arcmin and Dec < --max-dec
  tools/lists/messier.csv   Messier objects with Dec < --max-dec
  tools/lists/caldwell.csv  Caldwell objects with Dec < --max-dec (numbers and names from Wikipedia)

Positions and sizes come from OpenNGC (https://github.com/mattiaverga/OpenNGC, CC BY-SA 4.0); objects that are
not in OpenNGC are resolved through Sesame.
tools/lists/clusters.csv is written by hand; its positions are resolved by name at run time.

    python tools/make_lists.py
"""
import argparse, csv, io, os
import requests

OPENNGC = "https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files/NGC.csv"
ADDENDUM = "https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files/addendum.csv"
CALDWELL = "https://en.wikipedia.org/w/index.php?title=Caldwell_catalogue&action=raw"
GALAXY_TYPES = {"G", "GPair", "GTrpl", "GGroup"}


def sexa(s, hours):
    sign = -1 if s.startswith("-") else 1
    d, m, x = s.lstrip("+-").split(":")
    v = sign * (int(d) + int(m) / 60 + float(x) / 3600)
    return v * 15 if hours else v


def pretty(name):
    import re
    m = re.match(r"(NGC|IC)0*(\d+)(.*)$", name)
    suffix = m.group(3).strip() if m else ""
    if suffix.startswith("NED"):
        suffix = " " + suffix
    return f"{m.group(1)} {m.group(2)}{suffix}" if m else name


def fov_for(major, lo, hi):
    f = max(lo, min(hi, 1.8 * major))
    return round(f * 2) / 2


def row(r, lo, hi):
    major = float(r["MajAx"]) if r["MajAx"] else 0.0
    return {"id": r["Name"].lower().replace(" ", ""), "name": pretty(r["Name"]), "query": pretty(r["Name"]),
            "ra": f"{sexa(r['RA'], True):.5f}", "dec": f"{sexa(r['Dec'], False):+.5f}",
            "fov_arcmin": fov_for(major, lo, hi)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-size", type=float, default=2.0, help="NGC: minimum major axis (arcmin)")
    ap.add_argument("--max-dec", type=float, default=-28.0, help="keep Dec below this (deg)")
    args = ap.parse_args()

    text = requests.get(OPENNGC, timeout=60).text
    rows = {r["Name"]: r for r in csv.DictReader(io.StringIO(text), delimiter=";")}
    os.makedirs("tools/lists", exist_ok=True)
    cols = ["id", "name", "query", "ra", "dec", "fov_arcmin"]

    ngc = [r for r in rows.values()
           if r["Name"].startswith("NGC") and r["Type"] in GALAXY_TYPES and r["Dec"] and r["MajAx"]
           and float(r["MajAx"]) >= args.min_size and sexa(r["Dec"], False) < args.max_dec]
    ngc.sort(key=lambda r: r["Name"])
    with open("tools/lists/ngc.csv", "w", newline="") as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        for r in ngc:
            w.writerow(row(r, 4.0, 20.0))
    print(f"tools/lists/ngc.csv: {len(ngc)} galaxies")

    text = requests.get(ADDENDUM, timeout=60).text
    allrows = dict(rows)
    allrows.update({r["Name"]: r for r in csv.DictReader(io.StringIO(text), delimiter=";")})
    by_desig = {pretty(n): r for n, r in allrows.items() if r["Type"] != "Dup"}

    # Messier: OpenNGC carries the Messier number
    mess = sorted((r for r in allrows.values() if r.get("M") and r["Type"] != "Dup" and r["Dec"]
                   and sexa(r["Dec"], False) < args.max_dec), key=lambda r: int(r["M"]))
    with open("tools/lists/messier.csv", "w", newline="") as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        for r in mess:
            d = row(r, 5.0, 60.0)
            m = int(r["M"])
            common = r.get("Common names", "").split(",")[0].strip()
            d.update(id=f"m{m:03d}", name=f"M{m} ({pretty(r['Name'])}{', ' + common if common else ''})")
            w.writerow(d)
    print(f"tools/lists/messier.csv: {len(mess)} objects")

    # Caldwell: number, designation and common name from the Wikipedia table
    import re
    from astropy.coordinates import SkyCoord
    wiki = requests.get(CALDWELL, timeout=60, headers={"User-Agent": "seowony.github.io list builder"}).text
    untemplate = lambda t: re.sub(r"<ref[^>]*/>|<ref[^>]*>.*?</ref>|\{\{.*?\}\}", "", t)
    link = lambda t: re.sub(r"\[\[(?:[^|\]]*\|)?([^\]]*)\]\]", r"\1", t).replace("''", "").replace("&nbsp;", "").strip()
    cald = []
    for m in re.finditer(r"\|\s*\{\{hs\|\d+\}\}C(\d+)\s*\n\|\s*(.*?)\n\|\s*(.*?)\n", wiki):
        n, desig, common = int(m.group(1)), link(untemplate(m.group(2))), link(untemplate(m.group(3)))
        if desig in ("", "-"):
            continue                                  # objects with no catalogue number (e.g. the Coalsack)
        r = by_desig.get(desig)
        if r and r["Dec"]:
            d = row(r, 5.0, 60.0)
        else:
            try:
                c = SkyCoord.from_name(desig)
            except Exception:
                print(f"  C{n} {desig}: not resolved, skipped")
                continue
            d = {"query": desig, "ra": f"{c.ra.deg:.5f}", "dec": f"{c.dec.deg:+.5f}", "fov_arcmin": 30.0}
        if float(d["dec"]) >= args.max_dec:
            continue
        d.update(id=f"c{n:03d}", name=f"C{n} ({desig}{', ' + common if common else ''})", query=desig)
        cald.append(d)
    with open("tools/lists/caldwell.csv", "w", newline="") as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        w.writerows(cald)
    print(f"tools/lists/caldwell.csv: {len(cald)} objects")



if __name__ == "__main__":
    main()
