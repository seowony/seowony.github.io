#!/usr/bin/env python3
"""
Write the object lists used by tools/make_atlas.py.

  tools/lists/ngc.csv     NGC galaxies with major axis >= --min-size arcmin and Dec < --max-dec
  tools/lists/phangs.csv  the PHANGS-MUSE / PHANGS-JWST Cycle 1 galaxies with Dec < --max-dec

Positions and sizes come from OpenNGC (https://github.com/mattiaverga/OpenNGC, CC BY-SA 4.0).
tools/lists/clusters.csv is written by hand; its positions are resolved by name at run time.

    python tools/make_lists.py
"""
import argparse, csv, io, os
import requests

OPENNGC = "https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files/NGC.csv"
GALAXY_TYPES = {"G", "GPair", "GTrpl", "GGroup"}
PHANGS = ["NGC0628", "NGC1087", "NGC1300", "NGC1365", "NGC1385", "NGC1433", "NGC1512",
          "NGC1566", "NGC1672", "NGC2835", "NGC3351", "NGC3627", "NGC4254", "NGC4303",
          "NGC4321", "NGC4535", "NGC5068", "NGC7496", "IC5332"]


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

    with open("tools/lists/phangs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        kept = [n for n in PHANGS if sexa(rows[n]["Dec"], False) < args.max_dec]
        for n in kept:
            w.writerow(row(rows[n], 5.0, 20.0))
    print(f"tools/lists/phangs.csv: {len(kept)} of {len(PHANGS)} galaxies")


if __name__ == "__main__":
    main()
