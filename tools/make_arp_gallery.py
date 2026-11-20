#!/usr/bin/env python3
"""Kept for old commands: same as `python tools/make_atlas.py arp ...`."""
import os, sys
sys.argv.insert(1, "arp")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_atlas import main
main()
