#!/usr/bin/env python3
"""IGV session for the ISL1 time-course tracks (step 08c).
  python tracks_session.py <results/tracks>
Writes <tracks>/igv_session.xml with paths relative to it: per-day ISL1 means (one autoscale
group, so heights compare across days), per-day IgG means (own group), then the bed tracks:
consensus peaks per day, the tested sites, differential sites (DiffBind, background normalization).
"""
import os
import sys
from xml.sax.saxutils import quoteattr

T = sys.argv[1]
DAYS = ["D15", "D18", "D22"]
COL = {"D15": "76,120,168", "D18": "245,133,24", "D22": "84,162,75"}
GREY = "120,120,120"


def bw(path, name, color, group):
    return (f'    <Track attributeKey={quoteattr(name)} autoScale="true" autoscaleGroup="{group}" '
            f'clazz="org.broad.igv.track.DataSourceTrack" color="{color}" fontSize="10" id={quoteattr(path)} '
            f'name={quoteattr(name)} renderer="BAR_CHART" visible="true" windowFunction="mean"/>')


def bed(path, name):
    return (f'    <Track attributeKey={quoteattr(name)} clazz="org.broad.igv.track.FeatureTrack" displayMode="COLLAPSED" '
            f'fontSize="10" id={quoteattr(path)} name={quoteattr(name)} visible="true"/>')


entries = []
for d in DAYS:
    entries.append((f"groups/ISL1/{d}_Ctr_ISL1.mean.bw", bw(f"groups/ISL1/{d}_Ctr_ISL1.mean.bw", f"{d} ISL1", COL[d], 1)))
for d in DAYS:
    entries.append((f"groups/IgG/{d}_Ctr_IgG.mean.bw", bw(f"groups/IgG/{d}_Ctr_IgG.mean.bw", f"{d} IgG", GREY, 2)))
for d in DAYS:
    p = f"beds/{d}_Ctr_ISL1.macs3_ge2.bed"
    entries.append((p, bed(p, f"{d} peaks (macs3_ge2)")))
entries.append(("beds/sites.bed", bed("beds/sites.bed", "tested sites (401 bp)")))
for c in ["D18_vs_D15", "D22_vs_D18", "D22_vs_D15"]:
    p = f"beds/diff_{c}.bed"
    entries.append((p, bed(p, f"diff {c.replace('_', ' ')}")))
missing = [p for p, _ in entries if not os.path.exists(os.path.join(T, p))]
if missing:
    sys.exit("missing: " + ", ".join(missing))
res = "\n".join(f'    <Resource path={quoteattr(p)} type="{"bw" if p.endswith(".bw") else "bed"}"/>' for p, _ in entries)
with open(os.path.join(T, "igv_session.xml"), "w") as fh:
    fh.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
             '<Session genome="hg38" locus="chr5:51370000-51410000" version="8">\n'
             f'  <Resources>\n{res}\n  </Resources>\n  <Panel name="DataPanel">\n'
             + "\n".join(t for _, t in entries) + "\n  </Panel>\n</Session>\n")
print(f"igv_session.xml: {len(entries)} tracks")
