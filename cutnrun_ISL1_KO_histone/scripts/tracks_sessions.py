#!/usr/bin/env python3
"""IGV session files for the browser tracks (scripts/23c).

Usage: python3 tracks_sessions.py <results/tracks>

One session per mark (ISL1, H3K27ac, H3K27me3) in <tracks>/, with paths relative to it:
replicate-mean signal per group (same autoscale group per mark, so heights are comparable), the
matching IgG means, then the bed tracks (peaks, differential regions). Ctr blue, KO red, IgG grey.
"""
import os
import sys
from xml.sax.saxutils import quoteattr

T = sys.argv[1]
DAYS = ["D9", "D15", "D18", "D22"]
BLUE, RED, GREY = "31,119,180", "214,39,40", "120,120,120"


def bw_track(path, name, color, group):
    return (f'    <Track attributeKey={quoteattr(name)} autoScale="true" autoscaleGroup="{group}" '
            f'clazz="org.broad.igv.track.DataSourceTrack" color="{color}" fontSize="10" id={quoteattr(path)} '
            f'name={quoteattr(name)} renderer="BAR_CHART" visible="true" windowFunction="mean"/>')


def bed_track(path, name):
    return (f'    <Track attributeKey={quoteattr(name)} clazz="org.broad.igv.track.FeatureTrack" displayMode="COLLAPSED" '
            f'fontSize="10" id={quoteattr(path)} name={quoteattr(name)} visible="true"/>')


def session(mark, locus, entries):
    res = "\n".join(f'    <Resource path={quoteattr(p)} type="{"bw" if p.endswith(".bw") else "bed"}"/>' for p, _ in entries)
    tracks = "\n".join(t for _, t in entries)
    out = os.path.join(T, f"igv_session_{mark}.xml")
    with open(out, "w") as fh:
        fh.write(f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<Session genome="hg38" locus="{locus}" version="8">\n'
                 f'  <Resources>\n{res}\n  </Resources>\n  <Panel name="DataPanel">\n{tracks}\n  </Panel>\n</Session>\n')
    print(f"{out}: {len(entries)} tracks")


def have(p):
    return os.path.exists(os.path.join(T, p))


# ISL1: control means per day, their IgG means, both kept peak sets
e = []
for d in DAYS:
    p = f"groups/ISL1/{d}_Ctr_ISL1.mean.bw"
    if have(p): e.append((p, bw_track(p, f"{d} Ctr ISL1", BLUE, 1)))
for d in DAYS:
    p = f"groups/IgG/{d}_Ctr_ISL1_IgG.mean.bw"
    if have(p): e.append((p, bw_track(p, f"{d} IgG (ISL1)", GREY, 1)))
for s in ["macs3_ge2", "idr_any_pair", "macs3_2of2"]:
    for d in DAYS:
        p = f"beds/ISL1/{d}_Ctr_ISL1.{s}.bed"
        if have(p): e.append((p, bed_track(p, f"{d} ISL1 {s}")))
session("ISL1", "chr5:51300000-51500000", e)

for mark, locus in [("H3K27ac", "chr5:51300000-51500000"), ("H3K27me3", "chr5:51000000-52000000")]:
    e = []
    for d in DAYS:
        for line, col, lab in [("H3CTR", BLUE, "Ctr"), ("ISL1KO", RED, "KO")]:
            p = f"groups/{mark}/{d}_{line}_{mark}.mean.bw"
            if have(p): e.append((p, bw_track(p, f"{d} {lab} {mark}", col, 1)))
    for d in DAYS:
        for line, lab in [("H3CTR", "Ctr"), ("ISL1KO", "KO")]:
            p = f"groups/IgG/{d}_{line}_IgG.mean.bw"
            if have(p): e.append((p, bw_track(p, f"{d} {lab} IgG", GREY, 1)))
    for d in DAYS:
        for line, lab in [("H3CTR", "Ctr"), ("ISL1KO", "KO")]:
            p = f"beds/{mark}/{d}_{line}_{mark}.peaks.bed"
            if have(p): e.append((p, bed_track(p, f"{d} {lab} {mark} peaks")))
    for t in ["KO_vs_Ctr_D15", "KO_vs_Ctr_D18", "KO_vs_Ctr_D22", "KO_vs_Ctr_shared",
              "time_H3CTR_D18_vs_D15", "time_H3CTR_D22_vs_D18", "time_H3CTR_D22_vs_D15"]:
        p = f"beds/{mark}/diff_{t}.bed"
        if have(p): e.append((p, bed_track(p, f"{mark} diff {t} (red up, blue down)")))
    session(mark, locus, e)
