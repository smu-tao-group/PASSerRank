#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math, glob, pickle
import numpy as np
from tqdm import tqdm
from utils.count_atoms import count_atoms


DISTANCE_THRESHOLD = 10

pdb_diverse = pickle.load(open("../data/clean_data/pdb_info.pkl", "rb"))

# nearest pocket distance for each pdb
min_dists = []
selected_pdbs = []

for pdb, info in tqdm(pdb_diverse.items()):
    modulator, mod_id, chain, residues, _ = info.values()
    protein = open(f"../data/pdbs/{pdb}.pdb", "r").readlines()

    # ligand center of mass
    lig_x, lig_y, lig_z, lig_cnt = 0, 0, 0, 0
    for line in protein:
        if (
            line[:6] == "HETATM" and modulator == line[17:20].strip()
            and line[21] == chain and mod_id == line[22:26].strip()
        ):
            lig_x += float(line[30:38])
            lig_y += float(line[38:46])
            lig_z += float(line[46:54])
            lig_cnt += 1

    # drop if no ligand atom found
    if lig_cnt == 0:
        continue

    lig_x /= lig_cnt
    lig_y /= lig_cnt
    lig_z /= lig_cnt

    # collect all pocket
    pocket_dir = f"../data/pockets/{pdb}_out/pockets/"
    pocket_names = glob.glob(pocket_dir + "*.pdb")
    pocket_names = sorted(
        pocket_names,
        key = lambda x : int(x.split("pocket")[-1].split("_")[0])
    )

    # find the nearest pocket
    dists = []
    counts = []
    for idx, pocket_name in enumerate(pocket_names):
        pocket = open(pocket_name, "r").readlines()
        poc_x, poc_y, poc_z, poc_cnt = 0, 0, 0, 0

        for line in pocket:
            if line[:4] == "ATOM":
                poc_x += float(line[30:38])
                poc_y += float(line[38:46])
                poc_z += float(line[46:54])
                poc_cnt += 1

        # drop if no pocket atom
        if poc_cnt == 0:
            continue

        poc_x /= poc_cnt
        poc_y /= poc_cnt
        poc_z /= poc_cnt
        dist = math.sqrt(
            (poc_x - lig_x) ** 2 + (poc_y - lig_y) ** 2 +
            (poc_z - lig_z) ** 2
        )

        dists.append(dist)
        counts.append(count_atoms(pocket_name, residues))

    # drop if there are less than 2 pockets found (abnormal)
    if len(dists) <= 2:
        continue

    dist_min_idx = np.argmin(dists)
    min_dists.append(dists[dist_min_idx])
    selected_pdbs.append(pdb)


print("Number of PDBs: ", len(min_dists))

# for plotting
min_dists = np.array(min_dists)
pickle.dump(min_dists, open("../plots/min_dists.pkl", "wb"))

pdb_distance = {}
for idx, pdb in enumerate(selected_pdbs):
    if min_dists[idx] > DISTANCE_THRESHOLD:
        continue
    pdb_distance[pdb] = pdb_diverse[pdb]
    pdb_distance[pdb]["distance"] = min_dists[idx]


print(
    f"Number of PDBs with less than {DISTANCE_THRESHOLD}A: ",
    len(pdb_distance)
)
pickle.dump(
    pdb_distance, open("../data/clean_data/pdb_distance.pkl", "wb")
)
