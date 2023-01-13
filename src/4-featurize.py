#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob, os, pickle, math
from tqdm import tqdm
import numpy as np
from utils.pocket_feature import pocket_feature


pdb_diverse = pickle.load(open("../data/clean_data/pdb_diverse.pkl", "rb"))
labels = []
features = []

for pdb, info in tqdm(pdb_diverse.items()):
    modulator, mod_id, chain, residues, _, _ = info.values()
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
        key=lambda x: int(x.split("pocket")[-1].split("_")[0])
    )

    # find the nearest pocket
    dists = []

    # collect pocket features and labels
    cur_features = pocket_feature(f"../data/pockets/{pdb}_out/{pdb}_info.txt")
    selected_idxs = []

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
        selected_idxs.append(idx)

    # drop if there are less than 2 pockets found
    if len(dists) <= 2:
        continue

    features.append([cur_features[idx] for idx in selected_idxs])

    dist_min_idx = np.argmin(dists)
    cur_labels = [0] * len(dists)
    cur_labels[dist_min_idx] = 1
    labels.append(cur_labels)


# summarize
total_labels = sum([len(item) for item in labels])
positive_labels = sum([sum(item) for item in labels])
print(
    "total of %d pockets, with %d positive labels accounting for %.2f%%"
    % (total_labels, positive_labels, positive_labels / total_labels * 100)
)

# clear history
os.system("rm -r ../data/classification/*")

# dump data
pickle.dump(labels, open("../data/classification/labels.pkl", "wb"))
pickle.dump(features, open("../data/classification/features.pkl", "wb"))
