#!/usr/bin/env python
# -*- coding: utf-8 -*-

def count_atoms(pocket_dir, residues):
    pocket = open(pocket_dir, "r").readlines()
    count = 0

    for line in pocket:
        if line[:4] == "ATOM":
            res_name = line[17:20]
            res_id = line[22:26].strip()
            chain_id = line[21]

            if [chain_id, res_name, res_id] in residues:
                count += 1

    return count
