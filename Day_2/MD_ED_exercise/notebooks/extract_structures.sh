#!/bin/bash
# Extract a .gro structure for each row of selected_structures.txt
# (cluster, time, proj1, proj2, dG) via gmx trjconv -dump.
set -e

group=1  # output group index passed to trjconv (1 = Protein)

mkdir -p low_en_structures/basin_0 low_en_structures/basin_1 low_en_structures/basin_2

tail -n +2 selected_structures.txt | while read -r basin time proj1 proj2 dG; do
    out="low_en_structures/basin_${basin}/conf_${time}_basin_${basin}.gro"
    echo "$group" | gmx trjconv -f prod_c.xtc -s prod.tpr -dump "$time" -o "$out"
    echo "wrote $out"
done

# Combine each basin's structures into one multi-frame file for VMD, same as
# `cat conf_*.gro > selected_basin_N.gro` run from inside each basin folder.
for n in 0 1 2; do
    (cd "low_en_structures/basin_${n}" && cat conf_*.gro > "selected_basin_${n}.gro")
done
