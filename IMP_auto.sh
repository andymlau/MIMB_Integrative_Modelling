#!/bin/bash

alias python='/data/users/shared/software/anaconda/anaconda2/bin/python'

#module load imp/2.7.0
module load pymol

python modeling.py

mkdir ./output_atomistic;

pymol -c ./scripts/CG2pdb.py

python clustering.py

mkdir ./kmeans_100_5/atomistic_cluster.0;
mkdir ./kmeans_100_5/atomistic_cluster.1;
mkdir ./kmeans_100_5/atomistic_cluster.2;
mkdir ./kmeans_100_5/atomistic_cluster.3;
mkdir ./kmeans_100_5/atomistic_cluster.4;

pymol -c ./scripts/CG2pdb_CG0.py
pymol -c ./scripts/CG2pdb_CG1.py
pymol -c ./scripts/CG2pdb_CG2.py
pymol -c ./scripts/CG2pdb_CG3.py
pymol -c ./scripts/CG2pdb_CG4.py

#rm -r ./kmeans_100_5/cluster*
