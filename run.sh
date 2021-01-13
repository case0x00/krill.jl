#!/bin/bash

mu=0.012277471


# CR3BP
# pre L1-keyhole
#python3 core/krill.py $MU -0.75 0.0 0.0 0.0
# L1 keyhole
python3 core/krill.py $mu 0.309 0.077 -0.726 1.555


# twobody:
#python3 core/krill.py $MU 0 41797000 3070 0
