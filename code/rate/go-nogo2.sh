#! /bin/sh
#
# test.sh
# Copyright (C) 2019 rkim <rkim@salk.edu>
#
# Distributed under terms of the MIT license.
#

for i in {1..5}; do \
  python main_sequential.py --gpu 0 --gpu_frac 0.20 --n_trials 10000 --mode train \
  --N 100 --P_inh 0.20 --som_N 0 --task go-nogo2 --gain 1.5 --P_rec 0.20 \
  --act sigmoid --loss_fn l2 --apply_dale True --decay_taus 4 25 --prob 0.20 \
  --output_dir ../../ 
done


