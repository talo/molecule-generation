#!/bin/bash

time_hash=$(date +%s)
python3 /home/minhtoo/molecule-generation/optimise/run.py \
    2>&1 | tee /home/minhtoo/molecule-generation/logs/100parts_10steps_5res_5poses_8exh_v1_top100PRTs_${time_hash}.log