#!/bin/bash
# 立体异构体枚举运行脚本 - 多进程并行版本

cd /inspire/hdd/project/chemicalreaction/misixuan-CZXS24220243/qb/github/2daam2d/depth0_sample1/enumerate_steoro

python enumerate_stereo_parallel.py \
    /inspire/hdd/project/chemicalreaction/misixuan-CZXS24220243/qb/github/2daam2d/depth0_sample1/sampled_depth0_sample1_aug.csv \
    /inspire/hdd/project/chemicalreaction/misixuan-CZXS24220243/qb/github/2daam2d/depth0_sample1/enumerate_steoro/sampled_depth0_sample1_aug_stereo.csv \
    --max-isomers 128 \
    --num-workers 64 \
    --debug
