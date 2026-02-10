#!/usr/bin/env bash
set -euo pipefail

source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh
conda activate r

ROOT=/work/sph-zhaor/analysis/ckm
cd ${ROOT}

mkdir -p ${ROOT}/code ${ROOT}/jobs/r ${ROOT}/jobs/sh ${ROOT}/jobs/log \
         ${ROOT}/data/fg_prot ${ROOT}/data/fg_met ${ROOT}/res/tmp_fg

Rscript ${ROOT}/code/fg_setup.R
bash ${ROOT}/jobs/submit_all.sh
