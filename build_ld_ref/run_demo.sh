#!/usr/bin/env bash

build_ref_dir=`dirname $0`

bfile=${build_ref_dir}/demo/test
ld_ref_prefix=${build_ref_dir}/demo/ld_ref/test

# build LD ref
Rscript ${build_ref_dir}/1_split_chroms.R ${ld_ref_prefix} ${bfile}.bim --min 100 --max 500
bash ${build_ref_dir}/2_cal_ld.sh ${bfile} ${ld_ref_prefix} 50 1000000 | bash
bash ${build_ref_dir}/3_build_ld_ref.sh ${ld_ref_prefix} | bash
rm ${ld_ref_prefix}*.log ${ld_ref_prefix}*.nosex ${ld_ref_prefix}*.ld

# test LD ref
LD_path=`dirname $ld_ref_prefix`
Rscript -e "library(HDL);data(gwas1.example);data(gwas2.example);LD.path='${LD_path}';HDL.rg(gwas1.example, gwas2.example, LD.path, Nref=504, eigen.cut=0.01)"
