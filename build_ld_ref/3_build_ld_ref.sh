#!/usr/bin/env bash

outprefix=$1
if [ $# -ne 1 ];then
  echo "bash $0 <ld_ref_path/ld_ref_name>"
  exit 1
fi

src_dir=`dirname $0`

for ld_file in `ls ${outprefix}*.ld`; do
  echo "Rscript ${src_dir}/utils/eigen.R ${ld_file}";
done
