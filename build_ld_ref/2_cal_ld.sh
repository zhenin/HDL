#!/usr/bin/env bash

bfile=$1
outprefix=$2
bandwidth=${3:-500}      # bandwith (number of SNPs) for LD calculation, default=500
ld_window=${4:-100000}   # window size (kb) for LD calculation, default=1000000 (whole segment)
chroms=${5:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}

if [ $# -lt 2 ];then
  echo "bash $0 <path/to/bfile> <ld_ref_path/ld_ref_name> [bandwidth [ld_window [chroms]]]"
  exit 1
fi

src_dir=`dirname $0`

case `uname -s` in
  Linux*)     os=linux;;
  Darwin*)    os=macos;;
  *)          os=linux
esac

plink=${src_dir}/utils/plink_1.9_${os}

set -f
chroms=(${chroms//,/ })

cmd_prefix="Rscript ${src_dir}/utils/cal_ld.R --bandwidth ${bandwidth} --ld-window ${ld_window} --print-cmd"

for chrom in ${chroms[@]}; do
  ${cmd_prefix} ${bfile} ${outprefix} ${chrom} ${plink};
done
