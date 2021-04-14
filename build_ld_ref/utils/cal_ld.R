#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(require(dplyr)))
suppressWarnings(suppressMessages(require(tidyr)))
suppressWarnings(suppressMessages(require(argparser)))
suppressWarnings(suppressMessages(require(data.table)))

ld.cmd <- function(bfile, ldprefix, chrom, i.segment, from.snp, to.snp,
                   bandwidth=500, ld.window.kb=100000, plink.path='plink'){
  paste0(plink.path, ' --silent --bfile ', bfile, ' --r --ld-window ', bandwidth+1,
         ' --ld-window-kb ', as.integer(ld.window.kb, 0),
         ' --from ', from.snp, ' --to ', to.snp,
         ' --out ', ldprefix, '_chr', chrom, '.', i.segment, '_')
}

write.bim <- function(chrom, i.segment, from.snp, to.snp, ldprefix, bim){
  from.pos <- filter(bim, snp==from.snp)$pos
  to.pos <- filter(bim, snp==to.snp)$pos
  bim %>%
    filter(pos<=to.pos, pos>=from.pos) %>%
    write.table(., file=paste0(ldprefix, '_chr', chrom, '.', i.segment, '_.bim'),
                sep=' ', quote=F, row.names=F, col.names=F)
}

load.seg.split.bim <- function(ldprefix, chrom, bim.file){
  nsnps.list.file <- paste0(ldprefix, '_snp_counter_array.RData')
  snps.name.list.file <- paste0(ldprefix, '_snp_list_array.RData')
  if(!file.exists(nsnps.list.file)){
    stop(paste0('[file missing] : ', nsnps.list.file,
                '. Please check the prefix of LD reference or run the first (chromosomes spliting) step correctly.'))
  }
  if(!file.exists(snps.name.list.file)){
    stop(paste0('[file missing] : ', snps.name.list.file,
                '. Please check the prefix of LD reference or run the first (chromosomes spliting) step correctly.'))
  }
  load(nsnps.list.file)
  load(snps.name.list.file)
  
  idx <- which(names(nsnps.list)==chrom)
  if(idx == 0){
    stop(paste0('[chromosome missing] : ', chrom,
                '. You should `cat` all the required .bim files into a single .bim file, and re-run the first step.'))
  }
  
  bim <- fread(bim.file, header=F,
               col.names=c('chr', 'snp', 'gpos',
                           'pos', 'ref', 'alt')) %>%
    mutate(chr=as.character(chr)) %>%
    filter(chr==chrom)
  
  e <- do.call('sum', nsnps.list[1:idx])
  s <- e - sum(nsnps.list[[chrom]]) + 1
  snps <- snps.name.list[s:e]
  
  from.snps <- c()
  to.snps <- c()
  s <- 1
  segs <- nsnps.list[[chrom]]
  for(i.segment in 1:length(segs)){
    e <- s + segs[i.segment] - 1
    from.snps <- c(from.snps, snps[s])
    to.snps <- c(to.snps, snps[e])
    write.bim(chrom, i.segment, snps[s], snps[e], ldprefix, bim) # write .bim
    s <- e + 1
  }
  return(data.frame(chrom=chrom,
                    from.snp=from.snps,
                    to.snp=to.snps,
                    nsnps=nsnps.list[[chrom]],
                    i.segment=1:length(from.snps)))
}

args <- arg_parser('Calculate LD.') %>%
  add_argument('bfile', help='path to the plink bfile', type='character') %>%
  add_argument('ldprefix', help='ld_ref_path/ld_ref_name', type='character') %>%
  add_argument('chrom', help='chromsome', type='character') %>%
  add_argument('plink_path', help='path to plink software', type='character') %>%
  add_argument('--print-cmd', help='print plink ld calculation commands', flag=TRUE) %>%
  add_argument('--bandwidth', help='bandwith (# of SNPs) for LD calculation, default=500',
               type='numeric', default=500) %>%
  add_argument('--ld-window',
               help='window size (kb) for LD calculation, default=1000000 (whole segment)', 
               type='numeric', default=1000000) %>%
  parse_args()


bfile <- args$bfile
ldprefix <- args$ldprefix
chrom <- args$chrom
plink.path <- args$plink_path
bandwidth <- round(args$bandwidth, 0)
ld.window.kb <- round(args$ld_window, 0)
print.cmd <- args$print_cmd


bim.file <- paste0(bfile, '.bim')

# load segment info & split bim
seg.info <- load.seg.split.bim(ldprefix, chrom, paste0(bfile, '.bim'))


# calculate LD
cmds <- seg.info %>%
  rowwise(.) %>%
  mutate(cmd=ld.cmd(bfile, ldprefix, chrom, i.segment, from.snp, to.snp, bandwidth, ld.window.kb, plink.path)) %>%
  ungroup() %>%
  # mutate(cmd=paste0(cmd, ' ; gzip -f ', ldprefix, '.ld')) %>%
  select(cmd) %>%
  pull()

if(print.cmd){
  # print ld cmd
  cat(paste0(paste0(cmds, collapse='\n'), '\n'))
}else{
  # run ld cmd
  for(cmd in cmds){
    system(cmd)
  }
}

