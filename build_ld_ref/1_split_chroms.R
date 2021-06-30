#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(require(dplyr)))
suppressWarnings(suppressMessages(require(tidyr)))
suppressWarnings(suppressMessages(require(data.table)))
suppressWarnings(suppressMessages(require(argparser)))

piece.length <- function(len, mean.len = 5000) {
  n.chr <- length(len)
  n.piece <- rep(NA, n.chr)
  piece.len <- c()
  piece.end <- c()
  for (i in 1:n.chr) {
    n.piece[i] <- max(1, round(len[i]/mean.len))
    len.i <- floor(len[i]/n.piece[i])
    piece.len <- c(piece.len, rep(len.i, n.piece[i] - 1), len[i] - len.i*(n.piece[i] - 1))
    piece.end <- c(piece.end, seq(len.i, len[i], len.i)[-n.piece[i]], len[i])
  }
  sL <- var(piece.len)
  mL <- mean(piece.len)
  return(list(cv = sL/mL, mean.length = mean.len, piece.len = piece.len, piece.ends = piece.end))
}

optim.length <- function(len, min.length = 5000, max.length = 15000) {
  cv <- 1e6
  cvs <- c()
  for (mlen in max.length:min.length) {
    pl <- piece.length(len, mean.len = mlen)
    if (pl$cv < cv) {
      out <- pl
      cv <- pl$cv
      cvs <- c(cvs, cv)
    }
  }
  out$cvs <- cvs
  return(out)
}

drop.region <- function(df, r){
  region <- rev(strsplit(r, '')[[1]])
  end <- which(region=='-')[1]
  start <- which(region==':')[1]
  cur.chrom <- paste0(region[length(region):(start+1)], collapse='')
  start <- as.numeric(paste0(region[(start-1):(end+1)], collapse=''))
  end <- as.numeric(paste0(region[(end-1):1], collapse=''))
  if(cur.chrom=='' | is.na(cur.chrom) | is.na(start) | is.na(end) | start > end){
    stop(paste0('Invalid chromosome region: ', r))
  }
  
  df.sub1 <- filter(df, chrom==cur.chrom, pos<start)
  df.sub2 <- filter(df, chrom==cur.chrom, pos>end)
  if(nrow(df.sub1)>0 & nrow(df.sub2) > 0){
    seg.grp1 <- tail(df.sub1, 1)$seg.grp
    seg.grp2 <- head(df.sub2, 1)$seg.grp
    if(seg.grp1 != 0 & seg.grp1 == seg.grp2){
      max.seg.grp <- max(c(df.sub1$seg.grp, df.sub2$seg.grp)) + 1
      df.sub2 <- mutate(df.sub2, seg.grp=ifelse(seg.grp==seg.grp2, max.seg.grp, seg.grp))
    }
  }
  
  df <- df %>%
    filter(chrom != cur.chrom) %>%
    rbind(., rbind(df.sub1, df.sub2))
  return(df)
}

split.bim <- function(outprefix, bim, min.l=1000, max.l=3000, regions=NULL,
                      chrom.ord=as.character(c(1:100, 'X', 'Y', 'MT'))){
  dup.rsids <- bim %>%
    count(rsid, name='count') %>%
    filter(count>1) %>%
    select(rsid) %>%
    pull()

  if(length(dup.rsids) >= 1){
    dup.rsids <- paste0(dup.rsids, collapse=', ')
    msg <- paste0('Variant identifiers in the bim file MUST BE UNIQUE. Please check variants:\n',
                  dup.rsids, '\n\n')
    stop(msg)
  }

  df <- mutate(bim, chrom=as.character(chrom), seg.grp=1)
  if(!is.null(regions)){
    regions <- strsplit(regions, ';')[[1]]
    for(region in regions){
      df <- drop.region(df, region)
    }
  }

  df <- df %>%
    mutate(chrom=factor(chrom, levels=chrom.ord)) %>%
    arrange(chrom, pos) %>%
    add_count(chrom, seg.grp, name='count')

  len <- df %>%
    distinct(chrom, seg.grp, count) %>%
    select(count) %>%
    pull()

  m.len <- sum(len) / 50
  min.length <- round(m.len * 0.8)
  max.length <- min(round(m.len * 1.2), 20000)
  if(min.length > max.length) min.length <- round(max.length/1.2*0.8, 0)

  m.len <- (max.length + min.length) / 2
  
  if(m.len < min.l | min.l > 20000){
    msg <- paste0('We recommend max <= 20000 for eigen decomposition efficiency,',
                  ' and roughly splitting all chromosomes into 50 segments or more for the standard error estimation in the HDL.')
    warning(msg)
  }

  if(max.l <= min.length | min.l >= max.length){
    min.length <- min.l
    max.length <- max.l
  }else{
    min.length <- max(min.l, min.length)
    max.length <- min(max.l, max.length)
  }
  
  optim.res <- optim.length(len, min.length, max.length)
  
  df <- data.frame(e=optim.res$piece.len) %>%
    mutate(e=cumsum(e),
           s=lag(e),
           s=ifelse(is.na(s), 1, s+1)) %>%
    left_join(mutate(df, e=1:n()), ., by='e') %>%
    select(-seg.grp, -count) %>%
    
    arrange(desc(e)) %>%
    mutate(i=cumsum(!is.na(s))) %>%
    arrange(chrom, pos) %>%
    group_by(chrom, i) %>%
    mutate(s=min(pos),
           e=max(pos),
           len=n()) %>%
    ungroup() %>%
    select(chrom, rsid, pos, s, e, len) %>%
    arrange(chrom, pos)
  
  snps.name.list <- pull(df, rsid)
  
  df <- distinct(df, chrom, s, e, len)
  cat('Length of all segments ( average =',
      mean(df$len), '& standard deviation =', sd(df$len), '): \n\n')
  cat(df$len, '\n\n')
  chroms <- unique(df$chrom)
  nsnps.list <- list()
  for (cur.chrom in chroms){
    nsnps.list[[cur.chrom]] <- df %>%
      filter(chrom==cur.chrom) %>%
      select(len) %>%
      pull()
  }
  
  out.dir <- dirname(outprefix)
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  save(nsnps.list, file=paste0(outprefix, '_snp_counter_array.RData'))
  save(snps.name.list, file=paste0(outprefix, '_snp_list_array.RData'))
}

args <- arg_parser('Split chromosomes into segments.') %>%
  add_argument('ldprefix', help='ld_ref_path/ld_ref_name', type='character') %>%
  add_argument('bim', help='path to .bim file of ALL chromsomes to be included in your LD reference panel.', type='character') %>%
  add_argument('--min', help='min average number of SNPs in a segment', type='numeric', default=3000) %>%
  add_argument('--max', help='max average number of SNPs in a segment', type='numeric', default=20000) %>%
  add_argument('--exclude', help='exclude chromosome regions separated by `;` and quoted with `"`, e.g. "chromA:startA-endA;chromB:startB-endB;chromC:startC-endC;..."        or set as "keep" to avoid excluding', type='character', default='"6:25643792-33429951"') %>%
  parse_args()

ldprefix <- args$ldprefix
bim.file <- args$bim
min.l <- args$min
max.l <- args$max
regions <- gsub('^"|"$', '', args$exclude)
if(regions=='keep') regions <- NULL

if(min.l > max.l) stop('value of --min option is greater than value of --max option.')

bim <- fread(bim.file, select=c(1, 2, 4), col.names=c('chrom', 'rsid', 'pos'))
split.bim(ldprefix, bim, min.l, max.l, regions)