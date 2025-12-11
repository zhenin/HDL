load("Figure4.rda")

nejm.colors <- function (n = 8, alpha = 1){
  pal <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
           "#6F99AD", "#FFDC91", "#EE4C97")
  acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
  return(paste0(pal, acode)[1:n])
}

swap_allele <- function(x) {
  require(dplyr)
  
  x.1 <- x[x$Z >= 0,]
  x.2 <- x[x$Z < 0,]
  
  x.2  <- x.2 %>% rename(A1 = A2, A2 = A1) %>% mutate(Z = -Z)
  x.res <- rbind(x.1,x.2)
  x.res <- x.res[x$SNP,on = .(SNP)] %>% na.omit
  
  return(x.res)
}

plot_Z <- function(data, cate, leadsnp, xlab = TRUE, ylab = TRUE, lege, nosig = FALSE,sigp = 0.05, sigp1 = 5e-8, 
                   xmin = NA, xmax = NA, ymin = NA, ymax = NA,y_breaks = NULL, 
                   abline = TRUE, plotmar = margin(1, 5.5, 4, 5.5, "pt")) {
  
  require(ggplot2)
  # require(scales)
  p1 <- ggplot() +
    theme_classic() +
    labs(title = cate) +
    geom_point(data = data[SNP != leadsnp,],shape = 16,aes(x = position,y = Z,col = r2)) +
    geom_point(data = data[SNP == leadsnp,],shape = 18,aes(x = position,y = Z),col = nejm.colors()[2],size = 3) +
    scale_colour_gradient(high = nejm.colors()[1],low = nejm.colors()[6], limits = c(0,1)) +
    # scale_y_continuous(breaks= pretty_breaks()) +
    theme(title = element_text(size = 9,margin = margin()),
          legend.key.height = unit(.55, "cm"),
          legend.key.width = unit(.25, "cm"),
          legend.margin = margin(),
          legend.text = element_text(size = 7, margin = margin(l=2, unit = "pt")),
          legend.title = element_text(size = 9,margin = margin(b = 8, unit = "pt")),
          plot.title = element_text(margin = margin(t = 0,b = 2)),
          plot.margin = plotmar) 
  
  if (abline) {
    p1 <- p1 + geom_abline(slope = 0,intercept = 0,lty = 2,col = 'gray30')
  }
  
  if (!nosig) {
    sigZ <- sqrt(qchisq(sigp, df = 1, low = F))
    sigZ1 <- sqrt(qchisq(sigp1, df = 1, low = F))
    p1 <- p1 +geom_hline(yintercept = c(sigZ, -sigZ), lty=2, col='black') + 
      geom_hline(yintercept = c(sigZ1, -sigZ1), lty=2, col='black')
  }
  
  if (xlab == FALSE) {
    p1 <- p1 + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.line.x = element_blank())
  } else {
    p1 <- p1 + theme(axis.text.x = element_text(color = 'black'),
                     axis.ticks.x = element_line(color = 'black'))
  }
  
  if (ylab == FALSE) {
    p1 <- p1 + theme(axis.title.y = element_blank(),
                     axis.text.y = element_text(color = 'black',size = 11.5),
                     axis.ticks.y = element_line(color = 'black'))
  } else if (ylab == TRUE) {
    p1 <- p1 + theme(axis.title.y = element_text(angle = 0, vjust = 0.5,size = 14),
                     axis.text.y = element_text(color = 'black',size = 10),
                     axis.ticks.y = element_line(color = 'black'))
  }
  
  if (lege == FALSE) {
    p1 <- p1 + guides(col = 'none')
  } else {
    p1 <- p1 + labs(col = expression('LD R'^'2'))
    #p1 <- p1 + labs(col = expression('LD abs(R)'))
  }
  
  if (!is.null(y_breaks)) {
    p1 <- p1 + scale_y_continuous(breaks = y_breaks)
  }
  
  p1 <- p1 + coord_cartesian(xlim = c(xmin, xmax),ylim = c(ymin, ymax))
  
  return(p1)
}


grep_genes <- function(chr, region.start, region.end, gene_name, GRCh = 37) {
  require(dplyr)
  require(biomaRt)
  ## GRCh37 bp or GRCh38 bp
  if (GRCh == 37) {
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = GRCh)
  } else {
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  }
  genes <- getBM(
    attributes = c("external_gene_name", "ensembl_gene_id_version", "gene_biotype", "chromosome_name", "start_position", "end_position","strand"),
    mart = ensembl,
    filters = c("chromosome_name", "start", "end"), values = list(c(chr), region.start, region.end)
  )
  genes$draw_start_position <- ifelse(genes$strand > 0,genes$start_position,genes$end_position)
  genes$draw_end_position <- ifelse(genes$strand > 0,genes$end_position,genes$start_position)
  if (nrow(genes[genes$gene_biotype == "protein_coding", ]) == 0) {
    genes <- NULL
  } else if (!is.na(gene_name)) {
    genes <- genes[genes$gene_biotype == "protein_coding", ] %>% 
      filter(start_position > region.start, end_position < region.end) %>% 
      mutate(order = sample(1:nrow(.)),col = ifelse(external_gene_name == gene_name,'red','black'))
  } else {
    genes <- genes[genes$gene_biotype == "protein_coding", ] %>% 
      filter(start_position > region.start, end_position < region.end) %>% 
      mutate(order = sample(1:nrow(.)),col = 'black')
  }
  return(genes)
  
}

plot_gene <- function(genes = NULL,ylab = FALSE,xlab = TRUE,xlim_l = NULL,xlim_r = NULL,gene_size = 1,gap_size = 0.1,text_size = 20,arrow_size = 5,gap=2) {
  
  require(ggplot2)
  gene_num <- nrow(genes)
  minpos <- -gene_num
  maxpos <- gene_num
  #gap <- 2
  dis = abs(genes$draw_end_position - genes$draw_start_position)
  
  if (!is.null(genes)) {
    p <- ggplot() + theme_classic() + 
      #  geom_segment(data = genes, 
      #               aes(x = draw_start_position, xend = draw_end_position, y = maxpos - gap*order + gap_size/0.8 + 1, yend = maxpos - gap*order + gap_size/0.8 + 1),
      #               arrow = arrow(length = unit(arrow_size, "mm")),color = genes$col,size = gene_size*4) +
      #  geom_text(data = genes, aes(x = draw_end_position + (region.end - region.start)*0.15,
      #                              y = maxpos - gap*order*0.9 + 1,label = external_gene_name,fontface=3),
      #            size = gene_size*text_size,colour = genes$col) +
      geom_segment(data = genes, 
                   aes(x = draw_start_position, xend = draw_end_position, y = maxpos - gap*order + 1, yend = maxpos - gap*order + 1),
                   arrow = arrow(length = unit(arrow_size, "mm")),color = genes$col,size = gene_size*4) +
      geom_text(data = genes, aes(x = draw_start_position + (draw_end_position - draw_start_position)*4,
                                  y = maxpos - gap*order + 1,label = external_gene_name,fontface=3),
                size = gene_size*text_size,colour = genes$col) +
      theme(plot.margin = margin(0, 5.5, 5.5, 5.5, "pt")) +
      coord_cartesian(ylim = c(minpos, maxpos))
  } else {
    p <- NULL
    return(p)
  }
  
  if (any(c(!is.null(xlim_l),!is.null(xlim_r)))) {
    p <- p + coord_cartesian(xlim = c(xlim_l,xlim_r),ylim = c(minpos, maxpos))
  }
  
  if (ylab == FALSE) {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line.y = element_blank())
  } 
  
  if (xlab == FALSE) {
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.line.x = element_blank())
  } else if(xlab == TRUE) {
    p <- p + xlab(paste0('Chromosome ',unique(genes$chromosome_name),' position (Mbp)')) + 
      theme(axis.text.x = element_text(color = 'black',size = 12),
            axis.ticks.x = element_line(color = 'black'),
            axis.title.x = element_text(vjust = -.1,size = 13)) + 
      scale_x_continuous(labels = scales::number_format(scale = 1e-6,accuracy = 0.1))
  }
  
  return(p)
}



fld1name = "CELSR2"
fld2name = "LDL"
fld3name = "I25"



p1 <- plot_Z(data = CELSR2sub,cate = fld1name, leadsnp = leadsnp2, ylab = TRUE,nosig = TRUE,xlab = FALSE,
             lege = FALSE,xmin = region.start,xmax = region.end)

p2 <- plot_Z(data = LDLsub,cate = fld2name, leadsnp = leadsnp2, ylab = TRUE,nosig = TRUE,xlab = FALSE,
             lege = FALSE,xmin = region.start,xmax = region.end)

p3 <- plot_Z(data = I25sub,cate = fld3name,leadsnp = leadsnp2, xlab = FALSE,ylab = TRUE,
             nosig = TRUE,lege = TRUE,xmin = region.start,xmax = region.end)
require(patchwork)
layout <- "
A
A
B
B
C
C
"
p <- p1 + p2 + p3 + plot_layout(design = layout) 
p
