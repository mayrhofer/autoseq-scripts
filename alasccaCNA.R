#!/usr/bin/env Rscript

###########################################
# check command line
library(getopt)

#long, short, argmask, datatype, desk
#argmask 0=no arg, 1=req, 2=optional
args <- rbind(c("cnr", "b", 1, "character", "bin file from CNVkit"),
	      c("cns", "s", 1, "character", "segment file from CNVkit"),
	      c("germlinevcf", "g", 1, "character", "germline variant vcf file"),
	      c("somaticvcf", "m", 1, "character", "somatic variant vcf file"),
	      c("png", "p", 1, "character", "output png file"),
	      c("json", "o", 1, "character", "output json file"),
	      #c("genes.genePred", "g", 1, "character", "ensemble genePred file"),
	      c("chrsizes", "c", 1, "character", "chromosome sizes file"))


opts <- getopt(args)

if(is.null(opts$cnr)){
  stop("Must specify input bin file --cnr/-b.")
}
if(is.null(opts$cns)){
  stop("Must specify input segment file --cns/-s.")
}
if(is.null(opts$germlinevcf)){
  stop("Must specify germline vcf file --germlinevcf/-g.")
}
if(is.null(opts$somaticvcf)){
  stop("Must specify somatic vcf file --somaticvcf/-m.")
}
if(is.null(opts$png)){
  stop("Must specify output png file name --png/-p.")
}
if(is.null(opts$json)){
  stop("Must specify output json file name --json/-o")
}
if(is.null(opts$chrsizes)){
  stop("Must specify chrsizes file --chrsizes/-c.")
}

##############################

#Plot settings
colors_p <- colorRampPalette(c("#6600FF","#9900CC"),space="rgb")
colors_q <- colorRampPalette(c("#CC0099","#CC0000"),space="rgb")

# setup
library(data.table)
library(VariantAnnotation)
library(RJSONIO)

chrsizes <- fread(opts$chrsizes)
setnames(chrsizes, names(chrsizes), c("chr", "size"))
chrsizes$cumend <- cumsum( as.numeric( chrsizes$size) )
chrsizes$cumstart <- as.numeric(chrsizes$cumend) - as.numeric(chrsizes$size)
chrsizes$labelpos <- cumsum( as.numeric(chrsizes$size)) - chrsizes$size/2

mids_chrnames <- c(1:22, 'X', 'Y', 'MT')
mids_list <- c(1.25e+08, 93300000, 9.1e+07, 50400000, 48400000, 6.1e+07, 59900000,
               45600000, 4.9e+07, 40200000, 53700000, 35800000, 17900000, 17600000,
               1.9e+07, 36600000, 2.4e+07, 17200000, 26500000, 27500000, 13200000,
               14700000, 60600000, 12500000,NA)

chrsizes$mid <- mids_list[match(chrsizes$chr, mids_chrnames)]

g <- data.frame(label=c("ARID1A", "NRAS", "PIK3CA", "APC", "PIK3R1", "BRAF", "MYC", "AR", "PTEN","IGF2", "KRAS", "IRS2", "TP53", "PMS2", "BRCA1", "BRCA2", "MLH1"),
		chromosome=NA,start=NA,end=NA,stringsAsFactors = F)

genes=structure(list(label = c("ARID1A", "NRAS", "PIK3CA", "APC", "PIK3R1",
			       "BRAF", "MYC", "AR", "PTEN", "IGF2", "KRAS", "IRS2", "TP53",
			       "PMS2", "BRCA1", "BRCA2", "MLH1"), chromosome = c("1", "1", "3",
										 "5", "5", "7", "8", "X", "10", "11", "12", "13", "17", "7", "17",
										 "13", "3"), start = c(27022521L, 115247084L, 178866310L, 112043201L,
												       67584251L, 140433812L, 128748314L, 66763873L, 89623194L, 2150346L,
												       25358179L, 110406183L, 7571719L, 6012869L, 41196311L, 32889616L,
												       37034840L), end = c(27108601L, 115259515L, 178952497L, 112181936L,
															   67597649L, 140624564L, 128753680L, 66950461L, 89728532L, 2160204L,
															   25403854L, 110438914L, 7590868L, 6048737L, 41277468L, 32973809L,
															   37092337L)), .Names = c("label", "chromosome", "start", "end"
															   ), row.names = c(NA, -17L), class = "data.frame")
genes=genes[genes$label %in% c('BRAF','KRAS','NRAS','PIK3CA','PIK3R1','PTEN','IGF2'),]

genes$cumstart <- genes$start + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]
genes$cumend <- genes$end + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]
#genes=genes[genes$label %in% c('PTEN','IGF2'),]
pten=structure(list(start = c(89623194, 89623861, 89653781, 89685269,
			      89690802, 89692769, 89711874, 89717609, 89720650, 89725043),
		    end = c(89623860, 89624305, 89653866, 89685314, 89690846,
			    89693008, 89712016, 89717776, 89720875, 89731687)),
	       .Names = c("start", "end"),
	       row.names = c(NA, -10L), class = "data.frame")
igf2=structure(list(start = c(2170355, 2168795, 2156596, 2154746, 2150341),
		    end = c(2170833, 2169037, 2156759, 2154895, 2154453)),
	       .Names = c("start", "end"),
	       class = "data.frame", row.names = c(NA,-5L))








##Below: for looping.
# cns=sort(dir()[grep(pattern = '.cns',dir())])
# cnr=sort(dir()[grep(pattern = '.cnr',dir())])
# vcfgz=sort(dir()[grep(pattern = '.vcf.gz',dir())])
# for (sample in 1:length(cnr)) {
#   opts <- getopt(args,opt = c("--cnr", paste(getwd(),cnr[sample],sep='/'),
#                               "--cns", paste(getwd(),cns[sample],sep='/'),
#                               "--germlinevcf", paste(getwd(),vcfgz[sample*2-1],sep='/'),
#                               "--somaticvcf", paste(getwd(),vcfgz[sample*2],sep='/'),
#                               "--png", paste(sample,"plot.png",sep='.'),
#                               "--json", paste(sample,"json",sep='.'),
#                               #"--genelist", "ensembl_genes_v75_cleaned.genePred",
#                               "--chrsizes", "human_g1k_v37_decoy.chrsizes.txt") )
# comment this far to disable loop

segments <- fread(opts$cns,stringsAsFactors = F)
bins <- fread(opts$cnr,stringsAsFactors = F)

## Segment start and end pos
segments$cumstart <- NA
segments$cumend   <- NA
for(chr in chrsizes$chr){
  idx <- which( segments$chr == chr )
  cumchrsize <- chrsizes$cumend[which(chrsizes$chr == chr)] -
    chrsizes$size[which(chrsizes$chr == chr)]
  segments$cumstart[idx] <- segments$start[idx] + cumchrsize
  segments$cumend[idx] <- segments$end[idx] + cumchrsize
}
segments$centerpos <- segments$cumstart+(segments$cumend-segments$cumstart)/2

## Bin start and end pos
bins$cumstart <- NA
bins$cumend   <- NA
for(chr in chrsizes$chr){
  idx <- which( bins$chr == chr )
  cumchrsize <- chrsizes$cumend[which(chrsizes$chr == chr)] -
    chrsizes$size[which(chrsizes$chr == chr)]
  bins$cumstart[idx] <- bins$start[idx] + cumchrsize
  bins$cumend[idx] <- bins$end[idx] + cumchrsize
}
bins$centerpos <- bins$cumstart+(bins$cumend-bins$cumstart)/2

## Antitarget bins can be awfully low. Reduce their weight:
ix=bins$gene=='Background'
low=bins$log2 < -5
bins$weight[ix&low] <- 0.01


# Segment value correction
bins$segmented <- NA
for (i in 1:nrow(segments)) {
  ix=bins$chromosome==segments$chromosome[i] & bins$start>=segments$start[i] & bins$end<=segments$end[i]
  s=weighted.mean(x = bins$log2[ix],w = bins$weight[ix])
  bins$segmented[ix] <- segments$log2[i] <- s
}
# Baseline correction
baseline <- median(rep(segments$log2,round(segments$weight)))
segments$log2 <- segments$log2-baseline
bins$log2 <- bins$log2-baseline

# Smoothed bin values
bins$smoothed=NA
kernel <- c(dnorm(seq(0, 3, length.out=40)))
kernel <- c(kernel, rep(0, nrow(bins) - 2*length(kernel) + 1), rev(kernel[-1]))
kernel <- kernel / sum(kernel)
bins$smoothed <- Re(convolve(bins$log2, kernel))
bins$w=10
bins$w[bins$gene=='Background']=1
bins$pch=1
bins$pch[bins$gene=='Background']=16
bins$col='black'
bins$col[bins$gene=='Background']='grey'

vcf <- readVcf(opts$germlinevcf,genome = "GRCh37")
g <- geno(vcf)
chr <- pos <- rownames(g$DP)

for (i in 1:length(pos)) {
  temp <- strsplit(chr[i],':')[[1]]
  chr[i] <- as.character(temp[1])
  pos[i] <- as.numeric(strsplit(temp[2],'_')[[1]][1])
}
pos=as.numeric(pos)
alf <- data.frame(chromosome=chr, start=pos, end=pos, stringsAsFactors = F, cumstart=NA, cumend=NA)

if( ! all(is.na(alf$chromosome)) ) {
  
  for(chr in chrsizes$chr){
    idx <- which(alf$chromosome == chr)
    cumchrsize <- chrsizes$cumend[which(chrsizes$chr == chr)] -
      chrsizes$size[which(chrsizes$chr == chr)]
    alf$cumstart[idx] <- alf$start[idx] + cumchrsize
    alf$cumend[idx] <- alf$end[idx] + cumchrsize
  }
  
  alf$t <- as.numeric(g$AO[,2])/as.numeric(g$DP[,2])
  alf$n <- as.numeric(g$AO[,1])/as.numeric(g$DP[,1])
  alf$td <- as.numeric(g$DP[,2])
  alf$log2 <- NA
  delta <- 2e6; for (i in 1:nrow(alf)) {
    ix <- bins$chromosome==alf$chromosome[i] & bins$start>alf$start[i]-delta & bins$end<alf$start[i]+delta
    alf$log2[i]=median(bins$log2[ix],na.rm=T)
  }
  alf$ai=2*abs(alf$t-0.5)
  alf$col=NA
  for (i in 1:24) {
    ix=which(alf$chromosome==chrsizes$chr[i] & alf$start<chrsizes$mid[i])
    alf$col[ix]=colors_p(length(ix))
    ix=which(alf$chromosome==chrsizes$chr[i] & alf$start>chrsizes$mid[i])
    alf$col[ix]=colors_q(length(ix))
  }
  segments$ai <- NA; for (i in 1:nrow(segments)) {
    ix=alf$chromosome==segments$chromosome[i] & alf$start>=segments$start[i] & alf$end<=segments$end[i]
    if (sum(alf$td[ix],na.rm=T)>500) segments$ai[i] <- median(alf$ai[ix])#,w = alf$td[ix])
    #alf$ai[ix] <- median(alf$ai[ix])
  }
}


vcf2 <- readVcf(opts$somaticvcf,genome = "GRCh37")
g <- geno(vcf2)
chr <- pos <- rownames(g$DP)
for (i in 1:length(pos)) {
  temp <- strsplit(chr[i],':')[[1]]
  chr[i] <- as.character(temp[1])
  pos[i] <- as.numeric(strsplit(temp[2],'_')[[1]][1])
}
pos=as.numeric(pos)
salf <- data.frame(chromosome=chr,start=pos,end=pos,stringsAsFactors = F)
salf$cumstart <- NA
salf$cumend   <- NA
for(chr in chrsizes$chr){
  idx <- which(salf$chromosome == chr)
  cumchrsize <- chrsizes$cumend[which(chrsizes$chr == chr)] -
    chrsizes$size[which(chrsizes$chr == chr)]
  salf$cumstart[idx] <- salf$start[idx] + cumchrsize
  salf$cumend[idx] <- salf$end[idx] + cumchrsize
}

salf$t <- as.numeric(g$VD[,1])/as.numeric(g$DP[,1])
salf$t.altreads <- as.numeric(g$VD[,1])
salf$n <- as.numeric(g$VD[,2])/as.numeric(g$DP[,2])
salf$n.totreads <- as.numeric(g$DP[,2])

salf$type='other'
salf$type[isSNV(vcf2)]='snv'
salf$type[isDeletion(vcf2)]='del'
salf$type[isInsertion(vcf2)]='ins'

salf$gene=''
salf$aa=''
vep=info(vcf2)$CSQ
table=NULL

alasccaTX=c('ENST00000288602','ENST00000256078','ENST00000369535',
	    'ENST00000371953','ENST00000521381','ENST00000263967')

for (i in 1:length(vep)) try( {
  temp=NULL
  for (j in 1:length(vep[[i]])) {
    t2=c(i,strsplit(vep[[i]][j],'\\|')[[1]])
    tpos=strsplit(t2[9],'/')[[1]][1]
    taa=strsplit(t2[10],'/')[[1]][2]; if (is.na(taa)) { taa=''; tpos='' }
    t2=c(t2,paste(tpos, taa,sep=''))
    temp=rbind(temp,t2)
  }
  salf$gene[i]=paste(unique(temp[,26]),collapse=',')
  aa=unique(temp[temp[,4] %in% alasccaTX,ncol(temp)]); aa=aa[!aa %in% c('NA','-')]
  salf$aa[i]=paste(aa[aa!=''],collapse = ',')
  table=rbind(table,temp)
}, silent=T)

keep <- salf$t>0.05 & salf$n<0.05 & salf$t.altreads>3
salf <- salf[keep,]
salf$pch=rep(0,nrow(salf))
salf$pch[salf$type=='snv']=1
salf$pch[salf$type=='del']=6
salf$pch[salf$type=='ins']=2


## Process CNA results to file:
# assemble CNA results first
# Starting with PTEN binned values
marg <- 1e6 # margin btw PTEN and a control region
size <- 7e6 # size of control regions
# select PTEN by coordinate instead of by annotation
#bix <- bins$gene=='PTEN' # PTEN bins
bix <- bins$chromosome == "10" & bins$end > 89622870 & bins$start < 89731687
q10 <- bins$chromosome=='10' & bins$start >= chrsizes$mid[which(chrsizes$chr=='10')] # all 10q bins
# 10q control left of PTEN:
left <- q10 & !bix &
  bins$start >= min(bins$start[bix])-marg-size & # start of left control region
  bins$end <= min(bins$start[bix])-marg # end of left control region
# 10q control right of PTEN:
right <- q10 & !bix &
  bins$start >= max(bins$end[bix])+marg & # start of right control region
  bins$end <= max(bins$end[bix])+marg+size # end of right control region
# Also get the CBS segment(s) for PTEN:
#six <- grep('PTEN',segments$gene)
six <- which(segments$chromosome == "10" & segments$end > 89622870 & segments$start < 89731687)
sbix <- NULL; for (i in 1:length(six)) sbix[[i]] <- which(q10 & bins$start >= segments$start[six[i]] & bins$end <= segments$end[six[i]])

p.left <- wilcox.test(x = bins$log2[bix],y = bins$log2[left],alternative = 'less')$p.value # p value for del relative to left side control
p.right <- wilcox.test(x = bins$log2[bix],y = bins$log2[right],alternative = 'less')$p.value # right control
p.seg.left <- NULL; for (i in 1:length(sbix)) p.seg.left <- c(p.seg.left,wilcox.test(x = bins$log2[sbix[[i]]],y = bins$log2[left],alternative = 'less')$p.value) # PTEN segments to left ctrl
p.seg.right <- NULL; for (i in 1:length(sbix)) p.seg.right <- c(p.seg.right,wilcox.test(x = bins$log2[sbix[[i]]],y = bins$log2[right],alternative = 'less')$p.value) # PTEN segments to left ctrl
m.pten <- median(bins$log2[bix],na.rm = T)
m.left <- median(bins$log2[bix],na.rm = T)-median(bins$log2[left],na.rm = T)
m.right <- median(bins$log2[bix],na.rm = T)-median(bins$log2[right],na.rm = T)
m.seg.left <- NULL; for (i in 1:length(sbix)) m.seg.left <- c(m.seg.left,median(bins$log2[sbix[[i]]],na.rm = T)-median(bins$log2[left],na.rm = T))
m.seg.right <- NULL; for (i in 1:length(sbix)) m.seg.right <- c(m.seg.right,median(bins$log2[sbix[[i]]],na.rm = T)-median(bins$log2[right],na.rm = T))
l.seg <- segments$end[six]-segments$start[six]

# PTEN homozygous deletion reported if:
#   Whole PTEN or a small (<5M) PTEN-containing CBS segment is sign. below both left & right controls, with absolute diff at least 0.2 (for 13% drop in DNA, 25% typical purity and delta50%)
#   Thresholds: 0.001 -> base false positive rate of <1/1000 per control (lower due to two tests). With 5% TP, FDR<<2% as amplitude cutoff is also applied.
significance.threshold <- 1e-3
difference.threshold <- -0.2
all.pten.hom.loss <-
  p.left < significance.threshold &
  p.right < significance.threshold &
  m.left < difference.threshold &
  m.right < difference.threshold
seg.pten.hom.loss <- FALSE; if (l.seg < 5e6) for (i in 1:length(m.seg)) {
  seg.pten.hom.loss <-
    seg.pten.hom.loss |
    p.seg.left[i] < significance.threshold &
    p.seg.right[i] < significance.threshold &
    m.left[i] < difference.threshold &
    m.right[i] < difference.threshold
}
pten.hom.loss <- all.pten.hom.loss | seg.pten.hom.loss

# PTEN hemizygous loss and LOH:
six <- grep('PTEN',segments$gene); #seg(s) with PTEN
snps.near <- alf[alf$chromosome==10 & alf$start>=80e6 & alf$end<=100e6,]
snps.ptenseg <- alf[alf$chromosome==10 & alf$start>=min(segments$start[six]) & alf$end<=min(segments$end[six]),]
## call LOH if PTEN segment or vicinity AI is above 0.5 (corr to 50% cells) and with logr below .1
pten.loh <- FALSE
pten.loss <- FALSE

# If there's no SNP data, cowardly refuse to call pten loh or loss
if (any(!is.na(alf$chromosome))){
  pten.loh <-
    m.pten < 0.1 & (
	median(snps.near$ai,na.rm = T) > 0.5 |
	median(snps.ptenseg$ai,na.rm = T) > 0.5
    )
  ## call "LOSS" if PTEN logratio is below -0.4 (some 25% drop, corr to 50% of cells) and AI>.33
  pten.loss <-
    m.pten < -0.4 & (
      median(snps.near$ai,na.rm = T) > 0.33 |
      median(snps.ptenseg$ai,na.rm = T) > 0.33
    )
}

# write to JSON
call='NOCALL'
if (seg.pten.hom.loss|pten.hom.loss) call='HOMLOSS'
if (pten.loh|pten.loss) call='HETLOSS_or_LOH'
CNAlist=list(name='PTEN',call=call)
exportJson <- toJSON(CNAlist)
write(exportJson, opts$json)



#Screen configuration overview
#------------------------------------------------------------
#|      5                       |                           |
#-------------------------------|             4             |
#|      6                       |                           |
#------------------------------------------------------------
#|                               7                          |
#------------------------------------------------------------
#|                               8                          |
#------------------------------------------------------------
#|                               9                          |
#------------------------------------------------------------

png(filename = opts$png,width=11.7,height=8.3,units="in",res=300)
split.screen(figs=c(2,1))
split.screen(as.matrix(data.frame(
  left=c(0,0.25),
  right=c(0.2,.9635),
  bottom=c(0,0.097),
  top=c(1,0.903))),1)
split.screen(figs=c(2,1),3)
split.screen(as.matrix(data.frame(left=c(rep(0.015,3)),
				  right=c(rep(1,3)),
				  bottom=c(0.10,0.5,0.55),
				  top=c(0.5,0.55,1))),2)
split.screen(figs=c(4,6),4)



# ____       _   _   _
#/ ___|  ___| |_| |_(_)_ __   __ _ ___
#\___ \ / _ \ __| __| | '_ \ / _` / __|
# ___) |  __/ |_| |_| | | | | (_| \__ \
#|____/ \___|\__|\__|_|_| |_|\__, |___/
#                            |___/
#Settings

screen(1)
# Sample name and QC
plot(1,type='n',axes=F,xlab='',ylab='')
name=rev(strsplit(opts$cnr,'/')[[1]])[1]
name=strsplit(name,'-panel')[[1]][1]
mtext(name,3,padj=-4.5)
ix=bins$gene=='Background'
mapd1=round(median(abs(diff(bins$log2[!ix]))),2)
mapd2=round(median(abs(diff(bins$log2[ix]))),2)
readRatio=mapd2^2/mapd1^2; seqRatio=readRatio/4; ontar=round(100*seqRatio/(seqRatio+1))
simplePur=round(100*mean(salf$t))
tar=round(median(bins$log2[bins$gene!='Background']),2)
atar=round(median(bins$log2[bins$gene=='Background']),2)
mtext(paste(format(Sys.time(), "%F %H:%M:%S"),'',
	    '  SNPcov:',median(alf$td),
	    '  MAPD:',paste(mapd1,mapd2,sep='/'),
	    #'  medians:',paste(tar,atar,sep='/'),
	    '  Ontarget:',paste(ontar,'%',sep=''),
	    '  MeanSomatcAlleleFr:',paste(simplePur,'%',sep='')),
      3,padj=-4.5,adj=1,cex=0.75)

cex.axis <- .6
cex.mtext <- 1.5
cex.main <- 2
cex.text <- .7
axis1padj <- -2.35
axis2hadj <- 0.3
text1padj <- 1.7#1.7
text2padj<- -3
padj <- -3
lwd=1
color <- '#00000010'




#                                ____
# ___  ___ _ __ ___  ___ _ __   | ___|
#/ __|/ __| '__/ _ \/ _ \ '_ \  |___ \
#\__ \ (__| | |  __/  __/ | | |  ___) |
#|___/\___|_|  \___|\___|_| |_| |____/
#
#screen 5

screen(5)
#plot IGF2 only.
xlim=c(0e6,4e6)
ylim=c(-2,2)

par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
plot(1,type='n',xlim=xlim,ylim=ylim,axes=F,main='IGF2',cex.main=0.5,xaxt='n',yaxt='n',ylab=NA,xlab=NA)
rect(xleft = igf2$start,xright = igf2$end,ybottom = ylim[1],ytop = ylim[2],col = '#0000C030',border='#0000C030')
segments(
  y0=log2(c(.5,1,1.5,2)),
  x0=xlim[1],x1=xlim[2],
  col='grey',
  lwd=1,lty=1)
ix=bins$chromosome=='11' & bins$start>xlim[1] & bins$end<xlim[2]
points(x = (bins$start[ix]+bins$end[ix])/2,y = bins$log2[ix],pch=bins$pch[ix],col=bins$col[ix],cex=0.6)

col='#00C000CC'; ix=segments$chromosome=='11'
segments(x0=segments$start[ix],x1=segments$end[ix],
	 y0=segments$log2[ix],y1=segments$log2[ix],
	 col=col,
	 lwd=2)
axis(1,lwd=lwd,lend=1,cex.axis=cex.axis,padj=axis1padj,tck=-0.03)
axis(2,seq(-1,1,1),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
mtext('Log ratio',2,cex=cex.text,las=0,padj=text2padj)

# somatic muts
ix <- salf$chromosome==11
if (sum(ix)>0) points(x = salf$start[ix], y = rep(-1.7,sum(ix)),pch=salf$pch[ix], cex=0.6, col='#C00000CC')

# #plot snp cov vs alf (removed to fit IGF2)
# xlim=c(50,5000)
# ylim=(0:1)
#
# par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
# plot(1,type='n',xlim=xlim,ylim=ylim,axes=F,main='Heterozygous SNPs',cex.main=0.5,log='x')
# segments(x0=median(alf$td,na.rm = T),y0=0.1,x1=median(alf$td,na.rm = T),y1=.9,lwd=lwd,ylab='',xlab='',lty=3)
# points(alf$td,alf$t,cex=0.3,col='#00000080',xlim=xlim,ylim=ylim,pch=16,lwd=lwd,ylab='',xlab='')
#
# axis(1,c(100,500,2000),lwd=lwd,lend=1,cex.axis=cex.axis,padj=axis1padj,tck=-0.03)
# axis(2,seq(0,1,.25),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
# #mtext('Coverage',1,cex=cex.text,padj=text1padj)
# mtext('SNP allele ratio',2,cex=cex.text,las=0,padj=text2padj)
# text(500,0.05,labels=paste('Median coverage:',round(median(alf$td,na.rm = T)),sep=''),cex=0.5)

#                                 __
# ___  ___ _ __ ___  ___ _ __    / /_
#/ __|/ __| '__/ _ \/ _ \ '_ \  | '_ \
#\__ \ (__| | |  __/  __/ | | | | (_) |
#|___/\___|_|  \___|\___|_| |_|  \___/
#
#screen 6
screen(6)

#plot PTEN only.
xlim=c(89.6e6,89.75e6)
ylim=c(-2.5,1)

main='no loss'
if (CNAlist$call=='HETLOSS_or_LOH') main='LOH/loss'
if (CNAlist$call=='HOMLOSS') main='homozygous loss'
main=paste('PTEN:',main)

par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
plot(1,type='n',xlim=xlim,ylim=ylim,axes=F,main=main,cex.main=0.5,xaxt='n',yaxt='n',ylab=NA,xlab=NA)
rect(xleft = pten$start,xright = pten$end,ybottom = ylim[1],ytop = ylim[2],col = '#0000C030',border='#0000C030')
segments(
  y0=log2(c(.5,1,1.5,2)),
  x0=xlim[1],x1=xlim[2],
  col='grey',
  lwd=1,lty=1)
ix=bins$chromosome==10 & bins$start>xlim[1] & bins$end<xlim[2]
points(x = (bins$start[ix]+bins$end[ix])/2,y = bins$log2[ix],pch=bins$pch[ix],col='black',cex=0.6)
points(x = (bins$start[ix]+bins$end[ix])/2,y = bins$log2[ix],type='l',col='black',lwd=1)

col='#00C000CC'; ix=segments$chromosome=='10'
segments(x0=segments$start[ix],x1=segments$end[ix],
	 y0=segments$log2[ix],y1=segments$log2[ix],
	 col=col,
	 lwd=2)
axis(1,lwd=lwd,lend=1,cex.axis=cex.axis,padj=axis1padj,tck=-0.03)
axis(2,seq(-1,1,1),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
mtext('Log ratio',2,cex=cex.text,las=0,padj=text2padj)

# somatic muts
ix <- salf$chromosome==10
if (sum(ix)>0) {
  points(x = salf$start[ix], y = salf$t[ix],pch=salf$pch[ix], cex=0.6, col='#C00000CC')
  text(x = salf$start[ix],y=-2,labels = salf$aa[ix],cex=0.5,srt=90,col='#C00000CC')
}

#                                _  _
# ___  ___ _ __ ___  ___ _ __   | || |
#/ __|/ __| '__/ _ \/ _ \ '_ \  | || |_
#\__ \ (__| | |  __/  __/ | | | |__   _|
#|___/\___|_|  \___|\___|_| |_|    |_|
#
#screen 4
#allelefreq vs logR

screen(4)
par(mar=c(0,0,0,0))
plot(1,type='n',axes=F,xlab='',ylab='')
mtext('Log ratio',1,padj=text1padj,cex=cex.text)
mtext('Allelic imablance',2,padj=-3,cex=cex.text)

for (c in 1:24)
{
  cex=0.4
  screen(c+9) # for chrY, screen is 33.
  par(mar = c(0, 0, 0, 0))
  par(oma = c(0,0,0,0))
  par(mgp =c(1,0.5,0))
  xlim=c(-1.1,1.1)
  ylim=c(0,1)
  ix <- alf$chromosome %in% c(chrsizes[c]$chr,'X','Y')
  ixCurChr <-  alf$chromosome %in% chrsizes[c]$chr
  plot(alf$log2[!ix],alf$ai[!ix],xlim=xlim,ylim=ylim,lwd=lwd,axes=F,ylab='',xlab='',pch=16,
       col='#B0B0B030',cex=cex)
  points(alf$log2[ixCurChr],alf$ai[ixCurChr],col='#00800070',pch=16,cex=cex)

  d=density(bins$smoothed[bins$chromosome %in% as.numeric(1:22)])
  if (c==24) points(d$x,d$y/max(d$y),type='l')

  whole=log2(c(0.5,1,1.5,2))

  segments(
    x0=whole,x1=whole,
    y0=-0.05,y1=1.035,
    col='#D3D3D360',
    lwd=1)

  segments(
    x0=c(-2),x1=c(2),
    y0=c(1/3,1/2),y1=c(1/3,1/2),
    col='#D3D3D360',
    lwd=1)

  if(c == 23)
  {
    mtext("X", side = 3, line = -1, adj = 0.92, cex = 1)
  } else if(c == 24)
  {
    mtext("", side = 3, line = -1, adj = 0.92, cex = 1)
  } else
  {
    mtext(c, side = 3, line = -1, adj = 0.92, cex = 1)
  }

  if(c %in% 19:24) axis(side=1,cex.axis=0.5,at=seq(-1,1,.5),tck=-0.05,padj=-1.6,#las=3,
			#labels=c('-50%','0','+50%','+100%'),
			col='white',col.ticks='black',lend=1)
  if(c %in% c(1,7,13,19)) axis(side=2,cex.axis=0.6,tck=-0.05,
			       at=c(1/3,1/2),
			       labels=c('1:2','1:3'),
			       las=1,col='white',col.ticks='black',lend=1)
  if(c %in% c(6,12,18,24)) axis(side=4,cex.axis=0.6,tck=-0.05,
				at=c(1/3,1/2),
				labels=c('1:2','1:3'),
				las=1,col='white',col.ticks='black',lend=1)

  box(lwd=1)


}



#                                _ _
# ___  ___ _ __ ___  ___ _ __   / / |
#/ __|/ __| '__/ _ \/ _ \ '_ \  | | |
#\__ \ (__| | |  __/  __/ | | | | | |
#|___/\___|_|  \___|\___|_| |_| |_|_|
#
#screen 11  <- modified, now 9
#LogR plot
screen(9)

#Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
par(mar = c(0, 0, 0, 0))
par(oma = c(0,0,0,0))
par(mgp =c(1,0.5,0))
par(lend=1)

ymin = -2
ymax = 2
#Plot signal over complete genome
plot(NA,NA,#mpos[ix],mval[ix],
     pch=16,
     cex=0.3,
     main='',
     xlab = "",
     ylab = "",
     col = '#00000003',
     xaxt="n",
     axes=F,
     ylim = c(ymin,ymax),
     xlim = c(0,max(bins$cumend))
)


#Add genes as lines
segments(
  x0=(genes$cumstart+genes$cumend)/2,
  y0=-100,y1=100,
  col='#0000C020',
  lwd=2)


#add smoothed signal (OFF)
breaks=c(0,which(abs(diff(bins$start))>1e6),nrow(bins))
# for (i in 2:length(breaks)) {
#   ix=(breaks[i-1]+1):breaks[i]
#   points((bins$cumstart[ix]+bins$cumend[ix])/2,bins$smoothed[ix],type='l',col='#00C00060',lwd=3)
# }
# real signal (omitted above) on top of it
ix=bins$gene=='Background'
points((bins$cumstart[ix]+bins$cumend[ix])/2,bins$log2[ix], #offtargets
       pch=16,
       cex=0.6,
       #type='l',
       col = '#80808080'
)
points((bins$cumstart[!ix]+bins$cumend[!ix])/2,bins$log2[!ix], #ontargets
       pch=1,
       cex=0.7,
       #type='l',
       col = '#00000080'
)


seqminmax <- seq(ymin,ymax,by=0.5)
#Add axis to the left & right of signal
axis(side=2,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=0,las=1)
axis(side=4,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=max(bins$cumend),las=1)
mtext("Log ratio",side=2,line=0,cex=0.7,padj=1.15)

whole=(c(.5,1,1.5,2))
#Add grey segments
segments(
  y0=log2(whole),y1=log2(whole),
  x0=0,x1=max(bins$cumend),
  # col='#D3D3D360',
  col='#00000020',
  lwd=1)

#Add segment for each chromosome.
col='#00C000CC'
segments(x0=segments$cumstart,x1=segments$cumend,
	 y0=segments$log2,y1=segments$log2,
	 col=col,
	 lwd=2)

#Add a bar between chromosomes to distinguish them
segments(
  x0=chrsizes$cumend,
  y0=-100,y1=100,
  col='#00000099',
  lwd=1)



#                                _  ___
# ___  ___ _ __ ___  ___ _ __   / |/ _ \
#/ __|/ __| '__/ _ \/ _ \ '_ \  | | | | |
#\__ \ (__| | |  __/  __/ | | | | | |_| |
#|___/\___|_|  \___|\___|_| |_| |_|\___/
#
#screen 10 <- modified, now 8
#Genes
screen(8)

#Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
par(mar = c(0, 0, 0, 0))
par(oma = c(0,0,0,0))
par(mgp =c(1,0.5,0))

#Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,1),xlim=c(0,max(bins$cumend)))

chrsizes$labelpos[chrsizes$chr=='MT']=NA
text(x = chrsizes$labelpos,y = 0.5,labels = chrsizes$chr,cex=.5)

mtext(text='',side=2,las=1,line=-1.2)

# Gene names
text(x = (genes$cumstart+genes$cumend)/2,y = 0.5,labels = genes$label,srt=45,cex=0.4,col='#0000C0CC')



#                                 ___
# ___  ___ _ __ ___  ___ _ __    / _ \
#/ __|/ __| '__/ _ \/ _ \ '_ \  | (_) |
#\__ \ (__| | |  __/  __/ | | |  \__, |
#|___/\___|_|  \___|\___|_| |_|    /_/
#
#screen 9 <- modified, now 7

# BAF plot
screen(7)

#Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
par(mar = c(0, 0, 0, 0))
par(oma = c(0,0,0,0))
par(mgp =c(1,0.5,0))



#Index to avoid spos/sval values that are all the way at 0 or 1.
#ix = !(sval %in% c(0,1))

#Plot BAF over whole genome
plot(alf$cumstart,alf$t,
     pch=16,
     cex=0.6,
     cex.axis=1,
     main='',
     xlab = "",
     ylab = "",
     col = "#00000080",
     xaxt="n",
     axes=F,
     ylim = c(-0.1,1.1),
     xlim = c(0,max(bins$cumend))
)

#Add axis to the left,right and below of AI. The below axis is the chromosome numbers 1-24.
axis(side=2,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=0,las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
#axis(side=1,at=pre,pos=0,labels=c(seq(from=1,to=22),"X",'Y'),cex.axis=0.50,lty=0)#,tck=0,col.ticks='#00000000')
axis(side=4,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=max(bins$cumend),las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
mtext("SNP allele ratio",side=2,line=0,cex=0.7,padj = 1.15)

#Add a bar between chromosomes to distinguish them
segments(
  x0=chrsizes$cumend,
  y0=-100,y1=100,
  col='#00000099',
  lwd=1)

#Add genes as lines
segments(
  x0=(genes$cumstart+genes$cumend)/2,
  y0=-100,y1=100,
  col='#0000C020',
  lwd=2)

segments(
  y0=c(0,1/4,1/3,2/3,3/4,1),y1=c(0,1/4,1/3,2/3,3/4,1),
  x0=0,x1=max(bins$cumend),
  # col='#D3D3D360',
  col='#00000020',
  lwd=1)

## Add somatic mutations
points(salf$cumstart,salf$t,
     pch=salf$pch,
     cex=0.6,
     col='#C00000CC'
)
text(x = salf$cumstart,y=salf$t-0.07,labels = salf$aa,cex=0.6,srt=30,col='#C00000CC')

#Close all the opened split.screens and release the figure

close.screen(all.screens=T)



dev.off()


# } #used for loops
