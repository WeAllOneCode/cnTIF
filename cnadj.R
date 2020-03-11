#Read in genome-wide matrix (trans component) and copy number predicted by readDepth https://github.com/chrisamiller/readDepth
library(ggplot2); #for plotting with alpha, not required
library(plyr); #for round_any


libname <- 'HF2354'
bpbin  <- 50000
binsize=as.numeric(bpbin);



#--------------------------------------- functions --------------------------------------

splitMat <- function(GWmat, idxbin){
#given genome-wide matrix, normalized cis and trans separately
  n = nrow(GWmat);
  chrs = as.vector(unique(idxbin$CHR));
  chr_indexes = lapply(chrs, function(x){which(idxbin$CHR==x)});
  names(chr_indexes) = chrs;
  mat.cis = GWmat
  mat.tran = GWmat
  for ( ch in chrs){
      j=chr_indexes[[ch]];#indexes in matrix intrachr
      mat.cis[j,-j] = NA
      mat.tran[j,j] = NA
   }
   return(list(cis=mat.cis, tran=mat.tran));
}

getIndex <- function(idx, ch, pos1, pos2){
#given coordinates, return index row of idx
   ch = as.character(ch)
   idx1 = subset(idx, CHR==ch & Start < pos2 & End > pos1);
   return(as.numeric(rownames(idx1)));
}


transBinningCN <- function(cnvtable, otherBin, defaultCN=NA){
#estimate CN for otherBin
    colnames(cnvtable) = c("chrom", "loc.start", "loc.end", "seg.mean")
    nbin = nrow(otherBin)
    cnbin = numeric(nbin);
    for(i in 1:nbin){
        ch = as.character(otherBin[i,1]);
        s = otherBin[i,2]
        e = otherBin[i,3]
        scn = subset(cnvtable, chrom==ch  & loc.start <= e & loc.end >= s) #overlap to query is fine
        nr = nrow(scn);
        if(nr >0){
            cnbin[i] = round(mean(scn$seg.mean),3);
        }else{
            cnbin[i] = defaultCN;
        }
    }
    return(data.frame(otherBin[,1:3], cn=cnbin));
}

plotVector <- function(gwvec, idxbin, dotpos=NULL, tt='', ylabel='', vcol='darkblue', ymax=max(gwvec, na.rm=T), ymin=min(gwvec, na.rm=T), dotcol='#d400e2'){
#given a genome-side vector, plot on chr axis
#global var: ch.start, ch.end
    chrs = as.vector(unique(idxbin$CHR));
    chr_indexes = lapply(chrs, function(x){which(idxbin$CHR==x)});
    names(chr_indexes) = chrs;
    plot(gwvec, main=tt, ylab=ylabel, type='h', col=vcol, axes=FALSE, xlab='Chromosome', ylim=c(ymin,ymax), yaxs="i", xaxs="i");
    points(dotpos, gwvec[dotpos], col=dotcol, pch=16, cex=1);
    abline(v=c(ch.start, ch.end), lty=2, col='grey')
    axis(2)
    axis(1, at=(ch.start+ ch.end)/2, labels=chrs)
}

# ------------------------------ function end -------------------

idx = read.table('hg19_i50000.ind', header=T, row.names=paste0('INDEX', bpbin));
nbin = nrow(idx);
chrs = as.vector(unique(idx$CHR));
chr_indexes = lapply(chrs, function(x){which(idx$CHR==x)});
names(chr_indexes) = chrs;
regionNames = paste(sub('chr','',idx$CHR),idx$Start/binsize, sep='_');
ch.start = unlist(lapply(chr_indexes, head,1));
ch.end = unlist(lapply(chr_indexes, tail,1));

matfile = paste0(libname,paste0('_i',bpbin,'_rawPET0.mat')); #MODIFY matrix file name
message("Reading genome wide matrix: ", matfile)
gwmat = as.matrix(read.table(matfile)); #reading matrix; takes a long time (~30 minutes) for 50kb bins or 60739x60739
gwmat = gwmat[-1,-1]# the headers are bin order on each chromosome
colnames(gwmat) = rownames(gwmat) = regionNames; #takes some cpu time

n = nrow(gwmat);

ncis = sum(unlist(lapply(chr_indexes, length))^2); #total matrix elements in cis
ntran =  nrow(idx)^2 - ncis
m = splitMat(gwmat, idx)
matnorm.tran = (m$tran/sum(m$tran, na.rm=T))*ntran;

aveNtif = list(); #sums of normalized TIF
for (ch in chrs){
    j=chr_indexes[[ch]];#
    mat2  = matnorm.tran[j,-j]
    aveNtif[[ch]] = rowMeans(mat2, na.rm=T);
}

#CN-adjustment

###########################
#get blacklisted bins
blacklist = read.table('z3HF_ENCODE.blacklist.bed', stringsAsFactors=F);
blacklist = read.table('z3_encode_gap.blacklist.bed', stringsAsFactors=F);
blacklist = subset(blacklist, V1 %in% chrs)
iblack = NULL;
for ( i in 1:nrow(blacklist) ){
  iblack = c(iblack, getIndex(idx, blacklist[i,1], blacklist[i,2], blacklist[i,3]));
}
iblack = unique(sort(iblack));
###########################
ecdnaReg0 = read.table('ecDNA_regions.txt', stringsAsFactors = F);

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

### START libname here ###
ecdnaReg = subset(ecdnaReg0, V4==libname)
iec = NULL;
for ( i in 1:nrow(ecdnaReg) ){
  iec = c(iec, getIndex(idx, ecdnaReg[i,1], ecdnaReg[i,2], ecdnaReg[i,3]));
}
iec = unique(iec);
iexc = unique(sort(c(iblack, iec)));

lnum = substr(libname, 1,2)
bglabel=paste0("bg")
idlabel = rep(bglabel,nbin); idlabel[iec] = "EC";

aveNtif = unlist(aveNtif)

### Getting copy number ###
regionCN = read.table(paste0('cnv.',libname,'.bedgraph'), stringsAsFactors=F);
regionCN = read.table(paste0(libname,'.subtractBlacklist.CN.bdg'), stringsAsFactors=F);
bcn = transBinningCN(regionCN, idx)
bcntif = data.frame(bcn, tif=aveNtif);
bcntif$label = factor(idlabel);
ina = union(which(is.na(bcntif$cn)),iblack); #NA and blacklist
bcntif$cn[ina] = NA;
bcntif$tif[ina] = NA;
bcntif$roundCN = round_any(bcntif$cn,1);
###########################

obox = bcntif[-union(ina,iec),]; obox = subset(obox, cn>=0.5);
cnn = table(obox$roundCN);#monitor to remove groups with too small for boxplot
nr=nrow(obox)


ltif = lm(obox$tif ~ obox$cn); #pol2 line
tifcn.cor = cor.test(obox$cn, obox$tif, na.rm=T);
reglines = data.frame(ic=ltif$coefficients[1], m=ltif$coefficients[2], chip="Pol2")

#Scatter plot with fitted line
col5 = "#db72a8"
col.lr = 'skyblue'
col.ec = 'darkred'

plot(obox$cn, obox$tif, cex=2, pch=20, col=alpha(col5, 0.5), main=paste("Pol2",libname), ylab="average nTIF", xlab="CN");
abline(ltif, lty=2, lwd=2, col=col.lr)
legend("topright", legend=paste("Pearson's r:", round(tifcn.cor$estimate, 3)), lty=2, col=col.lr,lwd=2, bty="n")

#Scatter plot regression & display EC
x = bcntif$cn[iec]; ypol2 = bcntif$tif[iec];
plot(obox$cn, obox$tif, cex=2, pch=20, col=alpha(col5, 0.5), main=paste("Pol2",libname), ylab="average nTIF", xlab="CN", xlim=c(1,max(x)), ylim=c(0,max(ypol2)));
points(x, ypol2, pch=8, col=col.ec)
abline(ltif, lty=2, lwd=2, col=col.lr)
legend("topleft", pch=8, col=col.ec, legend="EC")


#genome wide profile of TIF
res = bcntif$tif - (ltif$coefficients[1] + ltif$coefficients[2]*bcntif$tif); res[res<0] = 0;

plotVector(bcntif$tif, idx, vcol = col5, tt=paste("Pol2",libname), ylabel = "average nTIF", ymax=1.1*max(bcntif$tif, na.rm=T));
points(iec, bcntif$tif[iec], pch=8, col=col.ec)

plotVector(res, idx, vcol = col5, tt=paste("Pol2 ChIA-PET",libname), ylabel = "cnadjTIF", ymax = 1.2*max(res,na.rm=T));
points(iec, res[iec], pch=8, col=col.ec);

#The result of this script is bcntif.
bcntif$cnadjTIF = res; #score of interchromosomal interaction adjusted by copy number variation
head(bcntif)
#cnadjTIF =< 0 might mean no enrichment of interchromosomal interaction for that bin.
write.table(bcntif, file=paste0(libname, '.cntif.txt'), quote=F, col.names=T, row.names=F, sep="\t")




