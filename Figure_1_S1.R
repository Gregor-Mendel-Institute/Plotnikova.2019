#Figure_1_S1.R

#Set filepaths according to your system
setwd('/Users/michael.nodine/Desktop/R_Desktop/embryonic_miRNAs/Plotnikova.2019/Figure_1_S1/')           #Home directory
fileRoot = '/Users/michael.nodine/Desktop/R_Desktop/embryonic_miRNAs/Plotnikova.2019/Figure_1_S1/'
quantRoot = '/Users/michael.nodine/Desktop/R_Desktop/embryonic_miRNAs/Plotnikova.2019/quant/'            #Where *tsv miRNA alignment files are located
annoRoot = '/Users/michael.nodine/Desktop/R_Desktop/embryonic_miRNAs/Plotnikova.2019/Ath_annotations/'   #Where annotation files are located (quantTable_140126.txt, )
resourceRoot = '/Users/michael.nodine/Desktop/R_Desktop/embryonic_miRNAs/Plotnikova.2019/resources/'

get.quant.by.sample <- function(sample,align) {
  sample.df = read.delim(paste(quantRoot,sample,'/',align,'/',sample,'.quant.on.miRNA_mature.tsv',sep=''), 
                         row.names=1, header=TRUE, stringsAsFactors=FALSE,strip.white=TRUE)
  return(sample.df)
}

get.GMR <- function(sample,align) {

  gmr.df = read.delim(paste(resourceRoot,'/emb.fb.lf.d1.GMR.tsv',sep=''), 
                      row.names=1, header=TRUE, stringsAsFactors=FALSE,strip.white=TRUE)
  GMR=gmr.df[sample,align]
  return(GMR)
}

get.dose.response <- function(sample,dil.factor,align) {

  #get reference table
  spike.ref=read.delim(paste(resourceRoot,'/quantTable_140126.txt',sep=''),
                       row.names=1, header=FALSE, skip=1, stringsAsFactors=FALSE,strip.white=TRUE)
  
  names(spike.ref)=c("amol","mlc.num")
  
  #get all detected spikes in sample
  ##get number of genome-matching reads for normalization
  GMR=get.GMR(sample,align)
  MGMR = GMR/1000000
  
  spike.root=paste(clusterRoot,sample,'/spikes/',sep="")
  
  spike.name.v = c("433","403","361","871","974","71","823","87")
  
  rownames(spike.ref)=spike.name.v
  
  spike.files=list.files(spike.root)
  
  rpm.v=c()
  
  for (spike in rownames(spike.ref)) {
    spike.file=paste("srna_spike.",spike,".match.tsv",sep="")
    
    #if (spike.file %in% spike.files == TRUE) {
    lineNum=readLines(paste(spike.root,spike.file,sep=''))
    
    if (length(lineNum) > 0) {
      spike.obs=read.delim(paste(spike.root,spike.file,sep=''),
                           row.names=1, header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE)
      #now normalize for number of genome-matching reads
      spike.obs.norm=spike.obs/MGMR
      
      rpm.v=c(rpm.v,sum(spike.obs.norm[,1]))
    }else {
      rpm.v=c(rpm.v,0)
    }
  }

  spike.ref.rpm=data.frame(spike.ref,RPM=rpm.v)
  
  spike.ref.rpm.new=data.frame(row.names=rownames(spike.ref.rpm), 
                               amol=spike.ref.rpm$amol * dil.factor, 
                               mlc.num=spike.ref.rpm$mlc.num * dil.factor * 1000000,
                               RPM=spike.ref.rpm$RPM)
  
  return(spike.ref.rpm.new)
  
}

get.5pMatrix.GMR <- function(sample,align) {


  data = read.delim(paste(quantRoot,sample,'/',align,'/',sample,'.norm.bed',sep=''),sep='\t', row.names=NULL, header=FALSE, stringsAsFactors=FALSE,strip.white=TRUE)
  names(data) = c('chromo','start','end','ID','RPM',"strand","seq","size")
  
  data.new=data.frame(seq=data$seq,size=data$size,fivePbase=substr(data$seq,1,1),RPM=data$RPM)
  
  fivePmatrix = matrix(1:(13 * 4) * 0, nrow=4, ncol=13, byrow=TRUE, dimnames=list(c('A','C','G','T'),18:30))
  
  for (i in 18:30) {
    for (base in c('A','C','T','G')) {
      fivePmatrix[base,as.character(i)] = sum(subset(data.new, size == i & fivePbase == base)$RPM)
    }
  }
  
  #normalize by proportion
  fivePmatrix.norm = fivePmatrix/sum(fivePmatrix) * 100
  
  dimnames(fivePmatrix) = list(c('A','C','G','U'),18:30)
  dimnames(fivePmatrix.norm) = list(c('A','C','G','U'),18:30)
  
  return(list(fivePmatrix,fivePmatrix.norm))
}

get.mir.fams <- function() {
  fam.df = read.delim(paste(annoRoot,'/annoFams_150828',sep=''), skip=1, 
                         row.names=1, header=TRUE, stringsAsFactors=FALSE,strip.white=TRUE)
  return(fam.df)  
}

################################################^^^FUNCTIONS^^^################################################

#Figure 1A and 1B; Supplemental Figure 1A to 1C

matrix.list.star.0=list(get.5pMatrix.GMR("500ng","star_0_mismatches"),get.5pMatrix.GMR("50ng","star_0_mismatches"),
                        get.5pMatrix.GMR("5ng","star_0_mismatches"),get.5pMatrix.GMR("1ng","star_0_mismatches"),
                        get.5pMatrix.GMR("0dot5ng","star_0_mismatches"))

plot.length <- function(matrix.list,align) {
  plotFile=paste(fileRoot,'graphs/dilution.series.length.distributions.',align,'.git.pdf',sep='')
  pdf(plotFile,useDingbats=FALSE, height=2.5, width=11)
  par(mfrow=c(1,5))
  
  ymax=500
  
  plot.length.sub <- function(matrix.df,amt) {
    col_v = c('green','blue','yellow','red')
    
    matrix.df.new=matrix.df/1000
    
    barplot(height=matrix.df.new, col = col_v, main="", xlab='Length (nt)', ylab='Reads per thousand', las=1, ylim=c(0,ymax))
    legend('topright', dimnames(matrix.df.new)[[1]], fill=col_v, box.lty=0)
    text(2.5, ymax, labels=amt, pos=1)
  }
  
  plot.length.sub(matrix.list[[1]][[1]],"500 ng")
  plot.length.sub(matrix.list[[2]][[1]],"50 ng")
  plot.length.sub(matrix.list[[3]][[1]],"5 ng")
  plot.length.sub(matrix.list[[4]][[1]],"1 ng")
  plot.length.sub(matrix.list[[5]][[1]],"0.5 ng")
  
  dev.off()
  
}

plot.length(matrix.list.star.0,"star_0_mismatches")


#Figure 1C and Supplemental Figure 1D to 1H

sample.v=c("bc_32","bc_34","500ng","50ng","5ng","1ng","0dot5ng")

get.master.mir.df <- function(align) {
  #align=""
  sample.mir.list = list()
  master.mir.df = data.frame(row.names=rownames(get.quant.by.sample("bc_32",align)))
  
  for (sample in sample.v) {
    mir.df=subset(get.quant.by.sample(sample,align), select=X20to24.norm)
    names(mir.df)=sample
    sample.mir.list[[sample]]=(mir.df)
    master.mir.df=cbind(master.mir.df,mir.df)
  }
  
  #group into families (more relevant to examine per family)
  mir.fams=get.mir.fams()
  
  #assign to master.mir.df
  master.mir.df.plus.fams=data.frame(master.mir.df,mir.fam=mir.fams[rownames(master.mir.df),'fam'])
  
  uni.fams=as.character(unique(master.mir.df.plus.fams$mir.fam))
  
  master.mir.df.fams=data.frame()
  
  for (fam in uni.fams) {
    #fam="miR160"
    
    ind.fam.df=subset(master.mir.df.plus.fams, mir.fam==fam)
    
    if (length(ind.fam.df[,1] > 1)) {
      ind.fam.row=data.frame(row.names=fam,t(colSums(ind.fam.df[,1:7])))
    }else {
      ind.fam.row=data.frame(row.names=ind.fam.df$mir.fam,ind.fam.df[,1:7])
    }
    
    master.mir.df.fams=rbind(master.mir.df.fams,ind.fam.row)
  }
  
  return(master.mir.df.fams) 
}

master.mir.df.fams.0 = get.master.mir.df("star_0_mismatches")

master.mir.df.fams=data.frame()

plot.repro <- function(outFile,rpm.thresh,master.df.fams) {

  plotFile=paste(fileRoot,'graphs/',outFile,sep='')
  pdf(plotFile,useDingbats=FALSE, height=3, width=16)
  par(mfrow=c(1,6))
  
  plot.dil.repro <- function(rpm.thresh,select,amt,master.df.fams) {

    dil.df.sub=subset(master.df.fams, master.df.fams[,select] >= rpm.thresh | master.df.fams[,"X500ng"] >= rpm.thresh, select=c("X500ng",select))
    
    dil.df.sub.log = log(dil.df.sub + 1, 10)
    
    ymax=max(dil.df.sub.log)
    
    plot(dil.df.sub.log[,2], dil.df.sub.log[,1], ylab="miRNA levels in 500 ng (RPM; log10)", main="miRNA levels",
         xlab=paste("miRNA levels in",amt,"ng (RPM; log10)"), las=1, pch=1, xlim=c(0,ymax), ylim=c(0,ymax))
    cor.value = round(cor.test(dil.df.sub.log[,1],dil.df.sub.log[,2], method="pearson")$estimate,2)
    text(x=0.5, y=ymax, labels=paste("R =",cor.value), pos=1)
    abline(a=0, b=1, lty=2)
  }
  
  plot.dil.repro(rpm.thresh,"bc_32","Biorep_2",master.df.fams)
  plot.dil.repro(rpm.thresh,"bc_34","Biorep_3",master.df.fams)
  plot.dil.repro(rpm.thresh,"X50ng","50",master.df.fams)
  plot.dil.repro(rpm.thresh,"X5ng","5",master.df.fams)
  plot.dil.repro(rpm.thresh,"X1ng","1",master.df.fams)
  plot.dil.repro(rpm.thresh,"X0dot5ng","0.5",master.df.fams)
  
  dev.off()
  
}

plot.repro("dilution.series.repro.fams.0rpm.star_0_mismatches.git.pdf",0,master.mir.df.fams.0)


#Figure 1D and Supplemental Figure 1I-1L

##Use number of molecules added prior to dilution: 5.26 ul * 1/(4 ul of mix set) = dilution factor = 1.315 * dil table = 

star.0.dr.list = list(get.dose.response("500ng",1.315,"star_0_mismatches"),get.dose.response("50ng",1.315,"star_0_mismatches"),
                    get.dose.response("5ng",1.315,"star_0_mismatches"),get.dose.response("1ng",1.315,"star_0_mismatches"),
                    get.dose.response("0dot5ng",1.315,"star_0_mismatches"))

plot.spike.dr <- function(dr.list,align) {
  plotFile=paste(fileRoot,'graphs/dilution.series.spike.ins.per.added.',align,'.pdf',sep='')
  pdf(plotFile,useDingbats=FALSE, height=2.5, width=11)
  par(mfrow=c(1,5))
  
  cor.value.v=c()
  
  col.v=c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef")
  
  mlc.num.v=c(dr.list[[1]]$mlc.num,dr.list[[2]]$mlc.num,dr.list[[3]]$mlc.num,dr.list[[4]]$mlc.num,dr.list[[5]]$mlc.num)
  mlc.num.v.log=log(mlc.num.v + 1,10)
  rpm.v=c(dr.list[[1]]$RPM,dr.list[[2]]$RPM,dr.list[[3]]$RPM,dr.list[[4]]$RPM,dr.list[[5]]$RPM)
  rpm.v.log=log(rpm.v + 1,10)
  
  ymin=min(mlc.num.v.log)
  xmin=0 #min(rpm.v.log)
  
  ymax=max(mlc.num.v.log)
  xmax=max(rpm.v.log)
  
  plot.dil.spikes <- function(df2,select,amt,color,add) {
    rpm.thresh=0
    
    dil.df.sub=subset(df2, RPM >= rpm.thresh, select=c(mlc.num,RPM))
    
    dil.df.sub.log = log(dil.df.sub + 1, 10)
    
    if (add==FALSE) {
      plot(dil.df.sub.log$RPM, dil.df.sub.log$mlc.num, xlab="Reads per million (log10)", main="Small RNA spike-ins",
           ylab="Number of molecules (log10)", las=1, pch=19, xlim=c(xmin,xmax), ylim=c(ymin,ymax), col=color)
    }else {
      points(dil.df.sub.log$RPM, dil.df.sub.log$mlc.num, pch=19, col=color)
    }
    
    #abline(v=1,lty=2)
    
    cor.value = round(cor.test(dil.df.sub.log$RPM,dil.df.sub.log$mlc.num, method="pearson")$estimate,2)
    #cor.value.v=c(cor.value.v,cor.value)
    text(0.75,ymax, labels=paste(amt,"ng\nR =",cor.value), pos=1)
    fit = lm(dil.df.sub.log$mlc.num ~ dil.df.sub.log$RPM)
    abline(fit, col="black", lty=2, lwd=1)
    
    
  }
  
  plot.dil.spikes(dr.list[[1]],"X500","500",col.v[1],FALSE)
  plot.dil.spikes(dr.list[[2]],"X50","50",col.v[2],FALSE)
  plot.dil.spikes(dr.list[[3]],"X5","5",col.v[3],FALSE)
  plot.dil.spikes(dr.list[[4]],"X1","1",col.v[4],FALSE)
  plot.dil.spikes(dr.list[[5]],"X05","0.5",col.v[5],FALSE)
  
  dev.off()
  
}

plot.spike.dr(star.0.dr.list,"star_0_mismatches.git")
