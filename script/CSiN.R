##########  setting up  ################

args = commandArgs(trailingOnly=TRUE) # a single argument of working folder
setwd(args[1])
neoantigen_object=list(error="No error")

save_and_exit<-function(neoantigen_object)
{
  save(neoantigen_object,file="neoantigen_object.RData")
  q()
}

# possible errors
if (any(grepl("error_no_hla_allele",list.files("."))))
{
  neoantigen_object$error="Alleles cannot be typed"
  save_and_exit(neoantigen_object)
}else if (any(grepl("error",list.files("."))))
{
  neoantigen_object$error="No valid neoantigen"
  save_and_exit(neoantigen_object)
}

# result files
neoantigen_file="./neoantigen_final.txt"
mutation_file="./somatic_mutation_tumor.txt"
clones_file="./clones.txt"
typing_file="./typing.txt"
typing_invalid_file="./typing.txt.invalid"

if (!file.exists(neoantigen_file)) 
{
  neoantigen_object$error="Job not completed"
  save_and_exit(neoantigen_object)
}

###########  read neoantigen file  #######################

neoantigen=read.table(neoantigen_file,stringsAsFactors = F,header=T,
  colClasses=c("ref"="character","alt"="character"))
neoantigen$mutation=paste(neoantigen$chr,neoantigen$pos,neoantigen$ref,neoantigen$alt)
neoantigen$gene=sub("\\:.*","",neoantigen$transcript)
neoantigen=neoantigen[order(neoantigen$mutation,neoantigen$gene,neoantigen$peptide,neoantigen$allele,
  -neoantigen$exon_exp),]
neoantigen=neoantigen[neoantigen$wildtype_match=="none_match",]
neoantigen=neoantigen[!duplicated(cbind(neoantigen$mutation,neoantigen$gene,neoantigen$peptide,
  neoantigen$allele)),]
neoantigen$vaf=neoantigen$tumor_var/(neoantigen$tumor_var+neoantigen$tumor_ref)
neoantigen=neoantigen[!sub("\\*.*","",neoantigen$allele) %in% c("HLA-DQA1","HLA-DRB3","HLA-DRB4","HLA-DRB5"),]
neoantigen_object$neoantigen=neoantigen

###########  read clones  #####################

if (file.exists(clones_file))
{
  clones=read.table(clones_file,sep="\t",header=T,stringsAsFactors = F)[,
    c("cloneId","cloneCount","cloneFraction","clonalSequence","allVHitsWithScore",
    "allDHitsWithScore","allJHitsWithScore","allCHitsWithScore","aaSeqCDR3")]
  if (dim(clones)[1]>500) {clones=clones[1:500,]}
  tr_a=clones[grepl("TRA",clones$allVHitsWithScore),]
  tr_b=clones[grepl("TRB",clones$allVHitsWithScore),]
  ig_heavy=clones[grepl("IGH",clones$allVHitsWithScore),]
  ig_light=clones[grepl("IGK",clones$allVHitsWithScore) | grepl("IGL",clones$allVHitsWithScore),]
    
  tr_a$cloneFraction=tr_a$cloneFraction/sum(tr_a$cloneCount)
  tr_b$cloneFraction=tr_b$cloneFraction/sum(tr_b$cloneCount)
  ig_heavy$cloneFraction=ig_heavy$cloneFraction/sum(ig_heavy$cloneCount)
  ig_light$cloneFraction=ig_light$cloneFraction/sum(ig_light$cloneCount)
  neoantigen_object$clones=list(tr_a=tr_a,tr_b=tr_b,ig_heavy=ig_heavy,ig_light=ig_light)
}else
{
  neoantigen_object$clones=NA
}
  
#######  mutation  ########################

mutations=read.table(mutation_file,sep="\t",header=T,stringsAsFactors = F,
  colClasses=c("Ref"="character","Alt"="character"))
mutations=mutations[mutations$Normal_alt/(mutations$Normal_alt+mutations$Normal_ref)<=0.02,]
mutations=mutations[mutations$Tumor_alt/(mutations$Tumor_alt+mutations$Tumor_ref)>=0.05,]
rownames(mutations)=paste(mutations$Chr,mutations$Start,mutations$Ref,mutations$Alt)
mutations$vaf=mutations$Tumor_alt/(mutations$Tumor_alt+mutations$Tumor_ref)
neoantigen_object$mutations=mutations

######  typed alleles  ###################

neoantigen_object$typed=read.table(typing_file,sep="\t",header=F,stringsAsFactors = F)
colnames(neoantigen_object$typed)=c("Class","Allele1","Allele2","Distance")
neoantigen_object$typed_invalid=readLines(typing_invalid_file)
neoantigen_object$typed_invalid=sub("\\.p'","",sub(".*HLA","HLA",neoantigen_object$typed_invalid))
valid_alleles=c("A","B","C","DRB1","DQB1") 
neoantigen_object$valid_type_n=sum(valid_alleles %in% neoantigen_object$typed$Class)*2-
  sum(sub("HLA-","",sub("\\*.*","",neoantigen_object$typed_invalid)) %in% valid_alleles)
if (neoantigen_object$valid_type_n<8) 
{
  neoantigen_object$error="Too few valid alleles with binding affinity data"
  save_and_exit(neoantigen_object)
}

########  plotting function  #####################

#gradient color bar
color.bar <- function(lut,nticks=2, ticks=seq(0,50,10), title='') {
  scale = (length(lut)-1)/50
  #  dev.new(width=1.75, height=5)
  plot(c(-1,10), c(0,50), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, at=seq(0,20,10),labels=seq(0,78,39),las=1,ticks,col.ticks='white',col='white',cex.axis=0.6)
  for (i in 1:(length(lut)-1)) {
    y = ((i-1)/scale)*0.4
    rect(-1.5,y,-1,y+1/scale, col=lut[i], border=NA)
  }
}

pie_<- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
                 init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                 col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = NA, col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    #lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
    #text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, adj = ifelse(P$x < 0, 1, 0), ...)
  }
  title(main = main, ...)
  invisible(NULL)
}

###########  CSiN  ##################

pdf(file=paste(args[1],"/CSiN.pdf",sep=""),width=4,height=4)
colfunc <- colorRampPalette(c("white","black"))
cutoffs=c(3,4,5,6,10,14,16)/8

neoantigen_object$CSiN=10*mean(sapply(rev(cutoffs),function(cutoff){
  # sub-CSiN score
  neoantigen=neoantigen_object$neoantigen[neoantigen_object$neoantigen$percentile_rank<=cutoff,]
  neoantigen_count=table(neoantigen$mutation)
  neoantigen_count=neoantigen_count/mean(neoantigen_count)
  mutations=neoantigen_object$mutations[unique(neoantigen$mutation),]
  mutations$vaf=mutations$vaf/mean(mutations$vaf)
  mutations=mutations[rank(-mutations$vaf)<500,]
  neoantigen_count=neoantigen_count[names(neoantigen_count) %in% rownames(mutations)]
  score=log(mean(neoantigen_count*mutations[names(neoantigen_count),"vaf"]))
  
  # sub-plot
  if (cutoff==2 || score!=0)
  {
    index=round(mutations[names(neoantigen_count),"vaf"]*10+2)
    index[index>30]=30
    pie_(x=neoantigen_count,col=colfunc(30)[index],
      radius=sum(cutoff>=cutoffs)/length(cutoffs),
      labels=ifelse(cutoff==2,'Cauchy-Schwarz index of Neoantigens',""))
  }

  # return
  par(new=TRUE)  
  score
}),na.rm=T)

color.bar(colorRampPalette(c("white","black"))(30))
dev.off()

cat("CSiN=",neoantigen_object$CSiN,"\n")
save_and_exit(neoantigen_object)
