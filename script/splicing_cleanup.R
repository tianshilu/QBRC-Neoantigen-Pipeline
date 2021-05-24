######  read arguments  #################

args=commandArgs(trailingOnly = TRUE) 
somatic=args[1]
output=args[2]
path=args[3]
debug=as.numeric(args[4])

#######  read neoantigen binding strength  ##################

source(paste(path,"/expressed_neoantigen_function.R",sep=""))
unlink(paste(somatic,"_error_no_valid_neoantigen.txt",sep=""))

if (debug==1) 
{
  system(paste("mkdir ",somatic,sep=""))
  system(paste("cp ",somatic,"*txt ",somatic,sep=""))
}

neoantigen=read_neoantigen(somatic)
if (is.null(neoantigen) || dim(neoantigen)[1]==0) 
{
  system(paste("echo \"1\" > ",somatic,"_error_no_valid_neoantigen.txt",sep=""))
  cat(paste("Warning: No valid neoantigens found!\n"))
  q()
}
if (is.null(neoantigen)) {stop("No putative neoantigens read!\n")}

#######  final cleanup  ###################

# keep columns related to splicing
neoantigen=neoantigen[,c(1:3,9:11)]
colnames(neoantigen)[1:3]=c("protein_id","start","end")
neoantigen$transcript_id=sub("\\.p[0-9]+","",neoantigen$protein_id,perl=T)
neoantigen$gene_id=sub("\\.[0-9]+$","",neoantigen$transcript_id,perl=T)

# read tumor out.gtf
gtf=read.table(paste(output,"/tumor/out.gtf",sep=""),stringsAsFactors=F,sep="\t")
gtf=gtf[gtf$V3=="transcript",]
gtf$transcript_id=sub(";.*","",sub(".*transcript_id\\s","",gtf$V9,perl=T))
gtf$related_transcript=sub(";.*","",sub(".*reference_id\\s","",gtf$V9,perl=T))
gtf$related_transcript[!grepl("reference_id",gtf$V9)]=""
gtf$FPKM=sub(";.*","",sub(".*FPKM\\s","",gtf$V9,perl=T))
gtf$gene_id=sub("\\.[0-9]+$","",gtf$transcript_id,perl=T)
rownames(gtf)=gtf$transcript_id

# FPKM
neoantigen$FPKM=as.numeric(gtf[neoantigen$transcript_id,"FPKM"])

# related reference transcripts
tmp=aggregate(gtf$related_transcript,by=list(gtf$gene_id),
  function(x) paste(x[1:min(5,length(x))],collapse=","))
tmp$x=sub("^,+","",tmp$x,perl=T)
tmp$x=sub(",+$","",tmp$x,perl=T)
tmp$x=gsub(",+",",",tmp$x,perl=T)
rownames(tmp)=tmp$Group.1
neoantigen$related_transcripts=tmp[neoantigen$gene_id,"x"]

# write to file
write.table(neoantigen[,c("gene_id","transcript_id","protein_id","start",
  "end","peptide","allele","percentile_rank","FPKM","related_transcripts")],
  file=paste(output,"/neoantigen_splicing_final.txt",sep=""),col.names=T,
  row.names=F,quote=F,sep="\t")
