######  read arguments  #################

args=commandArgs(trailingOnly = TRUE) 
exp_cutoff=as.numeric(args[1])
somatic=args[2]
output=args[3]
path=args[4]
expression_somatic=args[5]

#######  read neoantigen binding strength  ##################

source(paste(path,"/expressed_neoantigen_function.R",sep=""))
unlink(paste(somatic,"_error_no_valid_neoantigen.txt",sep=""))
neoantigen=read_neoantigen(somatic)
if (is.null(neoantigen) || dim(neoantigen)[1]==0) 
{
  system(paste("echo \"1\" > ",somatic,"_error_no_valid_neoantigen.txt",sep=""))
  cat(paste("Warning: No valid neoantigens found!\n"))
  q()
}
if (is.null(neoantigen)) {stop("No putative neoantigens read!\n")}
  
#######  read neocleotide mutation information  ###########

filtered=paste(somatic,"_filtered",sep="")
if (!file.exists(filtered)) {stop(paste(filtered,"doesn't exist!\n"))}
mutations=read.table(filtered,stringsAsFactors = F,
  colClasses=c("character","numeric","numeric","character","character",rep("numeric",4)))
rownames(mutations)=paste("mutation_",1:dim(mutations)[1],sep="")
neoantigen$chr=paste("chr",mutations[neoantigen$mutation,"V1"],sep="")
neoantigen$pos=mutations[neoantigen$mutation,"V2"]
neoantigen$ref=mutations[neoantigen$mutation,"V4"]
neoantigen$alt=mutations[neoantigen$mutation,"V5"]

##########  read expression data (RPKM)  ################

hla_genes=hla()

if (exp_cutoff>0) # if tumor RNA-Seq data are provided
{
  # read transcript and exon expression
  transcript_exp=read.table(paste(output,"/transcript.featureCounts",sep=""),
                            stringsAsFactors = F,header=T)[,c(1,6,7)]
  colnames(transcript_exp)[3]="expression"
  unlink(paste(output,"/transcript.featureCounts",sep=""))
  transcript_exp=rpkm(transcript_exp)
  
  exon_exp=read.table(paste(output,"/exon.featureCounts",sep=""),
                      stringsAsFactors = F,header=T)[,c(1,2,3,4,6,7)]
  colnames(exon_exp)[6]="expression"
  unlink(paste(output,"/exon.featureCounts",sep=""))
  exon_exp=rpkm(exon_exp)

  # integrate HLA expression
  neoantigen$HLA_exp=NA
  for (i in 1:dim(neoantigen)[1])
  {
    hla_gene=neoantigen$allele[i]
    hla_gene=sub("-","_",sub("\\*.*","",sub(".*\\/","HLA-",hla_gene)))
    keep=transcript_exp$Geneid %in% hla_genes[[hla_gene]]
    neoantigen[i,"HLA_exp"]=
      sum(transcript_exp$expression[keep]*transcript_exp$Length[keep])/
      sum(transcript_exp$Length[keep])*sum(keep)
  }
  
  # integrate transcript and exon level expression 
  transcript=sapply(strsplit(neoantigen$transcript,":"),function(x) x[2])
  neoantigen$transcript_exp=transcript_exp$expression[match(transcript,transcript_exp$Geneid)]
  
  neoantigen$exon_exp=NA
  for (i in 1:dim(neoantigen)[1])
  {
    keep=exon_exp$Geneid==transcript[i] & exon_exp$Start-20<neoantigen$pos[i] &
      exon_exp$End+20>neoantigen$pos[i]
    if (sum(keep)==0) {next} # transcript not found
    neoantigen$exon_exp[i]=exon_exp$expression[which(keep)[1]]
  }
  
  neoantigen=neoantigen[(!is.na(neoantigen$transcript_exp)) & (!is.na(neoantigen$exon_exp)),]
  neoantigen=neoantigen[neoantigen$transcript_exp>exp_cutoff & neoantigen$exon_exp>exp_cutoff,]
}else
{
  neoantigen$transcript_exp=neoantigen$exon_exp=neoantigen$HLA_exp=NA
}

if (dim(neoantigen)[1]==0) 
{
  system(paste("echo \"2\" > ",somatic,"_error_no_valid_neoantigen.txt",sep=""))
  cat(paste("Warning: No valid neoantigens found!\n"))
  q()
}

####### integrate RNA-seq somatic mutation results  #################

neoantigen$rna_normal_var=neoantigen$rna_normal_ref=
  neoantigen$rna_tumor_var=neoantigen$rna_tumor_ref=NA
if (expression_somatic!="NA")
{
  rna_somatic=read.table(expression_somatic,stringsAsFactors = F,sep="\t",header=T)
  for (i in 1:dim(neoantigen)[1])
  {
    j=which(rna_somatic$Chr==neoantigen$chr[i] & 
      abs(rna_somatic$Start-neoantigen$pos[i])<5)[1]
    neoantigen[i,c("rna_tumor_ref","rna_tumor_var","rna_normal_ref",
      "rna_normal_var")]=rna_somatic[j,c("Tumor_ref","Tumor_alt","Normal_ref",
      "Normal_alt")]
  }
}

############  write results  ##################

for (var in c("tumor_ref","tumor_var","normal_ref","normal_var",
  "rna_tumor_ref","rna_tumor_var","rna_normal_ref","rna_normal_var","transcript_exp",
  "exon_exp","HLA_exp"))
  {neoantigen[,var]=round(neoantigen[,var],1)}

write.table(neoantigen[,c("chr","pos","ref","alt","transcript","tumor_ref","tumor_var",
  "normal_ref","normal_var","rna_tumor_ref","rna_tumor_var",
  "transcript_exp","exon_exp","peptide","wildtype_match",
  "allele","HLA_exp","percentile_rank")],file=paste(somatic,"_final.txt",sep=""),
  sep="\t",quote=F,row.names = F)

if ("transcript_exp" %in% ls()) {save(transcript_exp,file=paste(somatic,"_rpkm.RData",sep=""))}
