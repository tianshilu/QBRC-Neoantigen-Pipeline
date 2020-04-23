#############  read data  ######################

args <- commandArgs(trailingOnly = TRUE)
input=args[1]
tumor_freq=args[2]
normal_freq=args[3]
output=args[4]
rpkm_cutoff=as.numeric(args[5])
code_path=args[6]
max_mutations=as.numeric(args[7])

unlink(output)
if (!file.exists(input)) {stop(paste(input,"doesn't exist!\n"))}
mutations=read.table(input,sep="\t",stringsAsFactors = F,header=T,
                     colClasses=c("Ref"="character","Alt"="character"))

#########  filter SNPs  ######################

# unique mutations in each patient
mutations$Variant=apply(mutations[,c("Chr","Start","Ref","Alt")],1,function(x) paste(x,collapse="_"))
mutations=mutations[!duplicated(mutations$Variant),]

# non-silent mutations
mutations=mutations[mutations$Func.refGene %in% c("exonic","exonic;splicing",
                                                  "splicing;exonic"),]
mutations=mutations[mutations$ExonicFunc.refGene %in% c("nonsynonymous SNV",
  "frameshift substitution","nonframeshift substitution","stoploss"),]

# filter SNPs by allele frequencies
mutations=mutations[mutations$Normal_alt/(mutations$Normal_alt+mutations$Normal_ref)<
  as.numeric(normal_freq),]
mutations=mutations[mutations$Tumor_alt/(mutations$Tumor_alt+mutations$Tumor_ref)>
  as.numeric(tumor_freq),]
if (dim(mutations)[1]==0) 
{ 
  warning("No somatic mutations in tumor sample are left after all filtering(1)!")
  q()
}

# if a gene has too many mutations in one patient, it is likely to be false or a real complicated varaint
# and will cause a problem in protein translation
tmp=table(mutations$Gene)
mutations=mutations[!mutations$Gene %in% names(tmp[tmp>=5]),]

# pre-filter by transcript expression, not accurate
if (rpkm_cutoff>0)
{
  source(paste(code_path,"/expressed_neoantigen_function.R",sep=""))
  transcript_file=sub("neoantigen$","transcript.featureCounts",input,perl=T)
  transcript_exp=read.table(transcript_file,stringsAsFactors = F,header=T)[,c(1,6,7)]
  colnames(transcript_exp)[3]="expression"
  transcript_exp=rpkm(transcript_exp)
  transcript_split=strsplit(mutations$AAChange.refGene,",")
  mutations$keep=T
  for (i in 1:dim(mutations)[1])
  {
    transcripts=
      sapply(1:length(transcript_split[[i]]),function(j) strsplit(transcript_split[[i]][j],":")[[1]][2])
    if (sum(transcript_exp[,1] %in% transcripts)<1 || 
        all(transcript_exp[transcript_exp[,1] %in% transcripts,3]<rpkm_cutoff))
      {mutations$keep[i]=F}
  }
  mutations=mutations[mutations$keep==T,]
}

if (dim(mutations)[1]==0 || dim(mutations)[1]>max_mutations) 
{ 
  warning(paste("No or too many somatic mutations (",dim(mutations)[1],
    ") in tumor sample are left after all filtering!",sep=""))
  q()
}

##########  write data  ################

mutations$Chr=sub("chr","",mutations$Chr) 
write.table(mutations[,c("Chr","Start","End","Ref","Alt","Tumor_ref","Tumor_alt","Normal_ref",
  "Normal_alt")],file=output,sep="\t",quote=F,row.names = F,col.names =F)
