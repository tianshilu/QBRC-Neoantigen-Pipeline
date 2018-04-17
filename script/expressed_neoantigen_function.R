rpkm<-function(x)
{
  if (sum(x$expression)==0) {stop("Error: gene expression count is 0!")}
  x$expression=x$expression/sum(x$expression)*10e6
  x$expression=x$expression/x$Length*1000
  x$Geneid=sub("\\.[0-9]+","",x$Geneid,perl=T)
  x
}

read_neoantigen<-function(somatic)
{
  hla=list("I"=list(type=c("A","B","C"),len=c(8:11)),
    "II"=list(type=c("DRB1","DRB3","DRB4","DRB5","DQB1"),len=c(15))) # length of neoantigen for HLA A/B/DRB/DQB
  
  neoantigen=c()
  for (class in names(hla))
  {
    for (type in hla[[class]][["type"]])
    {
      for (len in hla[[class]][["len"]])
      {
        for (allele in c(1,2))
        {
          file=paste(somatic,"_",type,"_",len,"_neoantigen_",allele,".txt",sep="")
          if (!file.exists(file)) {next}
          if (length(readLines(file))==0) # iedb cannot predict
          {
            unlink(file)
            next
          }
          cat(paste("Reading ",file,"\n"))
          tmp=read.table(file,header=T,stringsAsFactors = F,sep="\t")
          neoantigen=rbind(neoantigen,tmp)
          file.remove(file)
        }
      }
    }
  }
  
  neoantigen[!duplicated(neoantigen),]
}

hla<-function()
{
  list(HLA_A=c("NM_001242758","NM_002116","XM_005275331","XM_017030288","XR_430999"),
    HLA_B=c("NM_005514","XM_011514557","XR_926175"),
    HLA_C=c("NM_001243042","NM_002117"),
    HLA_DRB1=c("NM_001243965","NM_002124","XM_011547738"),
    HLA_DRB3=c("NM_022555"),
    HLA_DRB4=c("NM_021983"),
    HLA_DRB5=c("NM_002125"),
    HLA_DQB1=c("NM_001243961","NM_001243962","NM_002123"))
}
