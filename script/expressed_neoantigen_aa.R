######  read arguments  #################

# first argument: somatic mutation file 
# second argument: output directory
# third argument: path to the script folder

args=commandArgs(trailingOnly = TRUE) 
somatic=args[1]
output=args[2]
path=args[3]

#######  read neoantigen binding strength  ##################

source(paste(path,"/expressed_neoantigen_function.R",sep=""))
neoantigen=read_neoantigen(somatic)

#####  write results  ########################

write.table(neoantigen[,c("peptide","wildtype_match",
  "allele","percentile_rank")],file=paste(somatic,"_final.txt",sep=""),
  sep="\t",quote=F,row.names = F)
