# source("IGew_clinicdata6cl_alldata_NOalvmacro5.R")

rm(list=ls())

source("/Users/Raghu/Documents/asthma2/modified_code/dir_info.R")
source(libfile1)
#source(libfile2)

## get parameters #

# read clinical data #
clinic.file=file.path(data.dir, "dataset.txt")
clinic.data=read.delim(clinic.file)
data=apply(clinic.data[, 4:ncol(clinic.data)], 2, as.numeric)
col.data=colnames(data)

# get phenotype data #
cls=as.numeric(clinic.data[, "Cluster"])
phenotypes=cls
n.class.pheno=length(unique(phenotypes))

  infogain.data=discret.data=su.data=NULL
  for (i.var in 1:ncol(data)) {
  # i.var=1

    cat(i.var, "\n")

    data1=NULL
    cl.data1=length(unique(data[,i.var]))
    if (cl.data1>n.class.pheno) {
      data1=discret.equalwidth(data[,i.var], n.class.pheno)
    } else {
      data1=data[,i.var]
    }

    infogain.data=c(infogain.data, infogain(phenotypes, data1))
    su.data=c(su.data, infogain.normalized(phenotypes, data1))
    discret.data=rbind(discret.data, data1) 
  }

  idx.varig=1:ncol(data)
  idx.varig.reord=idx.varig[order(infogain.data[idx.varig], decreasing=T)]
  infogain.data.all=cbind(col.data[idx.varig.reord], infogain.data[idx.varig.reord])
  colnames(infogain.data.all)=c("Variable", "InfoGain")
  infogain.file.fp=file.path(res.dir, 'InfoGain.txt')
  write.table(infogain.data.all, infogain.file.fp, row.names=F,col.names=T, quote=F, sep="\t")

  idx.varsu=1:ncol(data)
  idx.varsu.reord=idx.varsu[order(su.data[idx.varsu], decreasing=T)]
  su.data.all=cbind(col.data[idx.varsu.reord], su.data[idx.varsu.reord])
  colnames(su.data.all)=c("Variable", "SU")
  su.file.fp=file.path(res.dir,"SuEW.txt")
  write.table(su.data.all, su.file.fp, row.names=F,col.names=T, quote=F, sep="\t")

  discret.data.all=cbind(col.data[idx.varig.reord], discret.data[idx.varig.reord, ])
  colnames(discret.data.all)=c("Variable", phenotypes)
  discret.data.file.fp=file.path(res.dir, "discretEW.txt")
  write.table(discret.data.all, discret.data.file.fp, row.names=F,col.names=T, quote=F, sep="\t")




