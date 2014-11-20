# source("IGew_clinicdata6cl_alldata_NOalvmacro5.R")

rm(list=ls())

source("/Users/Raghu/Documents/asthma2/modified_code/dir_info.R")
source(libfile1)
#source(libfile2)

## get parameters #

# read clinical data #
clinic.file=file.path(data.dir, "dataset.txt")
clinic.data=read.delim(clinic.file)
data=as.matrix(clinic.data[, 4:ncol(clinic.data)])
col.data=colnames(data)
overallres = res.dir

# get phenotype data #
cls=as.numeric(clinic.data[, "Cluster"])
originalphenotypes=cls
allclusters = unique(originalphenotypes)
discret.types = c(2,6)

for (type in discret.types) {
  a = gsub("%",type,"discret%")
  overallres.dir = file.path(overallres,a,"cluster%")

  for (cluster.var in allclusters) { #does each cluster analysis one at a time
    phenotypes = originalphenotypes #original clusters that is going to be modified 
    cat(cluster.var,"\n") #output which cluster is being worked on
    change = which(phenotypes > cluster.var) 
    phenotypes[change] = 0
    change = which(phenotypes < cluster.var) #anything that is not equal to the cluster we are interested in is going to be 0
    phenotypes[change] = 0
    change = which(phenotypes == cluster.var) #
    phenotypes[change] = cluster.var #change this equal to cluster.var when doing second case where discret breaks into 6 clusters/ change to be equal to 1 in first case of discret2
    
    if (type==2){
      clusterlabels=unique(phenotypes) #change to unique(originalphenotypes) for discret6 #change to unique(phenotypes) for discret2
    }
    if (type==6){
      clusterlabels=unique(originalphenotypes)
    }

    n.class.pheno=length(clusterlabels)


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


      res.dir=overallres.dir
      res.dir = gsub("%",cluster.var,res.dir) #change results to put into c
      idx.varig=1:ncol(data) #matrix of the sequence from 1 to the number of columns
      idx.varig.reord=idx.varig[order(infogain.data[idx.varig], decreasing=T)] #sorts the varibles from largest to smallest infogain #the largest are the variables that best edict cluster labels
      infogain.data.all=cbind(col.data[idx.varig.reord], infogain.data[idx.varig.reord]) # has the variables and infogain values from largest to smallest (Best to worst variables for estimating clusters)
      colnames(infogain.data.all)=c("Variable", "InfoGain") #set column headers
      filenameprefix = gsub("%",type,"%Infogain_.txt")
      filenamespecific = gsub("_",cluster.var,filenameprefix)
      infogain.file=filenamespecific#change begining of name of infile into InfogainEW
      infogain.file.fp=file.path(res.dir, infogain.file) #put results file in results directory
      write.table(infogain.data.all, infogain.file.fp, row.names=F,col.names=T, quote=F, sep="\t") #generate table in file
    
      #do the same thing for the normalized inforgain data 
      idx.varsu=1:ncol(data)
      idx.varsu.reord=idx.varsu[order(su.data[idx.varsu], decreasing=T)]
      su.data.all=cbind(col.data[idx.varsu.reord], su.data[idx.varsu.reord])
      colnames(su.data.all)=c("Variable", "SU")
      filenameprefix = gsub("%",type,"%SU_.txt")
      filenamespecific = gsub("_",cluster.var,filenameprefix)
      su.file=filenamespecific#change begining of name of infile into InfogainEW
      su.file.fp=file.path(res.dir, su.file)
      write.table(su.data.all, su.file.fp, row.names=F,col.names=T, quote=F, sep="\t")

      #make file for the matrix of varibales for each patient
      discret.data.all=cbind(col.data[idx.varig.reord], discret.data[idx.varig.reord, ])
      colnames(discret.data.all)=c("Variable", phenotypes)
      filenameprefix = gsub("%",type,"%DiscretEW_.txt")
      filenamespecific = gsub("_",cluster.var,filenameprefix)
      discret.data.file=filenamespecific#change begining of name of infile into InfogainEW
      discret.data.file.fp=file.path(res.dir, discret.data.file)
      write.table(discret.data.all, discret.data.file.fp, row.names=F,col.names=T, quote=F, sep="\t")

  }
}


