#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if ( length(args) < 2 ) {
    stop("I need a cluster Rds file and a threshold, maybe an output directory\n",
         call.=FALSE)
}
######### This script will cut a tree into groups
suppressWarnings(
    suppressPackageStartupMessages({
        library("cluster")
        library("MCMCpack")
        library("ape")
        library("reshape2")
    })
)

### this is a function to cut hierarchical cluster into groups
### Inputs are dendrogram (tree) and threshold of interest.
prune<-function(dendrogram,threshold) {
    g<-cutree(dendrogram,h=threshold)
    groupnumber<-length(table(g))
    print(paste("cutting hierarchy into",
                groupnumber,
                "groups:",args[1]),
          quote=FALSE)
    engorda<-paste('%0',nchar(groupnumber),"d",sep="")
    b<-noquote(paste("Group",sprintf(engorda,(1:groupnumber)),sep = "-"))
    namelists=lapply(split(names(g),g),sort)
    data.frame(sort(noquote(I(unlist(lapply(namelists,paste,collapse=","))))),
               row.names=noquote(b))
}
### done with the function

##### import tree
tree <- readRDS(args[1])

##### output directory
if ( is.na(args[3]) ) {
    outdir<-"Groups"
} else {
    outdir<-args[3]
}
if (!dir.exists(outdir)) {dir.create(outdir)}

##### output file:
if ( is.na(args[4]) ) {
    groupsfile<-sub("\\.diana-cluster.Rds",
                    paste("-",args[2],sep=""),
                    sub(".+/","",args[1]))
    groupsfile<-paste(groupsfile,"groups",sep=".")
} else {
    groupsfile<-args[4]
}

groupsfile<-paste(outdir,groupsfile,sep="/")

##### prune and write to file
thrd<-as.double(args[2])
groups<-prune(tree,thrd)
write.table(groups,
            file=groupsfile,
            row.names=T,col.names=F,quote=F,append=F)
