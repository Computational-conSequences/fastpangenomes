#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if ( length(args) < 1 ) {
    msg = paste("\n I need a file with distances by many methods, and, perhaps",
                " whether you want the cutoff point displayed (y/n), def: n",
               " ex:\n  findCutOffs.R SameTaxaTables/Enterobacterales.dist.bz2 y\n",sep="\n")
    stop(msg, call.=FALSE)
}
######### This script will cut a tree into groups

suppressWarnings(
    suppressPackageStartupMessages({
        library("tidyverse")
        library("cutpointr")
    })
)

## for plots:
#axisTitle <- 30
#axisTics  <- 25
#axisMargin <- 12
axisTitle <- 33
axisTics  <- 27.5
axisMargin <- 13.2
plotFonts = theme(
    plot.title   = element_text(size = 40),
    axis.title.x = element_text(size = axisTitle,margin=margin(t=axisMargin)),
    axis.text.x  = element_text(size = axisTics),
    axis.title.y = element_text(size = axisTitle,margin=margin(r=axisMargin)),
    axis.text.y  = element_text(size = axisTics)
    )


fullCuts<-function(chidos,gachos,encabezado) {
    print(paste("calculating",chidos,"vs",gachos,encabezado,"cutoffs"),quote=F)

    tablita<-subset(tablota, Same == chidos | Same == gachos)
    medida<-get(encabezado)

    muchos=multi_cutpointr(tablita,
                           x=medidas,
                           class=Same,
                           direction = "<=",
                           pos_class=chidos,
                           neg_class=gachos,
                           method = maximize_metric,
                           metric = medida,
                           break_ties = max
                           )
    
    cuttable<-select(muchos,
                     optimal_cutpoint,
                     all_of(encabezado),
                     acc,
                     sensitivity,
                     specificity,
                     AUC,
                     predictor
                     )
    
    dashedtbl<-cuttable %>% mutate_if(is.character,
                                      str_replace_all,
                                      pattern = "\\.",
                                      replacement="-")
    
    roundedtbl<-dashedtbl %>% mutate_if(is.numeric, sprintf, fmt = "%.6f")
    
    outfile = paste(paste(chidos,gachos,sep="-"),
                    encabezado,"cutoffs",sep=".")
    if (!dir.exists(outdir)) {dir.create(outdir)}
    outfile<-paste(outdir,outfile,sep="/")
    write.table(roundedtbl,
                file=outfile,
                quote=FALSE,
                sep="\t",
                row.names=FALSE)

}

singleGraph<-function(ms,chidos,gachos) {
    print(paste("plotting",ms,"ROC:",chidos,"vs",gachos),quote=F)
    tablita<-subset(tablota, Same == chidos | Same == gachos)
    cuts=cutpointr(tablita,
                   x=!!ms,
                   class=Same,
                   direction = "<=",
                   pos_class=chidos,
                   neg_class=gachos,
                   method = maximize_metric,
                   metric = F1_score,
                   break_ties = max
                   )
    AUC=paste("AUC:",sprintf("%.3f",cuts$AUC))
    OCP=paste("OCP:",sprintf("%.3f",cuts$optimal_cutpoint))
    if ( length(args) > 1 && args[2] == "y" ) {
        statlabel=paste(AUC,OCP,sep="\n")
        mymin = 0.35
        mymax = 0.65
    } else {
        statlabel=AUC
        mymin = 0.40
        mymax = 0.60
    }
    rocFile<-gsub("\\.","-",ms)
    if ( str_detect(ms,"ANI") ) {
        rocTitle<-ms
    } else if ( str_detect(ms,"dashing") ) {
        rocTitle<-"dashing"
    } else {
        rocTitle<-ms
    }
    outfile<-paste(rocFile,"pdf",sep=".")
    outfile<-paste(outdir,outfile,sep="/")
    roc<-plot_roc(cuts,display_cutpoint=F) +
        ggtitle(rocTitle) +
        geom_line(linewidth=2) +
        theme_light() +
        annotate("rect", xmin = 0.2, xmax = 0.8, ymin = mymin, ymax = mymax,
                 fill="white",colour="black") +
        annotate("text", x = 0.5, y = 0.5, label = statlabel, size = 14) +
        plotFonts
    ggsave(outfile, width=8, height=7)
}

outdir<-"CutoffTables"
print(paste("reading",args[1]),quote=F)
tablota<-read.table(args[1],sep="\t",head=T)
medidas<-grep("ANI|Signature|mash|dashing|fna",colnames(tablota),value=T,perl=T)
metricas<-c("accuracy","F1_score","sum_sens_spec")
for ( metrica in metricas ) {
    fullCuts("Species","Genus",metrica)
    if( any(tablota$Same=="Family") ) {
        fullCuts("Species","Family",metrica)
        fullCuts("Genus","Family",metrica)
    }
}

if ( length(args) > 1 && args[2] == "y" ) {
    print("displaying AUC and optimal cutoff in ROC plot",quote=F)
} else {
    print("only AUC will be displayed in ROC plot",quote=F)
}
for ( ms in medidas ) {
    singleGraph(ms,"Species","Genus")
}
