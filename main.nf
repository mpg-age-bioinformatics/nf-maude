#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.containers}" == "singularity" ]] ; 

      then

        cd ${params.image_folder}

        if [[ ! -f maude-7aa20cc.sif ]] ;
          then
            singularity pull maude-7a0ffa4.sif docker://index.docker.io/mpgagebioinformatics/maude:7aa20cc
        fi

    fi


    if [[ "${params.containers}" == "docker" ]] ; 

      then

        docker pull mpgagebioinformatics/maude:7aa20cc

    fi

    """

}

process promaude {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    val label
    val paired
    val control
    val treatment
    val paired
    val control_gene
  
  script:
    """
#!/usr/local/bin/Rscript
library(ggplot2)
library(MAUDE)
setwd("${params.output_maude}")
#file path
count_file_path = "${params.ouput_mageck_count}/counts.count.txt"
MAUDE_bin = "${params.facs}"
sgRNA_path = "${params.ctrl_guides}"
# load bin expression data, names of the bins should be the same as in the count table
bins = read.table(MAUDE_bin, row.names=NULL, sep="\\t", header=TRUE)
binStats = bins[!bins\$bin_type=="unsorted",]               #only sorted bins have fractions data
binStats = binStats[!is.na(binStats\$binStartQ),]
binStats\$fraction = binStats\$binEndQ - binStats\$binStartQ; #the fraction of cells captured is the difference in bin start and end percentiles
binStats\$binStartZ = qnorm(binStats\$binStartQ)
binStats\$binEndZ = qnorm(binStats\$binEndQ)
binStats\$Bin = binStats\$bin_type                           #MAUDE needs a column named "Bin"
# load count table, xxxxdouble check if it should be count or normalized count
count_file = read.table(count_file_path, header = T,check.names=FALSE)
count_file = count_file[,names(count_file) %in% c(names(count_file)[1:2],bins[,1])]
neg_sgRNA = read.table(sgRNA_path, header = F)
count_file\$isNontargeting = count_file[,1] %in% neg_sgRNA[,1]     #xxxdouble check sgRNA file format
count_file = unique(count_file) 
count_file_long = melt(count_file, id.vars = c("sgRNA","Gene","isNontargeting"))
count_file_long = merge(count_file_long,bins, by.x="variable", by.y = "samples")
#make long table back to wide table
count_file_good = data.frame(cast(count_file_long, 
                             sgRNA+Gene+isNontargeting+replicates ~ bin_type, value="value"))
#unsorted bins needs to labeled as "unsorted"
sortbins=unique(bins\$bin_type)
sortbins=sortbins[!sortbins%in%"unsorted"]
###guide Level Stats
guideLevelStats = findGuideHitsAllScreens(experiments = unique(count_file_good["replicates"]), 
                                          countDataFrame = count_file_good, binStats = binStats, 
                                          sortBins = sortbins , 
                                          unsortedBin = "unsorted", negativeControl = "isNontargeting")
write.table(guideLevelStats,file = "guideLevelStats.txt", row.names = F, sep = "\\t", quote = FALSE)
guideEffectsByRep = cast(guideLevelStats, 
                         sgRNA + isNontargeting + Gene ~ replicates, value="Z")
guideLevelStats\$replicates=as.factor(guideLevelStats\$replicates)
#plot
p = ggplot(guideLevelStats, aes(x=Z, colour=isNontargeting, linetype=replicates)) + geom_density()+
  theme_classic()+scale_y_continuous(expand=c(0,0)) + geom_vline(xintercept = 0)+
  xlab("Learned guide expression Z score"); 
pdf(paste0("guide-level_Zs.pdf"), width = 4, height = 4)
print(p)
dev.off()
#plot
if(length(bins\$replicates)>1){
  repl=unique(bins\$replicates)
  repl_pairs=combn(repl, 2, simplify = F) 
  
  for (i in 1:length(repl_pairs)) {
    
    x=repl_pairs[[i]][1]
    y=repl_pairs[[i]][2]
    p = ggplot(guideEffectsByRep[!guideEffectsByRep\$isNontargeting,], aes_string(x=x, y=y)) + 
      geom_point(size=0.3) + xlab(paste(x,"Z score")) + ylab(paste(y,"Z score")) + 
      ggtitle(sprintf("r = %f",cor(guideEffectsByRep[!guideEffectsByRep\$isNontargeting,x],
                                   guideEffectsByRep[!guideEffectsByRep\$isNontargeting,y])))+theme_classic(); 
    pdf(paste0("guideEffectsByRep",x,y,".pdf"), width = 4, height = 4)
    print(p)
    dev.off()
  }
}
###Gene Level Stats
GeneLevelStats = getElementwiseStats(experiments = unique(count_file_good["replicates"]), 
                    normNBSummaries = guideLevelStats, negativeControl = "isNontargeting",
                    elementIDs = "Gene")
write.table(GeneLevelStats,file = "GeneLevelStats.txt", row.names = F, sep = "\\t", quote = FALSE)
GeneLevelStatsByRep = cast(GeneLevelStats, Gene ~ replicates, value="meanZ")
if(length(bins\$replicates)>1){
  repl=unique(bins\$replicates)
  repl_pairs=combn(repl, 2, simplify = F) 
  
  for (i in 1:length(repl_pairs)) {
    
    x=repl_pairs[[i]][1]
    y=repl_pairs[[i]][2]
    p = ggplot(GeneLevelStatsByRep, aes_string(x=x, y=y)) + 
      geom_point(size=0.3) + xlab(paste(x,"Z score")) + ylab(paste(y,"Z score")) + 
      ggtitle(sprintf("r = %f",cor(GeneLevelStatsByRep[,x],GeneLevelStatsByRep[,y])))+theme_classic(); 
    pdf(paste0("GeneEffectsByRep",x,y,".pdf"), width = 4, height = 4)
    print(p)
    dev.off()
  }
}
sessionInfo()
    """
}

workflow images {
  main:
    get_images()
}


workflow {
    if ( 'facs' in params.keySet()  ) {
      if ( ! file("${params.output_maude}").isDirectory() ) {
        file("${params.output_maude}").mkdirs()
      }

      promaude( )
    }

}