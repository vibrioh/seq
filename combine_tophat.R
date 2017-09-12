
PANDA<-Sys.getenv(c("PANDA"))


sample.info<-read.delim("sample.txt", header=T, sep="\t")

samples<-as.character(sample.info$Sample)
org<-as.character(sample.info$Org)[1]
###################
#RefSeq Count Summary
###################

allData<-NULL
for(sample in samples){
    sample.data<-read.delim(paste(sample, "/genes.expr", sep=""), header=T, sep="\t")

  
  sample.data<-sample.data[,c("gene_id", "FPKM") ]
  if(is.null(allData)){
    allData<-sample.data
  }else{
    allData<-merge(allData, sample.data, by="gene_id", all.x=T, all.y=T)
  }
}
colnames(allData)<-c("RefSeq",as.character(samples))

if(org=="human"){
refflat<-read.delim(paste(PANDA, "/db/bwa/human/hg18_refMrna/refFlat_simple_annot_unique.txt", sep=""), header=T, sep="\t")
}else if(org=="mouse"){

refflat<-read.delim(paste(PANDA, "/db/bwa/mouse/mm9_refMrna/refFlat_simple_annot_unique.txt", sep=""), header=T, sep="\t")

}else if(org=="rat"){

refflat<-read.delim(paste(PANDA, "/db/bwa/rat/rn49_refMrna/refFlat_simple_annot_unique.txt", sep=""), header=T, sep="\t")
}




allData<-merge(refflat, allData, by="RefSeq")
sample.ref<-as.character(sample.info[as.character(sample.info$Type)=="Reference", "Sample"])
allData.norm<-allData[, samples]

allData.norm<-as.matrix(allData.norm)
allData.norm<-log2(allData.norm)

if(length(sample.ref)==1){
allData.ref<-allData.norm[, sample.ref]
allData.norm<-allData.norm-allData.ref
}else{
allData.ref<-allData.norm[, sample.ref]
allData.ref<-apply(allData.ref, 1, median, na.rm=T)
allData.norm<-allData.norm-allData.ref

}
colnames(allData.norm)<-paste("Log2Rat_", samples, sep="")
allData<-data.frame(allData, allData.norm)
write.table(allData, "Gene_count_summary1.txt", row.names=F, sep="\t", na="", quote=F)


################################################
# summary for transcript-level data
################################################
allData<-NULL
for(sample in samples){
    sample.data<-read.delim(paste(sample, "/transcripts.expr", sep=""), header=T, sep="\t")


  sample.data<-sample.data[,c("trans_id", "FPKM") ]
  if(is.null(allData)){
    allData<-sample.data
  }else{
    allData<-merge(allData, sample.data, by="trans_id", all.x=T, all.y=T)
  }
}
colnames(allData)<-c("RefSeq", samples)

if(org=="human"){
refflat<-read.delim("/data/pipeline_in/Genomes/bwa/human/hg18_refMrna/refFlat_simple_annot_unique.txt", header=T, sep="\t")
}else if(org=="mouse"){

refflat<-read.delim("/data/pipeline_in/Genomes/bwa/mouse/mm9_refMrna/refFlat_simple_annot_unique.txt", header=T, sep="\t")

}else if(org=="rat"){

refflat<-read.delim("/data/pipeline_in/Genomes/bwa/rat/rn4_refMrna/refFlat_simple_annot_unique.txt", header=T, sep="\t")

}

allData<-merge(refflat, allData, by="RefSeq")

}
write.table(allData, "Transcript_count_summary1.txt", row.names=F, sep="\t", na="", quote=F)


#############################
#GO enrichment 
#################################
GO_enrichment<-function(name, org){
    if (org=="human"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/human/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt",sep=""))
    }else if(org=="mouse"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/Mus/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))
    }else if(org=="rat"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/Rat/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))    
    }
    
    all[,"Symbol"]<-toupper(all[,"Symbol"])
    go.gene.num<-length(levels(factor(all$Symbol)))
    gene.list<-read.delim(name, header=T, sep="\t")
    gene.list$Symbol<-toupper(gene.list$Symbol)
    gene.list$Symbol<-gsub("^\\s+|\\s+$", "", gene.list$Symbol)
    gene.list<-data.frame(Symbol=levels(factor(gene.list$Symbol)))
    my.gene.go<-merge(gene.list, all)
    GOs<-levels(factor(my.gene.go$GO))
    my.gene.num<-length(levels(factor(my.gene.go$Symbol)))

    myall<-NULL
    for(i in 1:length(GOs)){

      mytmp.gene<-my.gene.go[as.character(my.gene.go$GO)==GOs[i],"Symbol" ]
      mytmp.name<-as.character(my.gene.go[as.character(my.gene.go$GO)==GOs[i],"GO" ])[1]
      mytmp.gene.num<-length(levels(factor(mytmp.gene)))
      tmp.gene<-all[as.character(all$GO)==GOs[i],"Symbol" ]
      tmp.gene.num<-length(levels(factor(tmp.gene)))
      mat<-matrix(c(mytmp.gene.num, (my.gene.num-mytmp.gene.num), tmp.gene.num, (go.gene.num-tmp.gene.num)), nrow=2)
      p<-fisher.test(mat)$p.value
      enrichment.ratio<-(mytmp.gene.num/(my.gene.num-mytmp.gene.num))/(tmp.gene.num/(go.gene.num-tmp.gene.num))
      myall<-rbind(myall, c(mytmp.name,  toString(as.character(mytmp.gene)), mytmp.gene.num,  my.gene.num, tmp.gene.num,go.gene.num, p, enrichment.ratio))
    }
    p.adj<-p.adjust(myall[,7], method="fdr")
    myall<-cbind(myall[, 1:7], p.adj, myall[, 8])
    colnames(myall)<-c("GO Term","Symbols", "count", "count.list", "count.term", "count.BP", "pvalue", "adjusted.p", "enrichment.ratio")
    myall_selected<-myall[intersect(which(as.numeric(myall[,"pvalue"])<0.05), which(as.numeric(myall[,"enrichment.ratio"])>1)),]
    if(length(which(as.numeric(myall_selected[,"count"])>2))<2){
        myall_selected<-as.data.frame(t(myall_selected[which(as.numeric(myall_selected[,"count"])>2),]))    
        name<-sub('.txt', '', name)
	write.table(myall_selected, paste(name, " GO.txt", sep=''), row.names=F, sep="\t", na="", quote=F)    
    }else{
        myall_selected<-myall_selected[which(as.numeric(myall_selected[,"count"])>2),]    
        myall_selected<-myall_selected[order(as.numeric(myall_selected[,"pvalue"]), decreasing=F), ]
        name<-sub('.txt', '', name)
        write.table(myall_selected, paste(name, " GO.txt", sep=''), row.names=F, sep="\t", na="", quote=F)
    
        ######clustering
        myall<-myall_selected[which(as.numeric(myall_selected[,"count"])>5),]
        if(length(which(as.numeric(myall_selected[,"count"])>5))>=2){
            myall_clustering<-NULL
            cluster_id<-1
            for(j in 1:nrow(myall)){
                if (length(which(myall[j,"GO Term"]==myall_clustering[, "GO Term"]))==0){
                    myall_clustering<-rbind(myall_clustering,cbind(cluster.id=cluster_id, t(myall[j,])))             
                    my.symbols<-strsplit(myall[j, "Symbols"], split=", ")[[1]]
                    my.count<-as.numeric(myall[j, "count"])
                    if(j<nrow(myall)){
                          for(k in (j+1):nrow(myall)){
                              if (length(which(myall[k,"GO Term"]==myall_clustering[, "GO Term"]))==0){
                                 its.symbols<-strsplit(myall[k, "Symbols"], split=", ")[[1]]
                                 its.count<-as.numeric(myall[k, "count"])
                                 if(length(intersect(my.symbols, its.symbols))/min(its.count, my.count)>0.5){
                                     myall_clustering<-rbind(myall_clustering, cbind(cluster.id=cluster_id, t(myall[k,])))
                                 } 
                              }else{next}
                          }
                          cluster_id<-cluster_id+1
                      }
                 }else{next}
            }
    
            cluster.factors<-factor(as.numeric(myall_clustering[,"cluster.id"]))
            ragged<-tapply(as.numeric(myall_clustering[,"pvalue"]), cluster.factors, median, na.rm=T)
        

            median.p<-NULL
            for(m in 1:length(levels(cluster.factors))){
                median.p<-rbind(median.p, as.data.frame(rep(ragged[m], length(which(as.numeric(myall_clustering[,"cluster.id"])==m)))))
            }
            myall_clustering<-cbind(myall_clustering, median.p)
            colnames(myall_clustering)[11]<-"median.clustering.pvalue"
            myall_clustering<-myall_clustering[order(as.numeric(myall_clustering[,"median.clustering.pvalue"]), decreasing=F), ]
        
            new.cluster.id<-as.factor(myall_clustering$median.clustering.pvalue)
            levels(new.cluster.id)<-paste("C", 1:length(levels(as.factor(myall_clustering$median.clustering.pvalue))), sep='')
            myall_clustering$cluster.id<-new.cluster.id
 
            write.table(myall_clustering, paste(name, " GO clustering.txt", sep=''), row.names=F, sep="\t", na="", quote=F)    
        }
    }
    
}



comparisons<-read.delim("comparisons.txt", header=T, sep="\t")

for (comp in 1:nrow(comparisons)){

    s2<-as.character(comparisons[comp, 1])
    s1<-as.character(comparisons[comp, 2])
    exp_diff<-read.delim(paste("known_Transcript_diff_", s2, "_", s1, "/gene_exp.diff", sep=''), header=T, sep="\t")
    exp_diff<-exp_diff[(exp_diff$q_value<0.05 & abs(exp_diff$log2.fold_change.)>0.58), ]
    exp_symbol<-merge(exp_diff, refflat, by.x="gene_id", by.y="RefSeq")[, 1:(ncol(exp_diff)+1)]
    exp_symbol_up<-exp_symbol[exp_symbol$log2.fold_change.>0, ]
    exp_symbol_up<-exp_symbol_up[order(exp_symbol_up$log2.fold_change., decreasing=T), ]
    exp_symbol_down<-exp_symbol[exp_symbol$log2.fold_change.<0, ]
    exp_symbol_down<-exp_symbol_down[order(exp_symbol_down$log2.fold_change., decreasing=F), ]

    exp_symbol_2f_up<-exp_symbol[exp_symbol$log2.fold_change.>1, ]
    exp_symbol_2f_up<-exp_symbol_2f_up[order(exp_symbol_2f_up$log2.fold_change., decreasing=T), ]
    exp_symbol_2f_down<-exp_symbol[exp_symbol$log2.fold_change.<(-1), ]
    exp_symbol_2f_down<-exp_symbol_2f_down[order(exp_symbol_2f_down$log2.fold_change., decreasing=F), ]

    
    name_up<-paste("known_Transcript_diff_", s2, "_", s1, "/gene_exp_", s2, "_", s1, "_selected_1.5fold_symbol_upregulated.diff", sep='')
    write.table(exp_symbol_up, name_up, row.names=F, sep="\t", na="", quote=F)    
    name_down<-paste("known_Transcript_diff_", s2, "_", s1, "/gene_exp_", s2, "_", s1, "_selected_1.5fold_symbol_downregulated.diff", sep='')
    write.table(exp_symbol_down, name_down, row.names=F, sep="\t", na="", quote=F)  

    name_2f_up<-paste("known_Transcript_diff_", s2, "_", s1, "/gene_exp_", s2, "_", s1, "_selected_2fold_symbol_upregulated.diff", sep='')
    write.table(exp_symbol_2f_up, name_2f_up, row.names=F, sep="\t", na="", quote=F)    
    name_2f_down<-paste("known_Transcript_diff_", s2, "_", s1, "/gene_exp_", s2, "_", s1, "_selected_2fold_symbol_downregulated.diff", sep='')
    write.table(exp_symbol_2f_down, name_2f_down, row.names=F, sep="\t", na="", quote=F)  
  

    GO_enrichment(name_up, org)    
    GO_enrichment(name_down, org)    

}








