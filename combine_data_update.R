

PANDA<-Sys.getenv(c("PANDA"))

norm<-function(n, data){
int0<-as.matrix(data[, (n+1):ncol(data)])
int<-int0
int[is.na(int)]<-0
#low<-apply(int, 1, function(x){return(length(which(x<10)))})
#int<-int[low<(ncol(data)-n),]


totals<-apply(int, 2, sum, na.rm=T)
#max_count<-max(totals)
int[int==0]<-1
int<-t(int)
int_new<-int*(1000000/totals)
int_new<-t(int_new)
#int_new[int_new==0]<-1

colnames(int_new)<-paste(colnames(int_new), "_norm", sep="")
#data<-data[low<(ncol(data)-n),]
data_new<-data.frame(data, int_new)
return(data_new)
}

sample.info<-read.delim("sample.txt", header=T, sep="\t")

samples<-sample.info$Sample
org<-sample.info$Org

###################
#RefSeq Count Summary
###################

allData<-NULL
for(i in 1:nrow(sample.info)){
  sample<-as.character(sample.info[i, "Sample"])
  org<-as.character(sample.info[i, "Org"])
  if(org=="human"){
    sample.data<-read.delim(paste(sample, ".hg18_exon.bwa.align.sam.sum.refseq.txt", sep=""), header=F, sep="\t")
  }else if( org=="mouse"){
    sample.data<-read.delim(paste(sample, ".mm9_exon.bwa.align.sam.sum.refseq.txt",sep=""), header=F, sep="\t") 
   }else if (org=="TB"){
    sample.data<-read.delim(paste(sample, ".h37rv_refMrna.bwa.align.sam.sum.txt",sep=""), header=F, sep="\t")

  }else if (org=="rat"){
    sample.data<-read.delim(paste(sample, ".rn4_exon.bwa.align.sam.sum.refseq.txt",sep=""), header=F, sep="\t")

  }


  colnames(sample.data)<-c("Symbol", "RefSeq", sample)
  sample.data<-sample.data[, 2:3]
  if(is.null(allData)){
    allData<-sample.data
  }else{

    allData<-merge(allData, sample.data, by="RefSeq", all.x=T, all.y=T)
  }

}
if(org=="TB"){
refflat<-read.delim(paste(PANDA, "/db/bwa/TB/h37rv_refMrna/refFlat.txt",sep=""), header=T, sep="\t")
colnames(refflat)<-c("Name", "RefSeq")
allData<-merge(refflat, allData, by="RefSeq")
write.table(allData, "RefSeq_count_summary1.txt", row.names=F, sep="\t", na="", quote=F)

}else{

if(org=="human"){
refflat<-read.delim(paste(PANDA, "/db/bwa/human/hg18_refMrna/refFlat_simple_annot_unique.txt",sep=""), header=T, sep="\t")
}else if(org=="mouse"){

refflat<-read.delim(paste(PANDA, "/db/bwa/mouse/mm9_refMrna/refFlat_simple_annot_unique.txt",sep=""), header=T, sep="\t")

}else if(org=="rat"){

refflat<-read.delim(paste(PANDA, "/db/bwa/rat/rn4_refMrna/refFlat_simple_annot_unique.txt",sep=""), header=T, sep="\t")

}


allData<-norm(1, allData)

allData<-merge(refflat, allData, by="RefSeq")


ts.length<-as.numeric(allData$Length)
samples<-paste(as.character(sample.info[, "Sample"]), "_norm", sep="")
sample.ref<-paste(as.character(sample.info[as.character(sample.info$Type)=="Reference", "Sample"]), "_norm", sep="")
allData.norm<-allData[, samples]

allData.norm<-as.matrix(allData.norm)
allData.norm.length<-allData.norm/ts.length*1000
colnames(allData.norm.length)<-paste(samples, "_PKB", sep="")
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
allData<-data.frame(allData,allData.norm.length, allData.norm)

write.table(allData, "RefSeq_count_summary1.txt", row.names=F, sep="\t", na="", quote=F)



############################
#exon_count summary
############################

allData<-NULL
for(i in 1:nrow(sample.info)){
  sample<-as.character(sample.info[i, "Sample"])
  org<-as.character(sample.info[i, "Org"])
  if(org=="human"){
    sample.data<-read.delim(paste(sample, ".hg18_exon.bwa.align.sam.sum.txt", sep=""), header=F, sep="\t")
  }else if( org=="mouse"){
    sample.data<-read.delim(paste(sample, ".mm9_exon.bwa.align.sam.sum.txt",sep=""), header=F, sep="\t") 
   }else if( org=="rat"){
    sample.data<-read.delim(paste(sample, ".rn4_exon.bwa.align.sam.sum.txt",sep=""), header=F, sep="\t")
   }


  ID<-paste(sample.data[, 3], sample.data[,4], sep=":")
  sample.data<-data.frame(ID=ID, sample.data[,7])
  
  colnames(sample.data)<-c("ID", sample)
 # sample.data<-sample.data[sample.data[,2]>10,]
  if(is.null(allData)){
    allData<-sample.data
  }else{

    allData<-merge(allData, sample.data, by="ID", all.x=T, all.y=T)
  }

}

allData<-norm(1, allData)

if(org=="human"){
exonAnnot<-read.delim(paste(PANDA, "/db/bwa/human/hg18_exon/hg18_exon_sum.txt", sep=""),header=F, sep="\t")
}else if(org=="mouse"){

exonAnnot<-read.delim(paste(PANDA, "/db/bwa/mouse/mm9_exon/mm9_exon_sum.txt",sep=""), header=F, sep="\t")

}else if(org=="rat"){

exonAnnot<-read.delim(paste(PANDA, "/db/bwa/rat/rn4_exon/rn4_exon_sum.txt",sep=""), header=F, sep="\t")

}


colnames(exonAnnot)<-c("Symbol", "RefSeq","ID","Strand", "Range","IDs")
ID<-paste(exonAnnot$ID, exonAnnot$Range, sep=":") 
exonAnnot<-data.frame(ID1=ID,number=1:nrow(exonAnnot), exonAnnot)
allData<-merge(exonAnnot, allData, by.x="ID1", by.y="ID")
allData<-allData[order(as.numeric(allData$number)),]

exonAnnot<-allData[, 1:ncol(exonAnnot)]
exonAnnot<-exonAnnot[, c("Symbol", "RefSeq","ID","Strand", "Range","IDs")]
allData.int<-allData[, (ncol(exonAnnot)+2):ncol(allData)]
allData<-data.frame(exonAnnot, allData.int)

samples<-paste(as.character(sample.info[, "Sample"]), "_norm", sep="")
sample.ref<-paste(as.character(sample.info[as.character(sample.info$Type)=="Reference", "Sample"]), "_norm", sep="")
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



write.table(allData, "Exon_count_summary1.txt", row.names=F, sep="\t", na="", quote=F)


############################
#splice junction summary
############################

allData<-NULL
for(i in 1:nrow(sample.info)){
  sample<-as.character(sample.info[i, "Sample"])
  org<-as.character(sample.info[i, "Org"])
  splice<-as.numeric(sample.info[i, "Splice_len"])
  if(org=="human"){
    sample.data<-read.delim(paste(sample, ".hg18_splice_", splice,".bwa.align.sam.sum.txt", sep=""), header=F, sep="\t")
  }else if( org=="mouse"){
    sample.data<-read.delim(paste(sample, ".mm9_splice_", splice,".bwa.align.sam.sum.txt",sep=""), header=F, sep="\t") 
   }else if( org=="rat"){
    sample.data<-read.delim(paste(sample, ".rn4_splice_", splice,".bwa.align.sam.sum.txt",sep=""), header=F, sep="\t")
   }


  ID<-paste(sample.data[, 3], sample.data[,4], sep=":")
  sample.data<-data.frame(ID=ID, sample.data[,7])
  colnames(sample.data)<-c("ID", sample)
  #sample.data<-sample.data[sample.data[,2]>10,]
  if(is.null(allData)){
    allData<-sample.data
  }else{

    allData<-merge(allData, sample.data, by="ID", all.x=T, all.y=T)
  }

}

allData<-norm(1, allData)

if(org=="human"){
spliceAnnot<-read.delim(paste(PANDA, "/db/bwa/human/hg18_splice_", splice, "/hg18_splice_", splice,"_sum.txt", sep=""), header=F, sep="\t")
}else if(org=="mouse"){

spliceAnnot<-read.delim(paste(PANDA, "/db/bwa/mouse/mm9_splice_", splice, "/mm9_splice_", splice,"_sum.txt", sep=""), header=F, sep="\t")

}else if(org=="rat"){

spliceAnnot<-read.delim(paste(PANDA, "/db/bwa/rat/rn4_splice_", splice, "/rn4_splice_", splice,"_sum.txt", sep=""), header=F, sep="\t")

}



colnames(spliceAnnot)<-c("Symbol", "RefSeq","ID","Strand", "Range","IDs")
ID<-paste(spliceAnnot$ID, spliceAnnot$Range, sep=":")
spliceAnnot<-data.frame(ID1=ID,number=1:nrow(spliceAnnot), spliceAnnot)
allData<-merge(spliceAnnot, allData, by.x="ID1", by.y="ID")
 
allData<-allData[order(as.numeric(allData$number)),]

spliceAnnot<-allData[, 1:ncol(spliceAnnot)]
spliceAnnot<-spliceAnnot[, c("Symbol", "RefSeq","ID","Strand", "Range","IDs")]
allData.int<-allData[, (ncol(spliceAnnot)+2):ncol(allData)]
allData<-data.frame(spliceAnnot, allData.int)

samples<-paste(as.character(sample.info[, "Sample"]), "_norm", sep="")
sample.ref<-paste(as.character(sample.info[as.character(sample.info$Type)=="Reference", "Sample"]), "_norm", sep="")
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



write.table(allData, "Splice_junction_count_summary1.txt", row.names=F, sep="\t", na="", quote=F)


}


################
summary HIV gene expression
################

if(org == "HIV"){

    HIV<-NULL
    for(i in 1:nrow(sample.info)){
      sample<-as.character(sample.info[i, "Sample"])
      HIV.sample<-read.delim(paste(sample, ".HIV_1_Gene.bwa.align.sam.sum.txt", sep=''), header=F, sep="\t")
      colnames(HIV.sample)<-c("ID", "Gene", sample)
      HIV.sample<-HIV.sample[, c("ID", sample)]
      if(i == 1){
          HIV<-HIV.sample
      }else{
          HIV<-merge(HIV, HIV.sample, all.x=T, all.y=T)
      }
    }

    sample.name<-as.character(sample.info[, "Sample"])
    HIV.int<-as.matrix(HIV[, sample.name])
    HIV.int[is.na(HIV.int)]<-1


    RNA_norm<-read.delim("norm.txt", header=T, sep="\t")
    int<-as.matrix(RNA_norm[, sample.name])
    sums<-apply(int, 2, sum, na.rm=T)
    scale<-max(sums)/sums
    HIV.int<-t(HIV.int)*scale
    HIV.int<-t(HIV.int)
    HIV<-data.frame(Gene=HIV$ID, HIV.int)

    write.table(HIV, "HIV_Genes.txt", row.names=F, sep="\t", quote=F, na="")
}









##################
#GO enrichment functions
############################
GO_enrichment<-function(name, org){
    if (org=="human"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/human/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))
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
    colnames(myall)<-c("GO.Term","Symbols", "count", "count.list", "count.term", "count.BP", "pvalue", "adjusted.p", "enrichment.ratio")
    myall_selected<-myall[intersect(which(as.numeric(myall[,"adjusted.p"])<0.05), which(as.numeric(myall[,"enrichment.ratio"])>1)),]
    if(class(myall_selected) == "character"){
        myall_selected<-t(myall_selected)
    }
    if(length(which(as.numeric(myall_selected[,"count"])>2))<2){
        myall_selected<-as.data.frame(t(myall_selected[which(as.numeric(myall_selected[,"count"])>2),]))    
        name<-sub('.txt', '', name)
	write.table(myall_selected, paste(name, " GO.txt", sep=''), row.names=F, sep="\t", na="", quote=F)    
    }else{
        myall_selected<-myall_selected[which(as.numeric(myall_selected[,"count"])>2),]    
        myall_selected<-myall_selected[order(as.numeric(myall_selected[,"adjusted.p"]), decreasing=F), ]
        name<-sub('.txt', '', name)
        write.table(myall_selected, paste(name, " GO.txt", sep=''), row.names=F, sep="\t", na="", quote=F)
    
        ######clustering
        myall<-myall_selected[which(as.numeric(myall_selected[,"count"])>5),]
        if(length(which(as.numeric(myall_selected[,"count"])>5))>=2){
            myall_clustering<-NULL
            cluster_id<-1
            for(j in 1:nrow(myall)){
                if (length(which(myall[j,"GO.Term"]==myall_clustering[, "GO.Term"]))==0){
                    myall_clustering<-rbind(myall_clustering,cbind(cluster.id=cluster_id, t(myall[j,])))             
                    my.symbols<-strsplit(myall[j, "Symbols"], split=", ")[[1]]
                    my.count<-as.numeric(myall[j, "count"])
                    if(j<nrow(myall)){
                          for(k in (j+1):nrow(myall)){
                              if (length(which(myall[k,"GO.Term"]==myall_clustering[, "GO.Term"]))==0){
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
            ragged<-tapply(as.numeric(myall_clustering[,"adjusted.p"]), cluster.factors, median, na.rm=T)
        

            median.p<-NULL
            for(m in 1:length(levels(cluster.factors))){
                median.p<-rbind(median.p, as.data.frame(rep(ragged[m], length(which(as.numeric(myall_clustering[,"cluster.id"])==m)))))
            }
            myall_clustering<-cbind(myall_clustering, median.p)
            colnames(myall_clustering)[11]<-"median.clustering.adjusted.pvalue"
            myall_clustering<-myall_clustering[order(as.numeric(myall_clustering[,"median.clustering.adjusted.pvalue"]), decreasing=F), ]
        
            new.cluster.id<-as.factor(myall_clustering$median.clustering.adjusted.pvalue)
            levels(new.cluster.id)<-paste("C", 1:length(levels(as.factor(myall_clustering$median.clustering.adjusted.pvalue))), sep='')
            myall_clustering$cluster.id<-new.cluster.id
 
            write.table(myall_clustering, paste(name, " GO clustering.txt", sep=''), row.names=F, sep="\t", na="", quote=F)    
        }
    }
    
}

######################################################
GO_enrichment_walk<-function(name, start, window, org){
    if (org=="human"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/human/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))
    }else if(org=="mouse"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/Mus/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))
    }else if(org=="rat"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/Rat/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))    
    }
    all[,"Symbol"]<-toupper(all[,"Symbol"])
    go.gene.num<-length(levels(factor(all$Symbol)))
    gene.list.whole<-read.delim(name, header=T, sep="\t")
    gene.list.whole$Symbol<-gsub("^\\s+|\\s+$", "", gene.list.whole$Symbol)
    gene.list.whole$Symbol<-toupper(gene.list.whole$Symbol)
    gene.list.whole<-data.frame(Symbol=levels(factor(gene.list.whole$Symbol)))
    my.gene.go<-merge(gene.list.whole, all)
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
      myall<-rbind(myall, c(mytmp.name,  toString(as.character(mytmp.gene)), mytmp.gene.num,  my.gene.num, tmp.gene.num,go.gene.num, p, -log10(p), enrichment.ratio))
    }
    p.adj<-p.adjust(myall[,7], method="fdr")
    myall<-cbind(myall[, 1:8], p.adj, myall[, 9])
    colnames(myall)<-c("GO Term","Symbols", "count", "count.list", "count.term", "count.BP", "pvalue", "nlog10.p", "adjusted.p", "enrichment.ratio")
    myall_whole_selected<-data.frame(GO=myall[which(as.numeric(myall[,"count"])>5), "GO Term"])

    sub_all<-merge(myall_whole_selected, all)
    sub_all<-sub_all[,c(2,1)]
    colnames(myall_whole_selected)<-"GO Term"


   
    ##############walk begins
    sub_all[,"Symbol"]<-toupper(sub_all[,"Symbol"])
    go.gene.num<-length(levels(factor(sub_all$Symbol)))    
    gene.list.whole<-read.delim(name, header=T, sep="\t")
    gene.list.whole<-gene.list.whole[order(abs(gene.list.whole$Log2Rat), decreasing=T),]       
    gene.list.whole$Symbol<-gsub("^\\s+|\\s+$", "", gene.list.whole$Symbol)

    if((nrow(gene.list.whole)-window)>0){
    
    if((nrow(gene.list.whole)-window)/start > floor((nrow(gene.list.whole)-window)/start)){
        walk<-floor((nrow(gene.list.whole)-window)/start)+2
    }else{
        walk<-floor((nrow(gene.list.whole)-window)/start)+1
    }
   
    go.merge<-NULL
    go.heatmap<-NULL 
    for(n in 1:walk){
    
        if(n==walk){
            gene.list<-gene.list.whole[(nrow(gene.list.whole)-299):nrow(gene.list.whole), ]
            start_pos<-(nrow(gene.list.whole)-299)
            end_pos<-nrow(gene.list.whole)
        }else{
            gene.list<-gene.list.whole[(start*(n-1)+1):(start*(n-1)+window), ]   
            start_pos<-(start*(n-1)+1)
            end_pos<-(start*(n-1)+window)                      
        }
        
        gene.list$Symbol<-toupper(gene.list$Symbol)
        gene.list<-data.frame(Symbol=levels(factor(gene.list$Symbol)))
        my.gene.go<-merge(gene.list, sub_all)
        GOs<-levels(factor(my.gene.go$GO))
        my.gene.num<-length(levels(factor(my.gene.go$Symbol)))

        myall<-NULL
        for(i in 1:length(GOs)){

          mytmp.gene<-my.gene.go[as.character(my.gene.go$GO)==GOs[i],"Symbol" ]
          mytmp.name<-as.character(my.gene.go[as.character(my.gene.go$GO)==GOs[i],"GO" ])[1]
          mytmp.gene.num<-length(levels(factor(mytmp.gene)))
          tmp.gene<-sub_all[as.character(sub_all$GO)==GOs[i],"Symbol" ]
          tmp.gene.num<-length(levels(factor(tmp.gene)))
          mat<-matrix(c(mytmp.gene.num, (my.gene.num-mytmp.gene.num), tmp.gene.num, (go.gene.num-tmp.gene.num)), nrow=2)
          p<-fisher.test(mat)$p.value
          enrichment.ratio<-(mytmp.gene.num/(my.gene.num-mytmp.gene.num))/(tmp.gene.num/(go.gene.num-tmp.gene.num))
          myall<-rbind(myall, c(mytmp.name, toString(as.character(mytmp.gene)), mytmp.gene.num,  my.gene.num, tmp.gene.num,go.gene.num, p, -log10(p), enrichment.ratio))
        }
        p.adj<-p.adjust(myall[,7], method="fdr")
        myall<-cbind(myall[, 1:8], p.adj, myall[, 9])
        colnames(myall)<-c("GO Term","Symbols", "count", "count.list", "count.term", "count.BP", "pvalue", "nlog10.p", "adjusted.p", "enrichment.ratio")
        myall_selected<-myall[intersect(which(as.numeric(myall[,"pvalue"])<0.05), which(as.numeric(myall[,"enrichment.ratio"])>1)),]
        myall_selected<-myall_selected[order(as.numeric(myall_selected[,"pvalue"]), decreasing=F), ]
   
        write.table(myall_selected, paste(sub('.txt', '', name), " position ", start_pos, "_", end_pos, " GO.txt", sep=''), row.names=F, sep="\t", na="", quote=F)
    
    
    
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
            colnames(myall_clustering)[12]<-"median.clustering.pvalue"
            myall_clustering<-myall_clustering[order(as.numeric(myall_clustering[,"median.clustering.pvalue"]), decreasing=F), ]
        
            new.cluster.id<-as.factor(myall_clustering$median.clustering.pvalue)
            levels(new.cluster.id)<-paste("C", 1:length(levels(as.factor(myall_clustering$median.clustering.pvalue))), sep='')
            myall_clustering$cluster.id<-new.cluster.id
        
            write.table(myall_clustering, paste(sub('.txt', '', name), " position ", start_pos, "_", end_pos, " GO clustering.txt", sep=''), row.names=F, sep="\t", na="", quote=F)
        }else{
            myall_clustering<-myall_clustering[0,]
            colnames(myall_clustering)<-c("cluster.id", colnames(myall_selected), "median.clustering.pvalue")
        }
        
        names(myall_clustering)[c(1, 3:ncol(myall_clustering))]<-paste(names(myall_clustering)[c(1, 3:ncol(myall_clustering))], ".position_", start_pos, "_", end_pos, sep='')
        if (n == 1 ){
            go.merge<-myall_clustering
            go.heatmap<-myall_clustering[, c(2,9)]
            go.gene.symbol<-myall_clustering[, c(2,3)]
        }else{
            go.merge<-merge(go.merge, myall_clustering, by="GO Term", all.x=T, all.y=T)
            go.heatmap<-merge(go.heatmap,myall_clustering[, c(2,9)], by="GO Term", all.x=T, all.y=T)
            go.gene.symbol<-merge(go.gene.symbol,myall_clustering[, c(2,3)], by="GO Term", all.x=T, all.y=T)
        }
        
    }
    xlabel<-colnames(go.heatmap)
    colnames(go.heatmap)<-gsub("nlog10.p.position", "p", xlabel)
    
    xlabel<-colnames(go.heatmap)
    colnames(go.gene.symbol)<-gsub("Symbols.position", "p", xlabel)


    write.table(go.merge, paste(sub('.txt', '', name), " GO merge.txt", sep=''), row.names=F, sep="\t", na="", quote=F)
    write.table(go.heatmap, paste(sub('.txt', '', name), " GO heatmap matrix.txt", sep=''), row.names=F, sep="\t", na="", quote=F)
    write.table(go.gene.symbol, paste(sub('.txt', '', name), " GO merge gene symbol.txt", sep=''), row.names=F, sep="\t", na="", quote=F)

    #cluster
    hm_matrix<-read.delim(paste(sub('.txt', '', name), " GO heatmap matrix.txt", sep=''), header=T, sep="\t")
    gs_matrix<-read.delim(paste(sub('.txt', '', name), " GO merge gene symbol.txt", sep=''), header=T, sep="\t")
    
    int<-hm_matrix[, 2:ncol(hm_matrix)]
    ordered.all<-NULL
    
    for( i in 1:ncol(int)){
    #i<-34
    if(i==1){
    tmp<-hm_matrix[!is.na(int[,i]),]
    
    }
    
    if(i>1){
    
    if(i==2){
    int.prior<-as.matrix(int[, 1:(i-1)], ncol=1)
    }else{
    int.prior<-as.matrix(int[, 1:(i-1)])
    
    }
    
    int.prior.not.na<-apply(int.prior, 1, function(x){return(length(which(!is.na(x)))==0)})
    
    tmp<-hm_matrix[!is.na(int[,i])&int.prior.not.na,]
    
    
    }
    
    
    int.sub<-as.matrix(tmp[, (i+1):ncol(tmp)])
    
    int.sub.weight<-apply(int.sub, 1, function(x){ n<-seq(0,(length(x)-1)); x[!is.na(x)]<-1;x<-x+n/length(x); return(sum(x, na.rm=T))}) 
    tmp<-cbind(tmp, int.sub.weight)
    tmp<-tmp[with(tmp, order(tmp[,"int.sub.weight"], -tmp[,(i+1)])), ]
    tmp<-tmp[, -which(colnames(tmp)%in%c("int.sub.weight"))]
    ordered.all<-rbind(ordered.all, tmp)
    
    
    }
    rownames(gs_matrix)<-gs_matrix$GO.Term
    gs_matrix<-gs_matrix[as.character(ordered.all$GO.Term),] 
    write.table(ordered.all, paste(sub('.txt', '', name), " GO heatmap matrix_ordered.txt", sep=''), na="", row.names=F, sep="\t", quote=F)
    write.table(gs_matrix, paste(sub('.txt', '', name), " GO merge gene symbol.txt", sep=''), na="", row.names=F, sep="\t", quote=F)
    

   #heatmap
   heat_matrix<-read.delim(paste(sub('.txt', '', name), " GO heatmap matrix_ordered.txt", sep=''), header=T, sep="\t")
   
   row.names(heat_matrix)<-heat_matrix$GO.Term
   go.heatmap<-heat_matrix[,2:ncol(heat_matrix)]
   go.heatmap[is.na(go.heatmap)]<-0 

   xlabel<-colnames(go.heatmap)
   colnames(go.heatmap)<-gsub("nlog10.p.position", "p", xlabel)
     
   library("gplots")
   range<-seq(0, 4, by=0.2)
   
   if(gene.list.whole$Log2Rat[1]>0){
       mycol<-colorpanel(n=(length(range)-1), low="grey",   high="red")
   }else{
       mycol<-colorpanel(n=(length(range)-1), low="grey",   high=colors()[258])
   }

   png(paste(sub('.txt', '', name), " GO heatmap.png",sep=''), width=3000, height=3000, res=650)
   heatmap.2(as.matrix(go.heatmap), Rowv=NA, Colv=NA, dendrogram="none", symkey=FALSE, density.info="none", trace="none", cexCol=0.3, cexRow=0.12, breaks=range, col=mycol)
   dev.off()
   
   }
      
     
}



Walk_distance<-function(name, window, org){
    if (org=="human"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/human/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))
    }else if(org=="mouse"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/Mus/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))
    }else if(org=="rat"){
        all<-read.delim(paste(PANDA, "/db/DAVIDKnowledgebase/Rat/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt", sep=""))    
    }
    all[,"Symbol"]<-toupper(all[,"Symbol"])
    go.gene.num<-length(levels(factor(all$Symbol)))
    gene.list.whole<-read.delim(name, header=T, sep="\t")
    gene.list.whole$Symbol<-gsub("^\\s+|\\s+$", "", gene.list.whole$Symbol)
    gene.list.whole$Symbol<-toupper(gene.list.whole$Symbol)
    gene.list.whole<-data.frame(Symbol=levels(factor(gene.list.whole$Symbol)))
    my.gene.go<-merge(gene.list.whole, all)
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
      myall<-rbind(myall, c(mytmp.name,  toString(as.character(mytmp.gene)), mytmp.gene.num,  my.gene.num, tmp.gene.num,go.gene.num, p, -log10(p), enrichment.ratio))
    }
    p.adj<-p.adjust(myall[,7], method="fdr")
    myall<-cbind(myall[, 1:8], p.adj, myall[, 9])
    colnames(myall)<-c("GO Term","Symbols", "count", "count.list", "count.term", "count.BP", "pvalue", "nlog10.p", "adjusted.p", "enrichment.ratio")
    myall_whole_selected<-data.frame(GO=myall[which(as.numeric(myall[,"count"])>5), "GO Term"])

    sub_all<-merge(myall_whole_selected, all)
    sub_all<-sub_all[,c(2,1)]
    colnames(myall_whole_selected)<-"GO Term"


   
    ##############walk begins
    sub_all[,"Symbol"]<-toupper(sub_all[,"Symbol"])
    go.gene.num<-length(levels(factor(sub_all$Symbol)))    
    gene.list.whole<-read.delim(name, header=T, sep="\t")
    gene.list.whole<-gene.list.whole[order(abs(gene.list.whole$Log2Rat), decreasing=T),]       
    gene.list.whole$Symbol<-gsub("^\\s+|\\s+$", "", gene.list.whole$Symbol)

    go.first<-NULL
    n=1
    while(n){
        gene.list<-gene.list.whole[n:(n+window-1), ]   
        
        gene.list$Symbol<-toupper(gene.list$Symbol)
        gene.list<-data.frame(Symbol=levels(factor(gene.list$Symbol)))
        my.gene.go<-merge(gene.list, sub_all)
        GOs<-levels(factor(my.gene.go$GO))
        my.gene.num<-length(levels(factor(my.gene.go$Symbol)))

        myall<-NULL
        for(i in 1:length(GOs)){

          mytmp.gene<-my.gene.go[as.character(my.gene.go$GO)==GOs[i],"Symbol" ]
          mytmp.name<-as.character(my.gene.go[as.character(my.gene.go$GO)==GOs[i],"GO" ])[1]
          mytmp.gene.num<-length(levels(factor(mytmp.gene)))
          tmp.gene<-sub_all[as.character(sub_all$GO)==GOs[i],"Symbol" ]
          tmp.gene.num<-length(levels(factor(tmp.gene)))
          mat<-matrix(c(mytmp.gene.num, (my.gene.num-mytmp.gene.num), tmp.gene.num, (go.gene.num-tmp.gene.num)), nrow=2)
          p<-fisher.test(mat)$p.value
          enrichment.ratio<-(mytmp.gene.num/(my.gene.num-mytmp.gene.num))/(tmp.gene.num/(go.gene.num-tmp.gene.num))
          myall<-rbind(myall, c(mytmp.name,  toString(as.character(mytmp.gene)), mytmp.gene.num,  my.gene.num, tmp.gene.num,go.gene.num, p, -log10(p), enrichment.ratio))
        }
        p.adj<-p.adjust(myall[,7], method="fdr")
	myall<-cbind(myall[, 1:8], p.adj, myall[, 9])
        colnames(myall)<-c("GO Term","Symbols", "count", "count.list", "count.term", "count.BP", "pvalue", "nlog10.p", "adjusted.p", "enrichment.ratio")
        myall_selected<-myall[intersect(which(as.numeric(myall[,"pvalue"])<0.05), which(as.numeric(myall[,"enrichment.ratio"])>1)),]
        myall_selected<-myall_selected[order(as.numeric(myall_selected[,"pvalue"]), decreasing=F), ]
     
    
        myall<-myall_selected[which(as.numeric(myall_selected[,"count"])>5), "GO Term"] 
        #myall<-apply(myall,2,as.character)
        if (n == 1 ){
            go.first<-myall
                n=n+5    
        }else{       
            if(length(intersect(go.first, myall))/min(length(go.first), length(myall))>0.8){
                n=n+5
            }else{break}
        }       
        
    }
    return(n-1)
}



#################################################################################
#DE_analysis
#################################################################################

mytest<-function(x,n){

x1<-x[1:n]
y1<-x[(n+1):length(x)]
x1<-x1[!is.na(x1)]
y1<-y1[!is.na(y1)]
if(length(x1)>2&&length(y1)>2){
return(t.test(x1, y1)$p.value)

}
else{return(1)}
}

get_logratio<-function(x, n){

return( median(x[1:n], na.rm=T)-median(x[(n+1):length(x)], na.rm=T))

}



norm<-function(data, n, p){
int0<-as.matrix(data[, (n+1):ncol(data)])
int<-int0
#q40<-quantile(int, p, na.rm=T)
int[is.na(int)|is.null(int)]<-1
low<-apply(int, 1, function(x){return(length(which(x<p)))})
#int<-int[low<ncol(int),]
totals<-apply(int, 2, sum, na.rm=T)
max_count<-max(totals)
int<-t(int)
int_new<-int*(max_count/totals)
int_new<-t(int_new)
int_new<-log2(int_new)
colnames(int_new)<-paste(colnames(int_new), "_norm", sep="")
data_new<-data.frame(data, int_new)
data_new<-data_new[low<(ncol(data)-n),]
return(data_new)
}

filtering<-function(data, n, minCount){
annot<-data[, 1:n]
int<-as.matrix(data[, (n+1):ncol(data)])
int[is.na(int)]<-0
data<-data.frame(annot, int)
low<-apply(int, 1, function(x){return(length(which(x<minCount)))})

data<-data[low<(ncol(data)-n),]
return(data)
}



####################
# define DE_analysis method
####################


DE_analysis<-function(s1,s2){

if(!file.exists(paste(s2, "_", s1, sep=""))){
 dir.create(paste(s2, "_", s1, sep=""), showWarnings=F)
}

g1<-as.character(samples[as.character(samples$Group)==s1, "Sample"])
g2<-as.character(samples[as.character(samples$Group)==s2, "Sample"])
data1<-data.frame(annot,geneCount[, c(g2, g1)])

data1<-norm(data1, ncol(annot), q50)
n<-length(g2)
int<-as.matrix(data1[, (ncol(annot)+length(g1)+length(g2)+1):ncol(data1)])
p<-apply(int, 1, mytest, n=n)

padj<-p.adjust(p, method="fdr")

Log2Rat<-apply(int, 1, get_logratio, n=n)

data1<-data.frame(data1,  p=p, p.adj=padj, Log2Rat=Log2Rat)
write.table(data1, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
data1.selected<-data1[abs(Log2Rat)>1,]
data1.selected<-data1.selected[order(data1.selected$Log2Rat),]
write.table(data1.selected, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_2fold.txt", sep=""), row.names=F, sep="\t", na="", quote=F)


##############
#DEGseq method
##############

int<-as.matrix(data1[, c(g1, g2)])
int[is.na(int)]<-0
colnames(int)<-c(g1, g2)
data1.simple<-data.frame(RefSeq=data1[, 1], int)
write.table(data1.simple, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_simple.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
geneExpFile<-paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_simple.txt", sep="")
layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE))

par(mar = c(2, 2, 2, 2))
DEGexp2(geneExpFile1 =geneExpFile, expCol1 = g2, groupLabel1 = s2, geneExpFile2 =geneExpFile, expCol2 = g1, groupLabel2 = s1, method = "MARS", outputDir=paste(s2, "_", s1, sep=""), sep="\t")
deg_result<-read.delim(paste(s2, "_", s1, "/output_score.txt", sep=""), header=T, sep="\t")
deg_result<-deg_result[, c(1, 6:10)]

p<-deg_result$p.value
p.adj<-p.adjust(p, method="fdr")
deg_result<-deg_result[p.adj<0.001,]
deg_result<-merge(data1, deg_result, by.x="RefSeq", by.y="GeneNames")
write.table(deg_result, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
deg_result.selected<-deg_result[abs(deg_result$Log2Rat)>1,]
deg_result.selected<-deg_result.selected[order(deg_result.selected$Log2Rat, decreasing=T),]
deg_result.selected_up<-deg_result.selected[deg_result.selected$Log2Rat>0, ]
deg_result.selected_down<-deg_result.selected[deg_result.selected$Log2Rat<0, ]
write.table(deg_result.selected, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_2fold.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
write.table(deg_result.selected_up, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_2fold_upregulated.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
write.table(deg_result.selected_down, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_2fold_downregulated.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
GO_enrichment(paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_2fold_upregulated.txt", sep=""), as.character(org[1]))
GO_enrichment(paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_2fold_downregulated.txt", sep=""), as.character(org[1]))


if (length(g1)==1 & length(g2)==1){    
    deg_result.selected2<-deg_result[abs(deg_result$Log2Rat)>0.58,]
    deg_result.selected2<-deg_result.selected2[order(deg_result.selected2$Log2Rat, decreasing=T),]
    deg_result.selected2_up<-deg_result.selected2[deg_result.selected2$Log2Rat>0, ]
    deg_result.selected2_down<-deg_result.selected2[deg_result.selected2$Log2Rat<0, ]
    deg_result.selected2_down<-deg_result.selected2_down[order(deg_result.selected2_down$Log2Rat, decreasing=F),]
    write.table(deg_result.selected2, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
    write.table(deg_result.selected2_up, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_upregulated.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
    write.table(deg_result.selected2_down, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_downregulated.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
    GO_enrichment(paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_upregulated.txt", sep=""), as.character(org[1]))
    GO_enrichment(paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_downregulated.txt", sep=""), as.character(org[1]))
}else{
    deg_result.selected2<-deg_result[(abs(deg_result$Log2Rat)>0.58 & deg_result$p<0.05),]
    deg_result.selected2<-deg_result.selected2[order(deg_result.selected2$Log2Rat, decreasing=T),]
    deg_result.selected2_up<-deg_result.selected2[deg_result.selected2$Log2Rat>0, ]
    deg_result.selected2_down<-deg_result.selected2[deg_result.selected2$Log2Rat<0, ]
    deg_result.selected2_down<-deg_result.selected2_down[order(deg_result.selected2_down$Log2Rat, decreasing=F),]
    write.table(deg_result.selected2, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_0.05p.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
    write.table(deg_result.selected2_up, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_0.05p_upregulated.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
    write.table(deg_result.selected2_down, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_0.05p_downregulated.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
    GO_enrichment(paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_0.05p_upregulated.txt", sep=""), as.character(org[1]))
    GO_enrichment(paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_1.5fold_0.05p_downregulated.txt", sep=""), as.character(org[1]))

}

deg_result.selected<-deg_result.selected[deg_result.selected$p.adj<0.05,]
write.table(deg_result.selected, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_2fold_0.05adjp.txt", sep=""), row.names=F, sep="\t", na="", quote=F)

}



DE_4fold<-function(s1, s2){


deg_result<-read.delim(paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_2fold.txt", sep=""), header=T)

deg_result<-deg_result[abs(deg_result$Log2Rat)>2,]
deg_result<-deg_result[order(deg_result$Log2Rat, decreasing=T),]
deg_result_selected<-deg_result[, c("RefSeq", "Log2Rat"),]
write.table(deg_result, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_4fold.txt", sep=""), row.names=F, sep="\t", na="", quote=F)
write.table(deg_result_selected, paste(s2, "_", s1, "/geneExp_",s2,"_", s1, "_filtered_DEGseq_4fold_Ingenuity.txt", sep=""), row.names=F, sep="\t", na="", quote=F)

}


library("DESeq")
library(DEGseq)
library(edgeR)
library(baySeq)
geneCount<-read.delim("RefSeq_count_summary1.txt", header=T, sep="\t")

if(org == "human"){
    annot<-geneCount[, 1:16]
    count<-as.matrix(geneCount[, 17:(17+nrow(sample.info)-1)])
}else{
    annot<-geneCount[, 1:17]
    count<-as.matrix(geneCount[, 18:(18+nrow(sample.info)-1)])
}
geneCount<-data.frame(annot, count)
q50<-100

geneCount.norm<-norm(geneCount, ncol(annot), q50)
write.table(geneCount.norm, "norm.txt", row.names=F, sep="\t", na="", quote=F)



###############################
#Comparisons
##################################
samples<-read.delim("sample.txt", header=T, sep="\t")
comparisons<-read.delim("comparisons.txt", header=T, sep="\t")

for (comp in 1:nrow(comparisons)){
    s2<-as.character(comparisons[comp, 1])
    s1<-as.character(comparisons[comp, 2])
    DE_analysis(s1,s2)
    DE_4fold(s1, s2)
}
