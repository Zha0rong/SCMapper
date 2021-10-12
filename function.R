SCMapper=function(Object_A,Object_B,
         annotation_name_A='seurat_clusters',annotation_name_B='seurat_clusters',
         Name_of_object_A,Name_of_object_B,
         numbersofPC=30,roundsofrandomization=1000,Normalization_method='LogNormalize',harmonize=T,p.val.threshold=0.05,prediction_plot_ncol=2,number_of_core=4) {
  require(cluster)
  require(reshape)
  require(ggplot2)
  require(Seurat)
  require(harmony)
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(lsa)

  
  Object_A@meta.data[['scmapper']]=paste(Name_of_object_A,
                                         '_',Object_A@meta.data[[annotation_name_A]],sep = '')
  Object_A@meta.data[['object']]=Name_of_object_A
  Object_B@meta.data[['scmapper']]=paste(Name_of_object_B,
                                         '_',Object_B@meta.data[[annotation_name_B]],sep = '')
  Object_B@meta.data[['object']]=Name_of_object_B
  
  Merged=merge(Object_A,Object_B)
  Merged[['percent.mt']]=NULL
  DefaultAssay(Merged)='RNA'
  if (Normalization_method=='SCT'){
    if ('percent.mt'%in%names(Merged@meta.data)) {
      Merged=SCTransform(Merged,conserve.memory = T,vars.to.regress = 'percent.mt',verbose = F)
      Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)
    }
    else {
      Merged=SCTransform(Merged,conserve.memory = T)
      Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)
    }
    print('Finished SCT Normalization')
  }
  if (Normalization_method=='LogNormalize'){
    Merged=NormalizeData(Merged,verbose = F)
    Merged=FindVariableFeatures(Merged,verbose = F)
    Merged=ScaleData(Merged,verbose = F)
    Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)
    print('Finished Log Normalization')
  }
  if (Normalization_method=='None'){
    Merged=FindVariableFeatures(Merged,verbose = F)
    Merged=ScaleData(Merged,verbose = F)
    Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)

    
    
  }
  
  if (harmonize) {
    if (Normalization_method=='SCT') {
      Merged=RunHarmony(Merged,group.by.vars = 'object',max.iter.harmony = 1000,
                        max.iter.cluster = 1000,verbose = F,assay.use='SCT')
      PCA.space=Merged@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
      
    }
    else {Merged=RunHarmony(Merged,group.by.vars = 'object',max.iter.harmony = 1000,
                      max.iter.cluster = 1000,verbose = F)
    PCA.space=Merged@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
    }
  }
  else {
    PCA.space=Merged@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
  }
  mapping.table=data.frame(cell=rownames(Merged@meta.data),cluster=Merged@meta.data[['scmapper']])
  Object_A_annotation=mapping.table$cluster[grepl(Name_of_object_A,mapping.table$cluster)]
  Object_A_annotation=as.character(unique(Object_A_annotation))
  Object_B_annotation=mapping.table$cluster[grepl(Name_of_object_B,mapping.table$cluster)]
  Object_B_annotation=as.character(unique(Object_B_annotation))
  zscore.matrix=matrix(rep(0,length(Object_A_annotation)*length(Object_B_annotation))
                       ,nrow = length(Object_A_annotation)
                       ,dimnames = list(Object_A_annotation,
                                        Object_B_annotation))
  zscore.matrix=zscore.matrix[order((rownames(zscore.matrix))),
                              order((colnames(zscore.matrix)))]
  pval.matrix=matrix(rep(0,length(Object_A_annotation)*length(Object_B_annotation))
                     ,nrow = length(Object_A_annotation)
                     ,dimnames = list(Object_A_annotation,
                                      Object_B_annotation))
  pval.matrix=pval.matrix[order((rownames(pval.matrix))),
                          order((colnames(pval.matrix)))]

  Detail_list=list()
  for (i in rownames(zscore.matrix)) {
    Detail_list[[i]]$Aggregated.randomization.results=data.frame(Distance=0,group='',X='',Y='')
  }
  meta.data=Merged@meta.data
  rm(Merged)
  ####Parallelization####
  cl <- makeCluster(number_of_core)
  registerDoSNOW(cl)
  
  Test_Candidates=paste(rep(rownames(zscore.matrix), each = length(colnames(zscore.matrix))), colnames(zscore.matrix), sep = "Object_42")
  iterations <- length(Test_Candidates)
  
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  
  opts <- list(progress = progress)
  
  x=foreach(i=Test_Candidates,.packages=c('ggplot2','lsa'), .options.snow = opts) %dopar% {
    a=as.character(unlist(strsplit(i,split = "Object_42")))[1]
    b=as.character(unlist(strsplit(i,split = "Object_42")))[2]
    
    euclidean <- function(x, y) sqrt(sum((x - y)^2))
    AnnotationA=mapping.table$cell[mapping.table$cluster==a]
    
    AnnotationA.centroid=apply(PCA.space[AnnotationA,],2,median)

    Aggregated.randomization.results=data.frame(Distance=0,group='',X='',Y='')
    AnnotationB=mapping.table$cell[mapping.table$cluster==b]
    
    AnnotationB.centroid=apply(PCA.space[AnnotationB,],2,median)
    
    true.distance=c()
    bootstrapping.distance=c()
    bootstrapping.sampling.distance=c()
    for (k in 1:roundsofrandomization) {
      sub_sampling_size=runif(1,min=0.5,max = 0.75)
      if (sub_sampling_size*length(AnnotationA)<=2) {
        sub_sampling_1=PCA.space[AnnotationA,]
        
        temp=PCA.space[sample(AnnotationA,size = 2,replace = T),]
        rownames(temp)=paste0(rownames(temp),'bootstrap')
        sub_sampling_1=rbind(sub_sampling_1,temp)
        
        sub_sampling_2=PCA.space[AnnotationB,]
        temp=PCA.space[sample(AnnotationB,size = 2,replace = T),]
        rownames(temp)=paste0(rownames(temp),'bootstrap')
        sub_sampling_2=rbind(sub_sampling_2,temp)
        
        sub_sampling_1=apply(sub_sampling_1,2,median)
        sub_sampling_2=apply(sub_sampling_2,2,median)
      }
      else {
        
        sub_sampling_1=sample(AnnotationA,size = sub_sampling_size*length(AnnotationA) ,replace = T)
        sub_sampling_2=sample(AnnotationB,size = sub_sampling_size*length(AnnotationB) ,replace = T)
        sub_sampling_1=apply(PCA.space[sub_sampling_1,],2,median)
        sub_sampling_2=apply(PCA.space[sub_sampling_2,],2,median)
      }

      sub_sampling_distance=euclidean(sub_sampling_1,sub_sampling_2)
      true.distance=c(true.distance,sub_sampling_distance)
      
      sampling_pool=rownames(PCA.space)
      sampling_pool.A=sampling_pool[sampling_pool%in%rownames(meta.data)[meta.data$object==Name_of_object_A]]
      sampling_pool.B=sampling_pool[sampling_pool%in%rownames(meta.data)[meta.data$object==Name_of_object_B]]
      
      sampling.subset.A=sample(sampling_pool.A,size = (length(AnnotationA)+length(AnnotationB))/2,replace = T)
      sampling.subset.B=sample(sampling_pool.B,size = (length(AnnotationA)+length(AnnotationB))/2,replace = T)
      sampling.subset.A.centroid=apply(PCA.space[sampling.subset.A,],2,median)
      
      sampling.subset.B.centroid=apply(PCA.space[sampling.subset.B,],2,median)

      sampling.distance=((euclidean(sampling.subset.B.centroid,AnnotationA.centroid))+
                           (euclidean(sampling.subset.A.centroid,AnnotationB.centroid)))/2
      sampling.distance=(euclidean(sampling.subset.B.centroid,AnnotationA.centroid))
      bootstrapping.sampling.distance=c(bootstrapping.sampling.distance,sampling.distance)
    }
    
    Comparison.results=data.frame(Distance=bootstrapping.sampling.distance,group=rep('Randomized',length(bootstrapping.sampling.distance)))
    Comparison.results=rbind(Comparison.results,data.frame(Distance=true.distance,group=rep(paste(a,b),length(true.distance))))
    Comparison.results$X=a
    Comparison.results$Y=b
    Aggregated.randomization.results=rbind(Aggregated.randomization.results,
                                                            Comparison.results)
    wilcox.results=wilcox.test(log10(true.distance+1),log10(bootstrapping.sampling.distance+1),alternative = 'less')
    zscore=median(true.distance)-median(bootstrapping.sampling.distance)
    zscore.matrix[a,b]=mean(true.distance)-mean(bootstrapping.sampling.distance)
    pval.matrix[a,b]=wilcox.results$p.value
    Aggregated.randomization.results=Aggregated.randomization.results[-1,]
    Aggregated.randomization.results=Aggregated.randomization.results
    
    Aggregated.Randomization.Results.plot=ggplot(Aggregated.randomization.results[Aggregated.randomization.results$group!='Randomized',],
                                                 aes(x=Distance,fill=Y))+geom_density(alpha=0.75)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
    Aggregated.randomization.results.plot=Aggregated.Randomization.Results.plot
    
    
    Results=list(a,b,zscore,wilcox.results,Aggregated.randomization.results,Aggregated.Randomization.Results.plot)
    return(Results)
  }
  close(pb)
  stopCluster(cl)
  
  for (i in 1:length(x)) {
    temp=x[[i]]
    zscore.matrix[temp[[1]],temp[[2]]]=temp[[3]]
    pval.matrix[temp[[1]],temp[[2]]]=temp[[4]]$p.value
    Detail_list[[temp[[1]]]]$Aggregated.randomization.results=rbind(Detail_list[[temp[[1]]]]$Aggregated.randomization.results,
                                                                    temp[[5]])

  }
  for (i in rownames(zscore.matrix)) {
    Detail_list[[i]]$Aggregated.randomization.results=Detail_list[[i]]$Aggregated.randomization.results[-1,]
    Aggregated.randomization.results=Detail_list[[i]]$Aggregated.randomization.results
    Aggregated.Randomization.Results.plot=ggplot(Aggregated.randomization.results[Aggregated.randomization.results$group!='Randomized',],
                                                 aes(x=Distance,fill=Y))+geom_density(alpha=0.75)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
    
    Detail_list[[i]]$Aggregated.Randomization.Results.plot=Aggregated.Randomization.Results.plot

  }
  rm(x)
  gc()

  ####
  #for (i in 1:nrow(zscore.matrix)) {
  #  AnnotationA=mapping.table$cell[mapping.table$cluster==rownames(zscore.matrix)[i]]

#    AnnotationA.centroid=apply(PCA.space[AnnotationA,],2,median)
#    Detail_list[[rownames(zscore.matrix)[i]]]=list()
#    Detail_list[[rownames(zscore.matrix)[i]]]$Aggregated.randomization.results=data.frame(Distance=0,group='',X='',Y='')
#    for (j in 1:ncol(zscore.matrix)) {
#      AnnotationB=mapping.table$cell[mapping.table$cluster==colnames(zscore.matrix)[j]]#

#      AnnotationB.centroid=apply(PCA.space[AnnotationB,],2,median)
      
 #     true.distance=c()
      #Bootstraping true distance
   #   bootstrapping.distance=c()
  #    bootstrapping.sampling.distance=c()
     # for (k in 1:roundsofrandomization) {
        # Booststrap
      #  sub_sampling_size=runif(1,min=0.5,max = 0.75)
      #  if (sub_sampling_size*length(AnnotationA)<=2) {
        #  sub_sampling_1=PCA.space[AnnotationA,]

        #  temp=PCA.space[sample(AnnotationA,size = 2,replace = T),]
        #  rownames(temp)=paste0(rownames(temp),'bootstrap')
         # sub_sampling_1=rbind(sub_sampling_1,temp)
          
         # sub_sampling_2=PCA.space[AnnotationB,]
        #  temp=PCA.space[sample(AnnotationB,size = 2,replace = T),]
         # rownames(temp)=paste0(rownames(temp),'bootstrap')
         # sub_sampling_2=rbind(sub_sampling_2,temp)

         # sub_sampling_1=apply(sub_sampling_1,2,median)
        #  sub_sampling_2=apply(sub_sampling_2,2,median)
       # }
      #  else {
          
        #  sub_sampling_1=sample(AnnotationA,size = sub_sampling_size*length(AnnotationA) ,replace = T)
       #   sub_sampling_2=sample(AnnotationB,size = sub_sampling_size*length(AnnotationB) ,replace = T)
       #   sub_sampling_1=apply(PCA.space[sub_sampling_1,],2,median)
       #   sub_sampling_2=apply(PCA.space[sub_sampling_2,],2,median)
      #  }

       # sub_sampling_distance=euclidean(sub_sampling_1,sub_sampling_2)
      #  true.distance=c(true.distance,sub_sampling_distance)
        
      #  sampling_pool=rownames(PCA.space)
      #  sampling_pool.A=sampling_pool[sampling_pool%in%rownames(meta.data)[meta.data$object==Name_of_object_A]]
      #  sampling_pool.B=sampling_pool[sampling_pool%in%rownames(meta.data)[meta.data$object==Name_of_object_B]]
        
      #  sampling.subset.A=sample(sampling_pool.A,size = (length(AnnotationA)+length(AnnotationB))/2,replace = T)
      #  sampling.subset.B=sample(sampling_pool.B,size = (length(AnnotationA)+length(AnnotationB))/2,replace = T)
      #  sampling.subset.A.centroid=apply(PCA.space[sampling.subset.A,],2,median)
        
      #  sampling.subset.B.centroid=apply(PCA.space[sampling.subset.B,],2,median)

      #  bootstrapping.sampling.distance=c(bootstrapping.sampling.distance,(euclidean(sampling.subset.B.centroid,AnnotationA.centroid)+euclidean(sampling.subset.A.centroid,AnnotationB.centroid))/2)
      #}
      
      #Comparison.results=data.frame(Distance=bootstrapping.sampling.distance,
      #                              group=rep('Randomized',length(bootstrapping.sampling.distance)))
      #Comparison.results=rbind(Comparison.results,
      #                         data.frame(Distance=true.distance,
      #                                    group=rep(paste(rownames(zscore.matrix)[i],colnames(zscore.matrix)[j]),length(true.distance))))
      #Comparison.results$X=rownames(zscore.matrix)[i]
      #Comparison.results$Y=colnames(zscore.matrix)[j]
      #Detail_list[[rownames(zscore.matrix)[i]]]$Aggregated.randomization.results=rbind(Detail_list[[rownames(zscore.matrix)[i]]]$Aggregated.randomization.results,
      #                                                                                 Comparison.results)
      #wilcox.results=wilcox.test(log10(true.distance+1),log10(bootstrapping.sampling.distance+1),alternative = 'less')
      #
      #Detail_list[[rownames(zscore.matrix)[i]]][[colnames(zscore.matrix)[j]]]=list(Comparison.results=Comparison.results,
      #                                                                             wilcox.results=wilcox.results)

      #zscore.matrix[i,j]=median(true.distance)-median(bootstrapping.sampling.distance)
      #pval.matrix[i,j]=wilcox.results$p.value
      
    #}
    #Detail_list[[rownames(zscore.matrix)[i]]]$Aggregated.randomization.results=Detail_list[[rownames(zscore.matrix)[i]]]$Aggregated.randomization.results[-1,]
    #Aggregated.randomization.results=Detail_list[[rownames(zscore.matrix)[i]]]$Aggregated.randomization.results
    
    #Aggregated.Randomization.Results.plot=ggplot(Aggregated.randomization.results[Aggregated.randomization.results$group!='Randomized',],
     #                                                    aes(x=Distance,fill=Y))+geom_density(alpha=0.75)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
    #Detail_list[[rownames(zscore.matrix)[i]]]$Aggregated.Randomization.Results.plot=Aggregated.Randomization.Results.plot
    
    
    
    
    #print(i/nrow(zscore.matrix))
  #}
  

  Overall.Aggregated.Randomization.Results=data.frame(Distance=0,group='',X='',Y='')
  for (comparison in names(Detail_list)) {
    Overall.Aggregated.Randomization.Results=rbind(Overall.Aggregated.Randomization.Results,
                                                   Detail_list[[comparison]]$Aggregated.randomization.results)
  }
  Overall.Aggregated.Randomization.Results=unique(Overall.Aggregated.Randomization.Results)

  
  
  Overall.Aggregated.Randomization.Results=Overall.Aggregated.Randomization.Results[-1,]
  collapsed.results.z=melt(zscore.matrix)
  collapsed.results.z=data.frame(collapsed.results.z)
  colnames(collapsed.results.z)=c('X','Y','ZSCORE')
  collapsed.results.p=melt(pval.matrix)
  collapsed.results.p=data.frame(collapsed.results.p)
  colnames(collapsed.results.p)=c('X','Y','Pval')
  collapsed.results=merge(collapsed.results.z,collapsed.results.p,by = c('X','Y'))
  for (i in 1:nrow(collapsed.results)) {
    if (collapsed.results$Pval[i]==0) {
      collapsed.results$Pval[i]=.Machine$double.xmin
    }
  }
  collapsed.results$padj=p.adjust(collapsed.results$Pval,method = 'BH')
  collapsed.results$logpval=-log10(collapsed.results$padj)
  collapsed.results$similarityscore=((-1)*collapsed.results$ZSCORE)
  collapsed.results$similarityscore=log10(abs(collapsed.results$similarityscore)+1)/(collapsed.results$similarityscore/abs(collapsed.results$similarityscore))
  collapsed.results$significany='Insignificant'
  for (i in 1:nrow(collapsed.results)) {
    if (collapsed.results$padj[i]<=p.val.threshold) {
      cross_validation=wilcox.test(log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[as.character(Overall.Aggregated.Randomization.Results$group)==paste(as.character(collapsed.results$X[i]),
                                                                                                                                     as.character(collapsed.results$Y[i]),sep = ' ')])+1),
                                   log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$X==collapsed.results$X[i]&Overall.Aggregated.Randomization.Results$group!='Randomized'])+1),alternative = 'less')
      if (cross_validation$p.value<=p.val.threshold){
        collapsed.results$significany[i]='Significant'
      }
    }
  }
  collapsed.results$prediction.similarity='Insignificant'
  for (i in 1:nrow(collapsed.results)) {
    if (collapsed.results$significany[i]=='Significant') {
      X=as.character(collapsed.results$X[i])
      #listofY=unique(as.character(collapsed.results$Y[collapsed.results$X==X]))
      #average=data.frame(Y=listofY,average=rep(0,length(listofY)))
      #for (j in 1:nrow(average)) {
      #  average$average[j]=mean(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$X==X&
      #                                                                            Overall.Aggregated.Randomization.Results$Y==average$Y[j]&
       #                                                                           Overall.Aggregated.Randomization.Results$group!='Randomized'])
      #}
      if (collapsed.results$similarityscore[i]==max(collapsed.results$similarityscore[collapsed.results$X==X])){
        candidates=collapsed.results$Y[i]
        
        cross_validation=wilcox.test(log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$X==X&Overall.Aggregated.Randomization.Results$Y==candidates&Overall.Aggregated.Randomization.Results$group!='Randomized'])+1),
                                     log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$Y==candidates&Overall.Aggregated.Randomization.Results$group!='Randomized'])+1),alternative = 'less')
        global_validation=wilcox.test(log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$X==X&Overall.Aggregated.Randomization.Results$Y==candidates&Overall.Aggregated.Randomization.Results$group!='Randomized'])+1),
                                     log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$group!='Randomized'])+1),alternative = 'less')
        
        if (cross_validation$p.value<=p.val.threshold&global_validation$p.value<=p.val.threshold) {
          collapsed.results$prediction.similarity[i]='Significant'
          
        }
        }
    }
  }
  Prediction.plot=ggplot(collapsed.results,aes(x=X,y=Y))+
    geom_point(aes(size=logpval,
                   colour=prediction.similarity))+theme(axis.text.x = element_text(angle=45,hjust=1),
                                                            axis.title.x = element_blank(),
                                                            axis.title.y = element_blank(),
                                                            panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(),
                                                            panel.background = element_blank(),
                                                            axis.line = element_line(colour = "black"))
  
  
  
  
  significancy.plot=ggplot(collapsed.results,aes(x=X,y=Y))+
    geom_point(aes(size=logpval,
                   colour=significany))+theme(axis.text.x = element_text(angle=45,hjust=1),
                                              axis.title.x = element_blank(),
                                              axis.title.y = element_blank(),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              panel.background = element_blank(),
                                              axis.line = element_line(colour = "black"))
  Overall.Aggregated.Randomization.Results.plot.1=ggplot(Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$group!='Randomized',],
                                                       aes(x=Distance,fill=Y))+geom_density(alpha=0.75)+facet_wrap(facets = ~X,ncol = prediction_plot_ncol)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  Overall.Aggregated.Randomization.Results.plot.2=ggplot(Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$group!='Randomized',],
                                                         aes(x=Distance,fill=X))+geom_density(alpha=0.75)+facet_wrap(facets = ~Y,ncol = prediction_plot_ncol)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  results=list(zscore.matrix=zscore.matrix,
               pval.matrix=pval.matrix,
               collapsed.results=collapsed.results,
               Similarity.significancy.plot=significancy.plot,
               Detail_list=Detail_list,
               Overall.Aggregated.Randomization.Results=Overall.Aggregated.Randomization.Results,
               Overall.Aggregated.Randomization.Results.plot.1=Overall.Aggregated.Randomization.Results.plot.1,
               Overall.Aggregated.Randomization.Results.plot.2=Overall.Aggregated.Randomization.Results.plot.2,
               Prediction.plot=Prediction.plot)
  return(results)
}





Aggregated.distribution.plotter=function(Overall.Aggregated.Randomization.Results,ncol=2) {
  library(ggplot2)
  Overall.Aggregated.Randomization.Results.plot=ggplot(Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$group!='Randomized',],
                                                       aes(x=Distance,fill=Y))+geom_density()+facet_wrap(facets = ~X,ncol = ncol)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  return(Overall.Aggregated.Randomization.Results.plot)
}





Subsetter=function(object,annotation_name,percentage,exclude=NULL) {
  #return subset cell names
  results=c()
  if (is.null(exclude)){
  for (i in unique(object@meta.data[[annotation_name]])) {
    temp=sample(rownames(object@meta.data)[object@meta.data[[annotation_name]]==i],size = length(rownames(object@meta.data)[object@meta.data[[annotation_name]]==i])*percentage,replace = F)
    results=c(results,temp)
  }
  }
  else {
    for (i in unique(object@meta.data[[annotation_name]])[!unique(object@meta.data[[annotation_name]])%in%exclude]) {
      temp=sample(rownames(object@meta.data)[object@meta.data[[annotation_name]]==i],size = length(rownames(object@meta.data)[object@meta.data[[annotation_name]]==i])*percentage,replace = F)
      results=c(results,temp)
    } 
  }
  return(results)
}