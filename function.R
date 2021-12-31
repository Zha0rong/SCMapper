SCMapper=function(Object_A,Object_B,
         annotation_name_A='seurat_clusters',annotation_name_B='seurat_clusters',
         Name_of_object_A,Name_of_object_B,
         Batch_A=NULL,Batch_B=NULL,
         numbersofPC=30,roundsofrandomization=1000,
         Normalization_method='LogNormalize',nFeature=2000,harmonize=F,
         p.val.threshold=0.05,prediction_plot_ncol=2,number_of_core=4,mt.feature.pattern=NULL,correct_by_expression_Matrix=T) {
  require(cluster)
  require(reshape)
  require(ggplot2)
  require(Seurat)
  require(harmony)
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(lsa)
  require(scater)
  require(muscat)
  source('utils.R')
  

  Object_A@meta.data[['scmapper']]=paste(Name_of_object_A,
                                         '_',Object_A@meta.data[[annotation_name_A]],sep = '')
  levels(Object_A@meta.data[['scmapper']])=paste(Name_of_object_A,
                                                 '_',levels(Object_A@meta.data[[annotation_name_A]]),sep = '')
  Object_A@meta.data[['object']]=Name_of_object_A
  Object_B@meta.data[['scmapper']]=paste(Name_of_object_B,
                                         '_',Object_B@meta.data[[annotation_name_B]],sep = '')
  levels(Object_B@meta.data[['scmapper']])=paste(Name_of_object_B,
                                                 '_',levels(Object_B@meta.data[[annotation_name_B]]),sep = '')
  Object_B@meta.data[['object']]=Name_of_object_B
  
  Merged=merge(Object_A,Object_B)
  filtering <- SingleCellExperiment(
    assays = list(counts = Merged@assays$RNA@counts), 
    colData = Merged@meta.data
  )
  
  filtering <- filtering[rowSums(counts(filtering) > 1) >= 1, ]
  Merged=CreateSeuratObject(counts(filtering),meta.data = data.frame(filtering@colData))
  rm(filtering)

  
  Merged$percent.mt=ifelse(is.null(mt.feature.pattern),yes = 0,no=PercentageFeatureSet(Merged,pattern = mt.feature.pattern))
  if (all(Merged[['percent.mt']]==0)) {
    Merged[['percent.mt']]=NULL
  }
  DefaultAssay(Merged)='RNA'
  if (Normalization_method=='SCT'){
    ncells=ifelse(ncol(Merged)*0.5<5000,yes = 5000,no=5000)
    
    if ('percent.mt'%in%names(Merged@meta.data)) {
      Merged=SCTransform(Merged,conserve.memory = T,vars.to.regress = 'percent.mt',
                         verbose = F,variable.features.n = nFeature,ncells = ncells)
      Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)

    }
    else {
      Merged=SCTransform(Merged,conserve.memory = T,variable.features.n = nFeature,ncells = ncells,verbose = F)
      Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)

    }
    if (correct_by_expression_Matrix) {
      genes.test=Merged[['SCT']]@var.features
      Expression.Correlation.Matrix=Expression_correlation(Merged=Merged,annotation_Name = 'scmapper',object_Name = 'object',Name_of_object_A = Name_of_object_A,
                                                           Name_of_object_B=Name_of_object_B,assay.use = 'SCT')
    }
    else {
      genes.test=0
      Expression.Correlation.Matrix=0
    }
    print('Finished SCT Normalization')
  }
  if (Normalization_method=='LogNormalize'){
    Merged=NormalizeData(Merged,verbose = F)
    Merged=FindVariableFeatures(Merged,verbose = F,nfeatures = nFeature)
    if (correct_by_expression_Matrix) {
    genes.test=Merged[['RNA']]@var.features
    Expression.Correlation.Matrix=Expression_correlation(Merged=Merged,annotation_Name = 'scmapper',object_Name = 'object',Name_of_object_A = Name_of_object_A,
                                                         Name_of_object_B=Name_of_object_B,assay.use = 'RNA')
    }
    else {
      genes.test=0
      Expression.Correlation.Matrix=0
    }
    
    if ('percent.mt'%in%names(Merged@meta.data)) {
      Merged=ScaleData(Merged,verbose = F,vars.to.regress = 'percent.mt')
      
    }
    else {
      Merged=ScaleData(Merged,verbose = F)
    }
    Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)

    print('Finished Log Normalization')
  }
  if (Normalization_method=='None'){
    Merged=FindVariableFeatures(Merged,verbose = F,nfeatures = nFeature)
    Merged=ScaleData(Merged,verbose = F)
    Merged=RunPCA(Merged,npcs = numbersofPC,verbose = F)
    if (correct_by_expression_Matrix) {
    genes.test=Merged[['RNA']]@var.features
    Expression.Correlation.Matrix=Expression_correlation(Merged=Merged,annotation_Name = 'scmapper',object_Name = 'object',Name_of_object_A = Name_of_object_A,
                                                         Name_of_object_B=Name_of_object_B,assay.use = 'RNA')
    }
    else {
      genes.test=0
      Expression.Correlation.Matrix=0
    }
    
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
    Detail_list[[i]]$ttest.results=list()
    
  }
  ####If harmonized####
  if (harmonize) {

    PCA.space=harmonize_within_sample(harmonize=harmonize,
                                      Normalization_method=Normalization_method,Merged,
                                      numbersofPC=numbersofPC,
                                      Batch_A=Batch_A,
                                      Batch_B=Batch_B,
                                      Name_of_object_A=Name_of_object_A,
                                      Name_of_object_B=Name_of_object_B)
    harmonized.space=PCA.space$harmonized.space
    original.space=PCA.space$original.space
    
    meta.data=Merged@meta.data
    DefaultAssay(Merged)='RNA'
    Merged=NormalizeData(Merged,verbose = F)
    if (correct_by_expression_Matrix){
    #Matrix=as.matrix(Merged@assays$RNA@data[genes.test,,,])
      Matrix=0
      
    }
    else {
      Matrix=0
    }
    
    rm(Merged)
    cl <- makeCluster(number_of_core)
    registerDoSNOW(cl)
    
    Test_Candidates=paste(rep(rownames(zscore.matrix), each = length(colnames(zscore.matrix))), colnames(zscore.matrix), sep = "Object_42")
    iterations <- length(Test_Candidates)
    
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    
    opts <- list(progress = progress)

    x=foreach(i=Test_Candidates,.packages=c('ggplot2','lsa','fpc','cluster'), .options.snow = opts) %dopar% {
      source('utils.R')
      
      Results=sampling_and_distance_test(i,mapping.table=mapping.table,harmonized.space=harmonized.space,original.space=original.space,roundsofrandomization=roundsofrandomization
                                         ,meta.data=meta.data,Name_of_object_A=Name_of_object_A,Name_of_object_B=Name_of_object_B,
                                         Expression.Correlation.Matrix=Expression.Correlation.Matrix,Matrix=Matrix,correct_by_expression_Matrix=correct_by_expression_Matrix)
      return(Results)
    }
  }
  ####If not harmonized####
  else {
    PCA.space=harmonize_within_sample(harmonize=harmonize,
                            Normalization_method=Normalization_method,Merged,
                            numbersofPC=numbersofPC,
                            Batch_A=Batch_A,
                            Batch_B=Batch_B,
                            Name_of_object_A=Name_of_object_A,
                            Name_of_object_B=Name_of_object_B) 

  meta.data=Merged@meta.data
  
  
  DefaultAssay(Merged)='RNA'
  Merged=NormalizeData(Merged,verbose = F)
  if (correct_by_expression_Matrix){
    #Matrix=as.matrix(Merged@assays$RNA@data[genes.test,,,])
    Matrix=0
    
  }
  else {
    Matrix=0
  }  

  rm(Merged)
  cl <- makeCluster(number_of_core)
  registerDoSNOW(cl)
  Test_Candidates=paste(rep(rownames(zscore.matrix), each = length(colnames(zscore.matrix))), colnames(zscore.matrix), sep = "Object_42")
  iterations <- length(Test_Candidates)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=Test_Candidates,.packages=c('ggplot2','lsa','fpc','cluster'), .options.snow = opts) %dopar% {
    source('utils.R')
    
    Results=sampling_and_distance_test(i,mapping.table=mapping.table,PCA.space=PCA.space,roundsofrandomization=roundsofrandomization,meta.data=meta.data
                                       ,Name_of_object_A=Name_of_object_A,Name_of_object_B=Name_of_object_B,
                                       Expression.Correlation.Matrix=Expression.Correlation.Matrix,Matrix=Matrix,correct_by_expression_Matrix=correct_by_expression_Matrix)
    return(Results)
  }
  }
  close(pb)
  stopCluster(cl)
  
  for (i in 1:length(x)) {
    temp=x[[i]]

    zscore.matrix[temp$a,temp$b]=temp$zscore
    pval.matrix[temp$a,temp$b]=temp$wilcox.results$p.value
    Detail_list[[temp$a]]$Aggregated.randomization.results=rbind(Detail_list[[temp$a]]$Aggregated.randomization.results,
                                                                    temp$Aggregated.randomization.results)
  }
  for (i in rownames(zscore.matrix)) {
    Detail_list[[i]]$Aggregated.randomization.results=Detail_list[[i]]$Aggregated.randomization.results[-1,]
    Aggregated.randomization.results=Detail_list[[i]]$Aggregated.randomization.results
    Aggregated.Results.plot=ggplot(Aggregated.randomization.results[Aggregated.randomization.results$group!='Randomized',],
                                   aes(x=Distance,fill=Y))+geom_density(alpha=0.75)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
    
    Detail_list[[i]]$Aggregated.Results.plot=Aggregated.Results.plot
    Aggregated.Randomization.Results.plot=ggplot(Aggregated.randomization.results,
                                                 aes(x=Distance,fill=group))+geom_density(alpha=0.75)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
    Detail_list[[i]]$Aggregated.Randomization.Results.plot=Aggregated.Randomization.Results.plot
    
  }
  rm(x)
  gc()
  
  
  ####Aggregated Results from randomization####
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
                                   (log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$Y==collapsed.results$Y[i]&Overall.Aggregated.Randomization.Results$group!='Randomized'])+1)
                                          ),alternative = 'less')
      if (cross_validation$p.value<=p.val.threshold){
        collapsed.results$significany[i]='Significant'
      }
    }
  }
  collapsed.results$prediction.similarity='Unmapped'
  collapsed.results$pval_list=''
  for (i in 1:nrow(collapsed.results)) {
    if (collapsed.results$significany[i]=='Significant') {
      X=as.character(collapsed.results$X[i])
      candidates=as.character(collapsed.results$Y[i])
      
      temp_X=Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$X==X
                                                      &Overall.Aggregated.Randomization.Results$Y==candidates
                                                      &Overall.Aggregated.Randomization.Results$group!='Randomized',]
      temp_all=Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$X==X
                                                      &Overall.Aggregated.Randomization.Results$group!='Randomized',]
      
      temp_Y=Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$Y==candidates
                                                      &Overall.Aggregated.Randomization.Results$group!='Randomized',]
      pval_list=c()
      median_list=c()
      for (idents in unique(temp_all$Y)[unique(temp_all$Y)!=candidates]) {
        temp=temp_all[temp_all$Y==idents,]
        test=wilcox.test(log10(temp_X$Distance+1),log10(temp$Distance+1),alternative = 'less')
        pval_list=c(pval_list,test$p.value)
      }
      for (idents in unique(temp_all$Y)) {
        temp=temp_all[temp_all$Y==idents,]
        median_list=c(median_list,mean(temp$Distance))
        
      }
      
      
      
      
      collapsed.results$pval_list[i]=paste(pval_list,collapse = ' ')
      
      if (mean(temp_X$Distance)==min(median_list)&all(pval_list<=p.val.threshold)){
        Test_reference=(log10(( Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$Y==candidates&Overall.Aggregated.Randomization.Results$X!=X&
                                                                                    Overall.Aggregated.Randomization.Results$group!='Randomized'  ])+1))
        Test_reference=sample(Test_reference,size = roundsofrandomization)
        cross_validation=wilcox.test(log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$X==X&Overall.Aggregated.Randomization.Results$Y==candidates&Overall.Aggregated.Randomization.Results$group!='Randomized'])+1),
                                Test_reference
                                      ,alternative = 'less')
        Test_reference=log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$group!='Randomized'])+1)
        Test_reference=sample(Test_reference,size = roundsofrandomization)
        
       global_validation=wilcox.test(log10(as.numeric(Overall.Aggregated.Randomization.Results$Distance[Overall.Aggregated.Randomization.Results$X==X&
                                                                                                      Overall.Aggregated.Randomization.Results$Y==candidates&Overall.Aggregated.Randomization.Results$group!='Randomized'])+1),
                                     Test_reference
                                            ,alternative = 'less')
        
        if (cross_validation$p.value<=p.val.threshold&global_validation$p.value<=p.val.threshold) {
          collapsed.results$prediction.similarity[i]='Mapped'
          
        }
      }
    }
  }

  ####Generate figures####
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
                                                       aes(x=log10(Distance+1),fill=Y))+geom_density(alpha=0.75)+facet_wrap(facets = ~X,ncol = prediction_plot_ncol)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  Overall.Aggregated.Randomization.Results.plot.2=ggplot(Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$group!='Randomized',],
                                                         aes(x=log10(Distance+1),fill=X))+geom_density(alpha=0.75)+facet_wrap(facets = ~Y,ncol = prediction_plot_ncol)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  ####Return Results####
  results=list(zscore.matrix=zscore.matrix,
               pval.matrix=pval.matrix,
               collapsed.results=collapsed.results,
               Similarity.significancy.plot=significancy.plot,
               Detail_list=Detail_list,
               Overall.Aggregated.Randomization.Results=Overall.Aggregated.Randomization.Results,
               Overall.Aggregated.Randomization.Results.plot.1=Overall.Aggregated.Randomization.Results.plot.1,
               Overall.Aggregated.Randomization.Results.plot.2=Overall.Aggregated.Randomization.Results.plot.2,
               Prediction.plot=Prediction.plot,Expression.Correlation.Matrix=Expression.Correlation.Matrix,
               PCA.space=PCA.space,
               mapping.table=mapping.table)
  return(results)
}

Single.comparison.plotter=function(Overall.Aggregated.Randomization.Results,X=NULL,Y=NULL) {
  if (!is.null(X)&!is.null(Y)){
    table=Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$X==X&Overall.Aggregated.Randomization.Results$Y==Y,]
    Plot=ggplot(table,
                                                           aes(x=Distance,fill=group))+geom_density(alpha=0.75)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  }
  else {
    return(list())
  }
  return(Plot)
  
  
  
}



Aggregated.distribution.plotter=function(Overall.Aggregated.Randomization.Results,ncol=2,exclude_randomized=F) {
  library(ggplot2)
  if (exclude_randomized) {
    Overall.Aggregated.Randomization.Results.plot=ggplot(Overall.Aggregated.Randomization.Results[Overall.Aggregated.Randomization.Results$group!='Randomized',],
                                                       aes(x=Distance,fill=Y))+geom_density()+facet_wrap(facets = ~X,ncol = ncol)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  }
  else {
    Overall.Aggregated.Randomization.Results.plot=ggplot(Overall.Aggregated.Randomization.Results,
                                                         aes(x=Distance,fill=group))+geom_density()+facet_wrap(facets = ~X,ncol = ncol)+theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
    
  }
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






