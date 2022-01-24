sampling_and_distance_test <- function(TestCandidates,mapping.table,harmonized.space=NULL,
                                       PCA.space=NULL,original.space=NULL,roundsofrandomization,meta.data,Name_of_object_A,Name_of_object_B
                                       ,Expression.Correlation.Matrix,Matrix,correct_by_expression_Matrix) {
  a=as.character(unlist(strsplit(TestCandidates,split = "Object_42")))[1]
  b=as.character(unlist(strsplit(TestCandidates,split = "Object_42")))[2]
  require(fpc)
  require(cluster)

  euclidean <- function(x, y) sqrt(sum((x - y)^2))

  if (!is.null(original.space)) {
    PCA.space=harmonized.space
  }
  AnnotationA=mapping.table$cell[mapping.table$cluster==a]
  AnnotationA.centroid=apply(PCA.space[AnnotationA,],2,median)
  Aggregated.randomization.results=data.frame(Distance=0,group='',X='',Y='')
  AnnotationB=mapping.table$cell[mapping.table$cluster==b]
  AnnotationB.centroid=apply(PCA.space[AnnotationB,],2,median)
  true.distance=c()
  bootstrapping.sampling.distance=c()
  for (k in 1:roundsofrandomization) {
    sub_sampling_size=runif(1,0.25,0.75)
    if (sub_sampling_size*length(AnnotationA)<=2) {
      sub_sampling_1=PCA.space[AnnotationA,]
      
      temp=PCA.space[sample(AnnotationA,size = 2,replace = T),]
      rownames(temp)=paste0(rownames(temp),'bootstrap')
      cells_1=rownames(temp)
      sub_sampling_1=rbind(sub_sampling_1,temp)
      sub_sampling_2=PCA.space[AnnotationB,]
      temp=PCA.space[sample(AnnotationB,size = 2,replace = T),]
      rownames(temp)=paste0(rownames(temp),'bootstrap')
      cells_2=rownames(temp)
      
      sub_sampling_2=rbind(sub_sampling_2,temp)
      sub_sampling_1=apply(sub_sampling_1,2,median)
      sub_sampling_2=apply(sub_sampling_2,2,median)
    }
    else {
      sub_sampling_1=sample(AnnotationA,size = sub_sampling_size*length(AnnotationA) ,replace = F)
      cells_1=sub_sampling_1
      sub_sampling_2=sample(AnnotationB,size = sub_sampling_size*length(AnnotationB) ,replace = F)
      cells_2=sub_sampling_2
      
      sub_sampling_1=apply(PCA.space[sub_sampling_1,],2,median)
      sub_sampling_2=apply(PCA.space[sub_sampling_2,],2,median)
    }
    sub_sampling_distance=euclidean(sub_sampling_1,sub_sampling_2)
    if (correct_by_expression_Matrix) {
      true.correlation=Expression.Correlation.Matrix[a,b]
    }
    else {
      true.correlation=1
    }
    
    true.distance=c(true.distance,sub_sampling_distance/true.correlation)
    sampling_space_A=rownames(meta.data)[meta.data$object==Name_of_object_A]
    sampling_space_B=rownames(meta.data)[meta.data$object==Name_of_object_B]
    sampling_A_centroid_cells=sample(sampling_space_A,size = sub_sampling_size*length(AnnotationA))
    sampling_B_centroid_cells=sample(sampling_space_B,size = sub_sampling_size*length(AnnotationB))
    sampling_A_centroid=apply(PCA.space[sampling_A_centroid_cells,],2,median)
    sampling_B_centroid=apply(PCA.space[sampling_B_centroid_cells,],2,median)
    if (correct_by_expression_Matrix) {
      
    #sub_sampling_1.to.sampling_B.correlation=Expression_correlation_Randomized(Matrix,cells_1,sampling_B_centroid_cells)
    #sub_sampling_2.to.sampling_A.correlation=Expression_correlation_Randomized(Matrix,cells_2,sampling_A_centroid_cells)
      sub_sampling_1.to.sampling_B.correlation=sample(as.numeric(Expression.Correlation.Matrix[a,]),size = 1)
      sub_sampling_2.to.sampling_A.correlation=sample(as.numeric(Expression.Correlation.Matrix[,b]),size = 1)
      #correlation=quantile(as.numeric(Expression.Correlation.Matrix))[4]
    }
    else {
      sub_sampling_1.to.sampling_B.correlation=1
      sub_sampling_2.to.sampling_A.correlation=1
    }
    sampling.distance=(((euclidean(sub_sampling_1,sampling_B_centroid)/sub_sampling_1.to.sampling_B.correlation)
                       +(euclidean(sub_sampling_2,sampling_A_centroid)/sub_sampling_2.to.sampling_A.correlation))/2)
    bootstrapping.sampling.distance=c(bootstrapping.sampling.distance,sampling.distance)
    }
    wilcox.results=wilcox.test(log10(true.distance+1),log10(bootstrapping.sampling.distance+1),alternative = 'less', paired = T)
    Comparison.results=data.frame(Distance=bootstrapping.sampling.distance,group=rep('Randomized',length(bootstrapping.sampling.distance)))
    #/Expression.Correlation.Matrix[a,b]
    Comparison.results=rbind(Comparison.results,data.frame(Distance=true.distance,group=rep(paste(a,b),length(true.distance))))
    Comparison.results$X=a
    Comparison.results$Y=b
    Aggregated.randomization.results=rbind(Aggregated.randomization.results,Comparison.results)
    zscore=median(true.distance)-median(bootstrapping.sampling.distance)
    Aggregated.randomization.results=Aggregated.randomization.results[-1,]
    Results=list(a=a,
                 b=b,
                 zscore=zscore,
                 wilcox.results=wilcox.results,
                 Aggregated.randomization.results=Aggregated.randomization.results)
  
  
  

  return(Results)
}





matrix_distance_test <- function(TestCandidates,mapping.table,harmonized.space=NULL,PCA.space=NULL,original.space=NULL,roundsofrandomization,meta.data,Name_of_object_A,Name_of_object_B) {
  a=as.character(unlist(strsplit(TestCandidates,split = "Object_42")))[1]
  b=as.character(unlist(strsplit(TestCandidates,split = "Object_42")))[2]
  require(fpc)
  require(cluster)
  euclidean <- function(x, y) sqrt(sum(((x) - (y))^2))
  
  if (!is.null(original.space)) {
    samplesize=0.5
    AnnotationA=mapping.table$cell[mapping.table$cluster==a]
    
    AnnotationA.matrix=harmonized.space[sample(AnnotationA,size = length(AnnotationA)*samplesize),]
    
    Aggregated.randomization.results=data.frame(Distance=0,group='',X='',Y='')
    AnnotationB=mapping.table$cell[mapping.table$cluster==b]
    
    AnnotationB.matrix=harmonized.space[sample(AnnotationB,size = length(AnnotationB)*samplesize),]
    sampling_pool=rownames(harmonized.space)
    sampling_pool.B=sampling_pool[sampling_pool%in%rownames(meta.data)[meta.data$object==Name_of_object_B]]
    sampling.subset.B=sample(sampling_pool.B,size = length(AnnotationB)*samplesize,replace = F)
    sampling.subset.B=harmonized.space[sampling.subset.B,]
    true.distance=c()
    bootstrapping.sampling.distance=c()
    for (i in 1:nrow(AnnotationA.matrix)) {
      for (j in 1:nrow(AnnotationB.matrix)) {
        true.distance=c(true.distance,euclidean(AnnotationA.matrix[i,],AnnotationB.matrix[j,]))
        bootstrapping.sampling.distance=c(bootstrapping.sampling.distance,euclidean(AnnotationA.matrix[i,],sampling.subset.B[j,]))
        
      }
    }
    
  }
  else {
    samplesize=0.5
    
    AnnotationA=mapping.table$cell[mapping.table$cluster==a]
    
    AnnotationA.matrix=PCA.space[sample(AnnotationA,size = length(AnnotationA)*samplesize),]
    
    Aggregated.randomization.results=data.frame(Distance=0,group='',X='',Y='')
    AnnotationB=mapping.table$cell[mapping.table$cluster==b]
    
    AnnotationB.matrix=PCA.space[sample(AnnotationB,size = length(AnnotationB)*samplesize),]
    sampling_pool=rownames(PCA.space)
    sampling_pool.B=sampling_pool[sampling_pool%in%rownames(meta.data)[meta.data$object==Name_of_object_B]]
    sampling.subset.B=sample(sampling_pool.B,size = length(AnnotationB)*samplesize,replace = F)
    sampling.subset.B=PCA.space[sampling.subset.B,]
    true.distance=c()
    bootstrapping.sampling.distance=c()
    for (i in 1:nrow(AnnotationA.matrix)) {
      for (j in 1:nrow(AnnotationB.matrix)) {
        true.distance=c(true.distance,euclidean(AnnotationA.matrix[i,],AnnotationB.matrix[j,]))
        bootstrapping.sampling.distance=c(bootstrapping.sampling.distance,euclidean(AnnotationA.matrix[i,],sampling.subset.B[j,]))
        
      }
    }
    
  }
  Comparison.results=data.frame(Distance=bootstrapping.sampling.distance,group=rep('Randomized',length(bootstrapping.sampling.distance)))
  Comparison.results=rbind(Comparison.results,data.frame(Distance=true.distance,group=rep(paste(a,b),length(true.distance))))
  Comparison.results$X=a
  Comparison.results$Y=b
  Aggregated.randomization.results=rbind(Aggregated.randomization.results,
                                         Comparison.results)
  
  wilcox.results=wilcox.test(log10(true.distance+1),(log10(bootstrapping.sampling.distance+1)),alternative = 'less', paired = T,correct = T)
  
  
  
  
  
  zscore=median(true.distance)-median(bootstrapping.sampling.distance)
  Aggregated.randomization.results=Aggregated.randomization.results[-1,]
  Results=list(a,b,zscore,wilcox.results,Aggregated.randomization.results)
  return(Results)
}





harmonize_within_sample <- function(harmonize,Normalization_method,Merged,numbersofPC,Batch_A,Batch_B,Name_of_object_A,Name_of_object_B) {
  require(harmony)
  if (harmonize) {
  if (Normalization_method=='SCT') {
    Merged=RunHarmony(Merged,group.by.vars = 'object',max.iter.harmony = 1000,
                      max.iter.cluster = 1000,verbose = F,assay.use='SCT')
    harmonized.space=Merged@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
    original.space=Merged@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
    
  }
  else {
    Merged=RunHarmony(Merged,group.by.vars = 'object',max.iter.harmony = 1000,
                      max.iter.cluster = 1000,verbose = F)
    harmonized.space=Merged@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
    original.space=Merged@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
    
  }
  
  if (!is.null(Batch_A)|!is.null(Batch_B)) {
    if (!is.null(Batch_A)) {
      PCA.A=subset(Merged,subset=object==Name_of_object_A)
      if (Normalization_method=='SCT') {
        PCA.A=RunHarmony(PCA.A,group.by.vars=Batch_A,assay.use='SCT',max.iter.harmony = 1000,max.iter.cluster = 1000,reduction='harmony',reduction.save='harmonized')
        Original.A=subset(Merged,subset=object==Name_of_object_A)
        
      }
      else {
        PCA.A=RunHarmony(PCA.A,group.by.vars=Batch_A,max.iter.harmony = 1000,max.iter.cluster = 1000,reduction='harmony',reduction.save='harmonized')
        Original.A=subset(Merged,subset=object==Name_of_object_A)
        
      }
    }
    if (!is.null(Batch_B)) {
      PCA.B=subset(Merged,subset=object==Name_of_object_B)
      if (Normalization_method=='SCT') {
        PCA.B=RunHarmony(PCA.B,group.by.vars=Batch_B,assay.use='SCT',max.iter.harmony = 1000,max.iter.cluster = 1000,reduction='harmony',reduction.save='harmonized')
        Original.B=subset(Merged,subset=object==Name_of_object_B)
        
      }
      else {
        PCA.B=RunHarmony(PCA.B,group.by.vars=Batch_B,max.iter.harmony = 1000,max.iter.cluster = 1000,reduction='harmony',reduction.save='harmonized')
        Original.B=subset(Merged,subset=object==Name_of_object_B)
        
        
      }
    }
    
    if (!is.null(Batch_A)&is.null(Batch_B)) {
      PCA.A=PCA.A@reductions$harmonized@cell.embeddings[,seq(1,numbersofPC)]
      Part_B=subset(Merged,subset=object==Name_of_object_B)@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
      colnames(PCA.A)=colnames(Part_B)
      harmonized.space=rbind(PCA.A,Part_B)
      
      Original.A=Original.A@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
      Part_B=subset(Merged,subset=object==Name_of_object_B)@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
      colnames(Original.A)=colnames(Part_B)
      original.space=rbind(Original.A,Part_B)
      
      
    }
    else if (is.null(Batch_A)&!is.null(Batch_B)) {
      Part_A=subset(Merged,subset=object==Name_of_object_A)@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
      PCA.B=PCA.B@reductions$harmonized@cell.embeddings[,seq(1,numbersofPC)]
      colnames(PCA.B)=colnames(Part_A)
      harmonized.space=rbind(Part_A,PCA.B)
      
      
      Original.B=Original.B@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
      Part_A=subset(Merged,subset=object==Name_of_object_A)@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
      colnames(Original.B)=colnames(Part_A)
      original.space=rbind(Part_A,Original.B)
      
      
      
      
    }
    else if (!is.null(Batch_A)&!is.null(Batch_B)) {
      PCA.A=PCA.A@reductions$harmonized@cell.embeddings[,seq(1,numbersofPC)]
      PCA.B=PCA.B@reductions$harmonized@cell.embeddings[,seq(1,numbersofPC)]
      harmonized.space=rbind(PCA.A,PCA.B)
      
      Original.A=Original.A@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
      Original.B=Original.B@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
      
      original.space=rbind(Original.A,Original.B)
      
    }
    
  }
  else {
    harmonized.space=Merged@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
    original.space=Merged@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
    
  }
    return(list(harmonized.space=harmonized.space,
                original.space=original.space))
  }
  else {
    if (!is.null(Batch_A)|!is.null(Batch_B)) {
      if (!is.null(Batch_A)) {
        PCA.A=subset(Merged,subset=object==Name_of_object_A)
        if (Normalization_method=='SCT') {
          PCA.A=RunHarmony(PCA.A,group.by.vars=Batch_A,assay.use='SCT',max.iter.harmony = 1000,max.iter.cluster = 1000)
          
        }
        else {
          PCA.A=RunHarmony(PCA.A,group.by.vars=Batch_A,max.iter.harmony = 1000,max.iter.cluster = 1000)
          
        }
      }
      if (!is.null(Batch_B)) {
        PCA.B=subset(Merged,subset=object==Name_of_object_B)
        if (Normalization_method=='SCT') {
          PCA.B=RunHarmony(PCA.B,group.by.vars=Batch_B,assay.use='SCT',max.iter.harmony = 1000,max.iter.cluster = 1000)
          
        }
        else {
          PCA.B=RunHarmony(PCA.B,group.by.vars=Batch_B,max.iter.harmony = 1000,max.iter.cluster = 1000)
          
        }
      }
      
      if (!is.null(Batch_A)&is.null(Batch_B)) {
        PCA.A=PCA.A@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
        Part_B=subset(Merged,subset=object==Name_of_object_B)@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
        colnames(PCA.A)=colnames(Part_B)
        PCA.space=rbind(PCA.A,Part_B)
      }
      if (is.null(Batch_A)&!is.null(Batch_B)) {
        Part_A=subset(Merged,subset=object==Name_of_object_A)@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
        PCA.B=PCA.B@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
        colnames(PCA.B)=colnames(Part_A)
        PCA.space=rbind(Part_A,PCA.B)
      }
      if (!is.null(Batch_A)&!is.null(Batch_B)) {
        PCA.A=PCA.A@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
        PCA.B=PCA.B@reductions$harmony@cell.embeddings[,seq(1,numbersofPC)]
        PCA.space=rbind(PCA.A,PCA.B)
      }
      
    }
    else {
      PCA.space=Merged@reductions$pca@cell.embeddings[,seq(1,numbersofPC)]
      
    }
    return(PCA.space)
  }
}


Expression_correlation <- function(Merged,annotation_Name,object_Name,Name_of_object_A,Name_of_object_B,assay.use,genes.test) {
  require(Matrix.utils)
  require(muscat)
  require(scater)
  DefaultAssay(Merged)='RNA'
  genes.test=rownames(Merged)
  Merged=NormalizeData(Merged,verbose = F)
  Merged <- SingleCellExperiment(
    assays = list(counts = Merged@assays$RNA@data[genes.test,,,]), 
    colData = Merged@meta.data
  )
  Expression.matrix=aggregateData(Merged,by = annotation_Name,fun = 'sum')
  
  
  results_matrix=assay(Expression.matrix)
  results_matrix=results_matrix[rowSums(results_matrix)>0,]
  results_matrix=cor(results_matrix,method = 'spearman')
  results_matrix=(results_matrix+1)/2
  results_matrix=results_matrix[grepl(Name_of_object_A,rownames(results_matrix)),
                                grepl(Name_of_object_B,colnames(results_matrix))]
  return(results_matrix)
}



Expression_correlation_Randomized <- function(Matrix,cellsA,cellsB) {
  require(Matrix.utils)
  require(muscat)
  require(scater)
  Matrix=Matrix[,colnames(Matrix)%in%c(cellsA,cellsB)]
  CellsA.sum=rowSums(Matrix[,colnames(Matrix)%in%c(cellsA)])
  CellsB.sum=rowSums(Matrix[,colnames(Matrix)%in%c(cellsB)])
  
  
  results=cor(CellsA.sum,CellsB.sum,method = 'spearman')
  results=(results+1)/2
  return(results)
}