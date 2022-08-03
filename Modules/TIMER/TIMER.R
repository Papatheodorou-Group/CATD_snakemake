## This code is adapted from from http://cistrome.org/TIMER/download.html.
## The method is described in Li et al. Genome Biology 2016;17(1):174., PMID 27549193.

## TIMER Pipeline for analyzing immune cell components in the tumor microenvironment

getPurityGenes <- function(dd,purity,thr.p=0.05,thr.c=0,mode='env'){
  tmp.dd=as.matrix(dd)
  tmp=lapply(rownames(tmp.dd),function(x)cor.test(as.numeric(tmp.dd[x,]),purity,method='s'))
  tmp.pp=sapply(tmp,function(x)x$p.value)
  tmp.cor=sapply(tmp,function(x)x$estimate)
  names(tmp.pp)=names(tmp.cor)=rownames(tmp.dd)
  if(mode=='env')vv=names(which(tmp.pp <=thr.p&tmp.cor < thr.c))
  if(mode=='tumor')vv=names(which(tmp.pp <=thr.p&tmp.cor > thr.c))
  return(vv)
}

##----- Combine TCGA gene expression profiles with the selected reference data, remove batch effect and aggregate samples of each immune category by taking the median -----##
RemoveBatchEffect <- function(dd,curated.ref.genes,curated.cell.types){
  library(sva)
  tmp.dd=as.matrix(dd)
  tmp.ss=intersect(rownames(tmp.dd),rownames(curated.ref.genes))
  N1=ncol(tmp.dd) # N1 is the number of mixture samples
  tmp.dd=cbind(tmp.dd[tmp.ss,],curated.ref.genes[tmp.ss,]) # Concatenate mixture profiles and reference profiles
  tmp.dd=as.matrix(tmp.dd)
  mode(tmp.dd)='numeric'
  N2=ncol(curated.ref.genes) # N2 is the number of samples in the reference
  tmp.batch=c(rep(1,N1),rep(2,N2)) # take the mixture samples and reference profiles as batch
  tmp.dd0=ComBat(tmp.dd,tmp.batch,c()) # uses ComBat function for batch correction
  dd.br=tmp.dd0[,1:N1] # mixture profiles after batch correction
  curated.ref.genes.br=tmp.dd0[,(N1+1):(N1+N2)] # reference profiles after batch correction
  tmp0=c()
  for(kk in unique(names(curated.cell.types))){
    tmp.vv=which(names(curated.cell.types)==kk)
    tmp0=cbind(tmp0,apply(curated.ref.genes.br[,tmp.vv],1,median,na.rm=T))
  } # takes the median value to collapse reference profile to have 1 profile per cell type
  curated.ref.genes.agg.br=tmp0
  colnames(curated.ref.genes.agg.br)=unique(names(curated.cell.types))
  #rownames(curated.ref.genes.agg.br)=rownames(curated.ref.genes.br)
  return(list(dd=dd.br,rr=curated.ref.genes.br,rrg=curated.ref.genes.agg.br))
}

RemoveOutliers <- function(vv, ref.dd, thr.q=0.99){
  ## removes upper thr.q quantile for every reference feature
  remove.vv=c()
  for(i in 1:ncol(ref.dd)){
    tmp=quantile(ref.dd[vv,i],thr.q)[1]
    tmp.vv=which(ref.dd[vv,i]>tmp)
    remove.vv=c(remove.vv,tmp.vv)
  }
  remove.vv=unique(remove.vv)
  return(vv[-remove.vv])
}

##----- Constrained regression method implemented in Abbas et al., 2009 -----##
getFractions.Abbas <- function(XX,YY,w=NA){
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(rep(0,ncol(XX)))
      tmp.XX=tmp.XX[,-ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0)return(rep(0,ncol(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX,YY,intercept=F) else tmp=lsfit(tmp.XX,YY,w,intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0,ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}

##----- function to process deconvolution method in batch -----##
BatchFractions <- function(XX,YYd){
  Fmat=c()
  for(i in 1:ncol(YYd)){
    YY=YYd[,i]
    tmp.F=getFractions.Abbas(XX,YY)
    #tmp.F=getFractions.Optim(XX,YY)
    Fmat=rbind(Fmat,tmp.F)
  }
  rownames(Fmat)=colnames(YYd)
  colnames(Fmat)=colnames(XX)
  return(Fmat)
}

TIMER_deconv <- function(mix,ref,curated.cell.types,sig){
  ##### Step 1: Batch correction
  tryCatch({
    tmp <- RemoveBatchEffect(mix,ref, curated.cell.types)
    dd.br=tmp$dd # mix profiles after batch correction
    curated.ref.genes.br=tmp$rr # regular ref profiles after batch correction
    curated.ref.genes.agg.br=tmp$rrg # collapsed ref
    ##### Step 2: Eliminate outlier genes
    id = RemoveOutliers(sig, curated.ref.genes.agg.br)
    ##### Step 3: Perform batch deconvolution
    XX=curated.ref.genes.agg.br[id,]
    YYd=dd.br[id,]
    Fmat=BatchFractions(XX,YYd)
    return(Fmat)}, error = function(e) message(e))
}
