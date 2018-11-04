## Input Variables
# y is a matrix of RNA-Seq count data. rows are genes/features and columns are samples.
# x1 is the design matrix. This should NOT include the intercept.
# x2 is not functional yet. Leave as NULL. If NULL it will add the intercept to the design.
# cIdx is a character vector of the gene/feature ID's for the negative controls.
# k is the number of unwanted factors to be used
# scIdx is not used here. Leave as NULL. A value of NULL has no effect on the algorithm.
# method is a character vector that signifies the RUV method to be implemented. I have
#only given you code for the "RUVg4poi" and "RUVg4v" methods. You can also implement
#the methods "RUV4" and "RUV2" since these are already existing methods in the 'ruv' R package.
# ret.counts is only relevant for RUV2 based methods at the moment.
# W is the matrix of unwanted factors
# groups is a character vector indicating which groups each of the samples belong to. It should
#be the same length as the number of samples and any samples belonging to the same group
#should have the same group indicator. Default NULL
# integer or character vector indicating which coefficients of the linear model are to be tested
#equal to zero. Values must be columns or column names of design. Default 1 indicating 1st column of x1.
# fdr is the FDR level one wishes to acheive as a fraction (not percent).
# pIdx is a character vector of the gene/feature IDs of the positive controls, if any. Default NULL
# deMethod indicates which DE method should be carried out, "edgeR" or "voom". This
#is independent of the specific RUV method used. In other words, any RUV method
#to obtain W can be paired with any DE method. Default NULL which corresponds to "voom"

## Functions
# log1() defines a link function, log(y+1), used for the RUVg4.poi() function
# RUVg4.poi() implements the RUV4 methodology but within a poisson-GLM framework
#The data, y, should be filtered and have undergone some form of library scaling
#The matrix of unwanted factors, W, is returned
# RUVg4.voom() implements the RUV4 methodology and is exactly the same method
#as the original RUV4() function in the 'ruv' package except that it
#calculates weights for each observation and uses these weights within the
#factor analysis step, ie. PCA step. The original RUV4() method does not use
#weigts or rather each observation gets equal weights.
#The matrix of unwanted factors, W, is returned
# getW() is a wrapper for all the RUV functions that I've been testing and includes
#the two above. In addition, it also includes the original RUV4 and RUV2 methods
#and RUVg which you can also run if you load the 'ruv' package and 'RUVSeq' package.
#The other RUV methods seen in the default value of the "methods" variable you
#will not be able to run.
#NOTE!!! The y input should be only the filtered data, ie. the data should NOT have
#undergone a library scaling. getW() does its own library scaling and feeds that into
#the RUV methods.
#The matrix of unwanted factors, W, is returned
# getDE() is a wrapper for testing DE. Either "voom" or "edgeR" can supplied to the
#variable "deMethod". It requires a W matrix which can be NULL if no normalization
#is wanted.
#NOTE!!! The y input should be only the filtered data, ie. the data should NOT have
#undergone a library scaling. getDE() does its own library scaling via the internal DE algorithms.
#The output is a list with the first element being the pvalues for each gene.
#The second and third elements are only relevant if you have positive controls and
#wish to know their ranks

##Note:
# If you wish to have the RUVg4.poi and RUVg4.voom functions return normalized counts
# to use for biomarker detection (ie. a global normalization) then you'll have to set that up
# yourself, however, at the bottom of each of these functions is commented out code that
# should get you the normalized counts. However, it would be good for you to make sure that
# part of the code functions correclty as I'm not sure how long that has been commented out and
# things might have changed since commenting it out.

#Also, if you wish to have getDE() return additional information besides pvalues then you'll have to
#monkey with the code to get that out.

#Link Function for RUVg4.poi, i.e. log(y+1)
log1 <- function()
{
  ## link
  linkfun <- function(y) log(y+1)
  
  ## inverse link
  linkinv <- function(eta)  exp(eta) - 1
  
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) exp(eta)
  
  valideta <- function(eta) TRUE
  link <- "log(y+1)"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),
            class = "link-glm")
  
}

RUVg4.poi <- function (y, x1, x2=NULL, cIdx, k)
{
  
  coFun <- function( x , y, off=NULL)
  {
    if( is.null(off) )
    {
      return( sapply(1:ncol(y),FUN=function(i)
      {
        fit = tryCatch(glm(y[,i]~x-1,family=poisson(link=log1())),error=function(e) lm(log(y[,i]+1)~x-1),silent=TRUE)
        coefs = coef(fit)
        return( coefs )
      } ))
    }
    else return( sapply(1:ncol(y),FUN=function(i)
    {
      fit = tryCatch(glm(y[,i]~x-1,offset = off[,i],family=poisson(link=log1())),error=function(e) lm(log(y[,i]+1)~x-1,offset = off[,i]),silent=TRUE)
      coefs = coef(fit)
      return( coefs )
    }))
  }
  
  if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(y))
  if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
  if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))
  
  Y <- t(y)
  Yc <- Y[,cIdx]
  m <- nrow(Y)
  n <- ncol(Y)
  nc <- ncol(Yc)
  
  D = cbind(x1 , x2)
  co = coFun(D , Y)
  Z = log(Y+1) - D %*% co
  
  svdWa <- svd(Z)
  W_o <- svdWa$u[, (1:k), drop = FALSE]
  
  D = cbind(D, W_o)
  co = coFun(D , Y)
  colnames(co) <- colnames(Y)
  ac = co[(ncol(D)-ncol(W_o)+1):(ncol(D)-ncol(W_o)+k),cIdx,drop=FALSE]
  x2c = co[(ncol(x1)+1):(ncol(x1)+ncol(x2)),cIdx,drop=FALSE] #this needs to be more general as well for X2 matrices
  
  OFF = W_o %*% ac + x2 %*% x2c
  Z = matrix(as.vector(Yc),nc=1,nr=m*nc)
  OFF = matrix(as.vector(OFF),nc=1,nr=m*nc)
  x1.expand = apply( x1 , 1 , FUN=function(vec) rep(vec,each=k))
  if( is.null(dim(x1.expand)) ) x1.expand = matrix(x1.expand,nr=k)
  D = do.call("rbind", lapply( 1:ncol(ac), FUN=function(i) t(x1.expand * ac[,i]) ) )
  co = coFun(D, Z, OFF)
  if( is.null(dim(co)) ) co=matrix(co,nr=k)
  
  W = W_o + do.call("+",lapply(1:ncol(x1),FUN=function(i) x1[,i,drop=FALSE] %*% t(co[((i-1)*k+1):(i*k),])))
  
         D = cbind(x1,x2, W)
         co = coFun(D, Y)
         a <- co[(ncol(D)-ncol(W)+1):(ncol(D)-ncol(W)+k),,drop=FALSE]
         Ynorm <- matrix(pmax(round(exp(log(Y+1) - W %*% a) - 1,0),0),nr=nrow(Y),ncol=ncol(Y))
         rownames(Ynorm) <- rownames(Y)
         colnames(Ynorm) <- colnames(Y)
  
  return( list(W,Ynorm) )
}

RUVg4.voom <- function (y, x1, x2=NULL, cIdx, k)
{
  require(limma)
  require(scde)
  
  coFun <- function( x , y, off=NULL, weights = NULL , retWeights=FALSE )
  {
    if( is.null(off) )
    {
      dge <- DGEList(counts=y)
      dge <- calcNormFactors(dge)
      v=voom(dge,x)
      weights=t(v$weights)
      fit=lmFit(v,x)
      coefs=t(fit$coefficients)
      if( retWeights )
        return( rbind( coefs , weights ) )
      else return( coefs )
    }
    else
    {
      fit = lm(y~x-1,weights=weights,offset=OFF)
      return( coef(fit) )
    }
  }
  
  if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(y))
  if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
  if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))
  
  dge <- DGEList(counts=y)
  dge <- calcNormFactors(dge)
  lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
  Y <- log2(t(y+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
  Yc <- Y[, cIdx]
  m <- nrow(Y)
  n <- ncol(Y)
  nc <- ncol(Yc)
  
  D = cbind(x1 , x2)
  co = coFun(D , y, retWeights = TRUE)
  weights = co[(ncol(D)+1):nrow(co),]
  co = co[1:ncol(D),]
  Z = Y - D %*% co
  
  tmp = bwpca(Z,matw=weights,npcs=k,center=FALSE)
  W_o = apply(tmp$scores,2,FUN=function(vec) vec/sqrt(sum(vec^2)))
  
  D = cbind(D, W_o)
  co = coFun(D , y, retWeights=TRUE)
  ac = co[(ncol(D)-ncol(W_o)+1):(ncol(D)-ncol(W_o)+k),cIdx,drop=FALSE]
  x2c = co[(ncol(x1)+1):(ncol(x1)+ncol(x2)),cIdx,drop=FALSE]
  weights = co[(ncol(D)+1):nrow(co),]
  
  OFF = W_o %*% ac + x2 %*% x2c
  Z = matrix(as.vector(Yc),nc=1,nr=m*nc)
  Weights = matrix(as.vector(weights[,cIdx]),nc=1,nr=m*nc)
  OFF = matrix(as.vector(OFF),nc=1,nr=m*nc)
  x1.expand = apply( x1 , 1 , FUN=function(vec) rep(vec,each=k))
  if( is.null(dim(x1.expand)) ) x1.expand = matrix(x1.expand,nr=k)
  D = do.call("rbind", lapply( 1:ncol(ac), FUN=function(i) t(x1.expand * ac[,i]) ) )
  co = coFun( D , Z , OFF , Weights )
  if( is.null(dim(co)) ) co=matrix(co,nr=k)
  
  W = W_o + do.call("+",lapply(1:ncol(x1),FUN=function(i) x1[,i,drop=FALSE] %*% co[((i-1)*k+1):(i*k)]))
  
      D = cbind(x1, x2, W)
      co = coFun(D , y)
      a <- co[(ncol(D)-ncol(W)+1):(ncol(D)-ncol(W)+k),,drop=FALSE]
      Ynorm <- Y - W %*% a
      ynorm <- pmax(round(2^Ynorm*(lib.size+1)/1e+06-0.5,0),0)
  
      return( list(W,Ynorm) )
}

getW <- function( counts , cIdx, k , scIdx = NULL , x1=NULL, x2=NULL,
                  method = c("RUV2","RUVg2","RUVs","RUVg2v","RUVg2poi","RUV4","RUVg4poi","RUVg4v","RUVg4nb","ALL") , ret.counts=FALSE )
{
  require(ruv)
  require(RUVSeq)
  require(limma)
  
  if( k==0 ) return( NULL )
  #For getting W:
  #RUV2: voom transform. no weights. std linear model to intercept. Same as RUVg2 except for the tranform
  #RUVg2: log-e+1 transform. no weights. std linear model to intercept.
  #RUVs: log-e+1 transform. no weights. std linear model to intercept.
  #RUVg2v: voom transform. voom weights. weighted linear model to intercept. Same as RUV2 above but uses weights in the PCA step.
  #RUVg2poi: log-e+1 transform. no weights. Poi-GLM fit to intercept.
  #RUV4: voom transform. no weights. std linear model.
  #RUVg4poi: log-e+1 transform. no weights. poisson GLM.
  #RUVg4v: voom transform. voom weights. weighted linear model.
  #RUVg4nb: log-e+1 tranform. no weights. negative binomial GLM.
  trcounts <- function( counts )
  {
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge)
    lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
    counts <- log2((counts+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
    return(counts)
  }
  method = match.arg(method,c("RUV2","RUVg2","RUVs","RUVg2v","RUVg2poi","RUV4","RUVg4poi","RUVg4v","RUVg4nb","ALL"))
  cIdx2 <- rownames(counts)%in%cIdx
  uq <- betweenLaneNormalization(counts, which="upper")
  trcts <- trcounts(counts)
  W = switch(method, RUV2 = RUV2(Y=t(trcts), X=x1, ctl=cIdx2, k=k, do_projectionplot=FALSE)$W,
             RUVg2 = RUVg(x=uq,cIdx=cIdx,k=k)$W,
             RUVs = RUVs(x=uq, cIdx=cIdx, k=k, scIdx=scIdx)$W,
             RUVg2v = RUVg2.voom(y=counts, x1=x1, x2=x2, cIdx=cIdx, k=k, ret.counts=ret.counts),
             RUVg2poi = RUVg2.poi( y=uq, x1=x1, cIdx=cIdx, k=k, ret.counts=ret.counts) ,
             RUV4 = RUV4(Y=t(trcts), X=x1, ctl=cIdx2, k=k)$W,
             RUVg4poi = RUVg4.poi(y=uq, x1=x1, cIdx=cIdx, k=k),
             RUVg4v = RUVg4.voom(y=counts, x1=x1, cIdx=cIdx, k=k),
             RUVg4nb = RUVg4.nb(y=uq, x1=x1, cIdx=cIdx, k=k) ,
             ALL = list("RUV2" = RUV2(Y=t(trcts), X=x1, ctl=cIdx2, k=k, do_projectionplot=FALSE)$W,
                        "RUVg2" = RUVg(x=uq,cIdx=cIdx,k=k)$W,
                        "RUVs" = RUVs(x=uq, cIdx=cIdx, k=k, scIdx=scIdx)$W,
                        "RUVg2v" = RUVg2.voom(y=counts, x1=x1, x2=x2, cIdx=cIdx, k=k, ret.counts=ret.counts),
                        "RUVg2poi" = RUVg2.poi( y=uq, x1=x1, cIdx=cIdx, k=k, ret.counts=ret.counts),
                        "RUV4" = RUV4(Y=t(trcts), X=x1, ctl=cIdx2, k=k)$W,
                        "RUVg4poi" = RUVg4.poi(y=uq, x1=x1, cIdx=cIdx, k=k),
                        "RUVg4v" = RUVg4.voom(y=counts, x1=x1, cIdx=cIdx, k=k),
                        "RUVg4nb" = RUVg4.nb(y=uq, x1=x1, cIdx=cIdx, k=k)))
  return(W)
}


getDE <- function(counts,x1,x2=NULL,W,groups=NULL,fdr=0.05,coef=1,pIdx=NULL,method=c("RUV2","RUVg2","RUVs","RUVg2v","RUVg2poi","RUV4","RUVg4poi","RUVg4v","RUVg4nb","ALL"), deMethod=NULL)
{
  require(edgeR)
  require(limma)
  method = match.arg(method,c("RUV2","RUVg2","RUVs","RUVg2v","RUVg2poi","RUV4","RUVg4poi","RUVg4v","RUVg4nb","ALL"))
  eRFun <- function(counts, x1, x2, W, groups, fdr, coef)
  {
    D=cbind(x1,x2,W)
    y <- DGEList(counts = counts, group = groups)
    y <- calcNormFactors(y, method = "upperquartile")
    y <- estimateGLMCommonDisp(y, D)
    y <- estimateGLMTagwiseDisp(y, D)
    fit <- glmFit(y, D)
    lrt <- glmLRT(fit, coef=coef)
    tbl <- topTags(lrt, n=Inf)$table
    pvals <- tbl$PValue
    names(pvals) <- rownames(tbl)
    de <- rownames(tbl[tbl$FDR<=fdr,])
    posDE=NULL
    posAll=NULL
    if( !is.null(pIdx) )
    {
      posDE <- match(pIdx,de)
      posDE[is.na(posDE)] <- Inf
      posDE <- sort(posDE)
      posAll <- match(pIdx,names(pvals))
      posAll[is.na(posAll)] <- Inf
      posAll <- sort(posAll)
      posAll[which(posDE==Inf)] <- -1*posAll[which(posDE==Inf)]
    }
    pvals <- pvals[rownames(counts)]
    return( list( "Pvals" = pvals , "PosDE" = posDE , "PosRank" = posAll ) )
  }
  vFun <- function(counts, x1, x2, W, fdr, coef)
  {
    D=cbind(x1,x2,W)
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge)
    v=voom(dge,D)
    fit=lmFit(v,D)
    fit2 = eBayes(fit)
    tbl <- topTable(fit2,coef=coef,number=nrow(counts),sort.by="P")
    pvals <- tbl$P.Value
    names(pvals) <- rownames(tbl)
    de <- rownames(tbl[tbl$adj.P.Val<=fdr,])
    posDE=NULL
    posAll=NULL
    if( !is.null(pIdx) )
    {
      posDE <- match(pIdx,de)
      posDE[is.na(posDE)] <- Inf
      posDE <- sort(posDE)
      posAll <- match(pIdx,names(pvals))
      posAll[is.na(posAll)] <- Inf
      posAll <- sort(posAll)
      posAll[which(posDE==Inf)] <- -1*posAll[which(posDE==Inf)]
    }
    pvals <- pvals[rownames(counts)]
    return( list( "Pvals" = pvals , "PosDE" = posDE , "PosRank" = posAll ) )
  }
  
  if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(counts))
  if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
  if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))
  
  if( deMethod=="edgeR" )
    out <- eRFun(counts=counts, x1=x1, x2=x2, W=W, groups=groups, fdr=fdr, coef=coef)
  else if( deMethod=="voom" )
    out <- vFun(counts=counts, x1=x1, x2=x2, W=W, fdr=fdr, coef=coef)
  
  return(out)
}


