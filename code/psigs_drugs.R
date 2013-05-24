#################################################
## Probabilistic gen signatures
##
## Benjamin Haibe-Kains
##
## R code under license Artistic-2.0
## 2013-05-24
#################################################


## remove all existing objects to start from a clean R session
rm(list=ls())

## directory containing all the (intermediate) analysis results
saveres <- file.path("saveres")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

## "bad" characters
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

#################################################
## libraries
#################################################

# library(Vennerable)
suppressPackageStartupMessages(require(WriteXLS)) || stop("Library WriteXLS is not available!")
suppressPackageStartupMessages(require(genefu)) || stop("Library genefu is not available!")
# library(GOSim)

## number of cores to use during normalization
suppressPackageStartupMessages(require(parallel)) || stop("Library parallel is not available!")
nbcore <- 8
availcore <- detectCores()
if(nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)

## load library mRMRe
## install the development version of the package if needed
# library(devtools)
# install_github("mRMRe", username="bhaibeka", ref="master")
suppressPackageStartupMessages(require(mRMRe)) || stop("Library mRMRe is not available!")


#################################################
## parameters
#################################################

## number of signatures
solnb <- 10

## signature length
solength <- 100

## method to use for the ensemble mRMR (either "exhaustive" or "bootstrap")
fs.method <- "exhaustive"

## which drug
drugnoi <- "drugid_PD0325901"

## which tissue
# tissuen <- list("all", "each", ...)
# tissuen <- list("each", "all")
tissuen <- list("LUNG"="lung")
## should tissue type be ignored?
tissue.ignore <- FALSE


#################################################
## load pharmacogenomic data
#################################################

load(file.path("data", "data_all.RData"))


#################################################
## select the tissue(s) of interest
#################################################

if(tissue.ignore) {
  tissue.all <- lapply(tissue.all, function(x) {
    xx <- rep("ignoretissues", length(x))
    names(xx) <- names(x)
    return(xx)
  })
} 
utissue <- sort(unique(as.character(unlist(tissue.all))))
names(utissue) <- toupper(gsub(badchars, "", utissue))
tissuex <- NULL
for(i in 1:length(tissuen)) {
  switch(tissuen[[i]],
    "all"={
      tissuex <- c(tissuex, list("ALLTISSUES"=utissue))
    },
    "each"={
      tt <- as.list(utissue)
      names(tt) <- names(utissue)
      tissuex <- c(tissuex, tt)
    },
    {
      tissuex <- c(tissuex, tissuen[i])
      if(!is.null(names(tissuen[i]))) { names(tissuex)[length(tissuex)] <- toupper(gsub(badchars, "", tissuen[[i]])) }
    }
  )
}
## check whether names are unique
if(sum(duplicated(names(tissuex))) > 0) { names(tissuex) <- genefu::rename.duplicate(x=names(tissuex)) }


#################################################
## select the tissue(s) of interest
#################################################

for(i in 1:length(tissuex)) {
  ## tissue of interest
  tissue.ix <- tissuex[[i]]
  message(sprintf("\tfor tissue type %s", names(tissuex)[i]))
  ## gene expression data
  datat.all <- mapply(function(x, y, z) {
    myx <- !is.na(y) & is.element(as.character(y), as.character(z))
    return(x[myx, , drop=FALSE]) }, x=data.all, y=tissue.all, MoreArgs=list(z=tissue.ix), SIMPLIFY=FALSE)
  ## drugpheno
  drugphenot.all <- mapply(function(x, y, z) {
    myx <- !is.na(y) & is.element(as.character(y), as.character(z))
    return(x[myx, , drop=FALSE]) }, x=drugpheno.all, y=tissue.all, MoreArgs=list(z=tissue.ix), SIMPLIFY=FALSE)
  ## sampleinfo
  sampleinfot.all <- mapply(function(x, y, z) {
    myx <- !is.na(y) & is.element(as.character(y), as.character(z))
    return(x[myx, , drop=FALSE]) }, x=sampleinfo.all, y=tissue.all, MoreArgs=list(z=tissue.ix), SIMPLIFY=FALSE)
  ## batch
  batcht.all <- mapply(function(x, y, z) {
    myx <- !is.na(y) & is.element(as.character(y), as.character(z))
    return(x[myx]) }, x=batch.all, y=tissue.all, MoreArgs=list(z=tissue.ix), SIMPLIFY=FALSE) 
  ## tissue
  tissuet.all <- mapply(function(x, z) {
    myx <- !is.na(x) & is.element(as.character(x), as.character(z))
    return(x[myx]) }, x=tissue.all, MoreArgs=list(z=tissue.ix), SIMPLIFY=FALSE)
  ## subset the data for specicifc tissue type
  if(tissue.ignore) {
    tissue.all <- lapply(tissue.all, function(x) {
      xx <- rep("ignoretissues", length(x))
      names(xx) <- names(x)
      return(xx)
    })
    tissue.ix <- "ignoretissues"
  }
    
  ## check if enough tissue types
  if(any(sapply(datat.all, nrow) > 3)) {
    ## create a directory to store tissue-specific analysis results
    saverest <- file.path(saveres, sprintf("res_%s", names(tissuex)[i]))
    if(!file.exists(saverest)) { dir.create(saverest, showWarnings=FALSE) } 
   
    myfnall <- file.path(saverest, sprintf("drug_%s_%s.RData", gsub("drugid_", "", drugnoi), names(tissuex)[i]))
    if(!file.exists(myfnall)) {
      ## for each dataset
      myfs <- ens.m <- ens.systime <- NULL
      for(kk in 1:length(datat.all)) {
        myfn <- file.path(saverest, sprintf("drug_%s_%s_%s.RData", names(data.all)[kk], gsub("drugid_", "", drugnoi), names(tissuex)[i]))
        if(!file.exists(myfn)) {
          ## set the number of core to use with mRMRe
          set.thread.count(nbcore)
          ## create mRMRe.Data object
          dataset <- data.frame(drugphenot.all[[kk]][ , drugnoi], datat.all[[kk]])
          colnames(dataset)[1] <- drugnoi
          mdata <- mRMR.data(data=dataset)
          ## ensemble mRMR feature selection
          ens.systimet <- system.time(ens.mt <- mRMR.ensemble(mdata, solution_count=solnb, feature_count=solength, target_indices=1, method=fs.method))
          myfst <- apply(solutions(ens.mt)[[1]], 2, function(x, y) { return(y[x]) }, y=featureNames(mdata))
          save(list=c("ens.mt", "myfst", "ens.systimet"), compress=TRUE, file=myfn)
        } else { load(myfn) }
        myfs <- c(myfs, list(myfst))
        ens.m <- c(ens.m, list(ens.mt))
        ens.systime <- c(ens.systime, ens.systimet)
      }
      names(myfs) <- names(ens.m) <- names(ens.systime) <- names(data.all)
      save(list=c("myfs", "ens.m", "ens.systime"), compress=TRUE, file=myfnall)
    } else { load(myfnall) }
    ## save the signatures
    xx <- NULL
    for(jj in 1:length(myfs)) {
      xx <- c(xx, list(data.frame("sol"=myfs[[jj]])))
    }
    names(xx) <- names(myfs)
    WriteXLS(x="xx", ExcelFileName=file.path(saverest, sprintf("drug_%s_%s.xls", gsub("drugid_", "", drugnoi), names(tissuex)[i])), BoldHeaderRow=TRUE)
  } else { message("\t\t-> Not enough cell lines...") }
  
}