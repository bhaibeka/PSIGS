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
saveres <- file.path("res_drug")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

## libraries
library(Vennerable)
library(WriteXLS)
library(genefu)
library(minet)
library(sideChannelAttack)
library(GOSim)

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
library(mRMRe)