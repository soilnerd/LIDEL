######################
#Main program to run post processing for chains
#and convergence diagnostics
#
#Nell Campbell
#created 10/14/15
######################

rm(list=ls())
library(deSolve)
library(gplots)
library(pscl)
library(rootSolve)
library(coda)
library(compiler)

#Model 3 post processing
  #chain 1
    source("C:/LIDEL/Model3/Chain1/LIDEL_post_processing.R")
  #chain 2
    source("C:/LIDEL/Model3/Chain2/LIDEL_post_processing2.R")
  #chain 3
    source("C:/LIDEL/Model3/Chain3/LIDEL_post_processing3.R")
