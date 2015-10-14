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

#Model 1 post processing
  #chain 1
    source("C:/LIDEL/Model1/Chain1/LIDEL_post_processing.R")
  #chain 2
    source("C:/LIDEL/Model1/Chain2/LIDEL_post_processing2.R")
  #chain 3
    source("C:/LIDEL/Model1/Chain3/LIDEL_post_processing3.R")

#Model 2 post processing
  #chain 1
    source("C:/LIDEL/Model2/Chain1/LIDEL_post_processing.R")
  #chain 2
    source("C:/LIDEL/Model2/Chain2/LIDEL_post_processing2.R")
  #chain 3
    source("C:/LIDEL/Model2/Chain3/LIDEL_post_processing3.R")

#Model 3 post processing
  #chain 1
    source("C:/LIDEL/Model3/Chain1/LIDEL_post_processing.R")
  #chain 2
    source("C:/LIDEL/Model3/Chain2/LIDEL_post_processing2.R")
  #chain 3
    source("C:/LIDEL/Model3/Chain3/LIDEL_post_processing3.R")

#Model 4 post processing
  #chain 1
    source("C:/LIDEL/Model4/Chain1/LIDEL_post_processing.R")
  #chain 2
    source("C:/LIDEL/Model4/Chain2/LIDEL_post_processing2.R")
  #chain 3
    source("C:/LIDEL/Model4/Chain3/LIDEL_post_processing3.R")

