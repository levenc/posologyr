library(rxode2)
library(posologyr)
library(testthat)
setRxThreads(1L) #for CRAN, following the advice of mattfidler

test_check("posologyr")
