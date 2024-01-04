library(rxode2)
library(posologyr)
library(testthat)
setRxThreads(2L) #for CRAN, following the advice of mattfidler

test_check("posologyr")
