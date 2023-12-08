data.table::setDTthreads(2) #throttle data.table for CRAN

library(testthat)
library(posologyr)

test_check("posologyr")
