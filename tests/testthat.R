# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(leapR)
data(ncipid)
data(shortlist)
data(longlist)

##download file once for all tests
pdata <- download.file("https://api.figshare.com/v2/file/download/56536217", method = "libcurl", destfile = "protData.rda")
load("protData.rda")
p <- file.remove("protData.rda")

test_check("leapR")
