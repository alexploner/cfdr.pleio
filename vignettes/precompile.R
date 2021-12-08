# Precompiled vignette that is
#  a. slow to build
#  b. downloads non-trivial amounts of data from other people's server
#
# This hack avoids the issues for package users while maintaining
# standard tool chain for buildding the packe with vignettes
#
# Must manually move image files from cfdr.pleio/figure to
#  cfd.pleio/vignettes/figure after knit
#
# # As seen at https://ropensci.org/blog/2019/12/08/precompute-vignettes/

library(knitr)
knit(input = "vignettes/Introduction.Rmd.orig", output = "vignettes/Introduction.Rmd")
