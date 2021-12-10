#initializing rptfast package
install.packages("devtools")
library(devtools)
#create package
create_package("/Users/sophiehoffman/Desktop/ggf.fast")
use_mit_license()
use_package("regentrans")
use_package("Rcpp")
