##### Install HDL #####
if(!"devtools" %in% rownames(installed.packages())){
    install.packages("devtools",dependencies=TRUE)
  }
library(devtools)
install_github("zhenin/HDL/HDL")


