##### Install HDL #####
if(!require("dplyr",character.only = TRUE)){
    install.packages("dplyr",dependencies=TRUE, lib = Sys.getenv("R_LIBS_USER"))
  }
install.packages("HDL", 
                 repos = NULL, 
                 lib = Sys.getenv("R_LIBS_USER"), 
                 type = "source")

