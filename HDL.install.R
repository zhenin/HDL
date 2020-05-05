##### Install HDL #####
if(!require("dplyr",character.only = TRUE)){
  try_error <- try(install.packages("dplyr", 
                                    dependencies=TRUE, 
                                    lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  if(!is.null(try_error)){
    try_error <- try(install.packages("dplyr", 
                                      dependencies=TRUE), silent = TRUE)
  }
}
try_error <- try(install.packages("HDL",
                                  repos = NULL,
                                  lib = Sys.getenv("R_LIBS_USER"),
                                  type = "source"), silent = TRUE)
if(!is.null(try_error)){
  try_error <- try(install.packages("HDL",
                                    repos = NULL,
                                    type = "source"), silent = TRUE)
}

