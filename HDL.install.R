##### Install HDL #####
if(!require("dplyr",character.only = TRUE)){
  try_error <- try(install.packages("dplyr", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("dplyr", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}
if(!require("data.table",character.only = TRUE)){
  try_error <- try(install.packages("data.table", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("data.table", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}
data.table.version <- packageVersion("data.table")
if(data.table.version < "1.12.1"){
  message("Package 'data.table' needs to updated...")
  try_error <- try(install.packages("data.table", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("data.table", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}
try_error <- try(install.packages("HDL",
                                  repos = NULL,
                                  type = "source"), silent = TRUE)
# if(!is.null(try_error)){
#   try_error <- try(install.packages("HDL",
#                                     repos = NULL,
#                                     lib = Sys.getenv("R_LIBS_USER"),
#                                     type = "source"), silent = TRUE)
# }

