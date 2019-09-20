HDL.package.name <- list.files()[grep(pattern = glob2rx("HDL*.tar.gz"), x = list.files())]
##### Install HDL #####
if(!require("dplyr",character.only = TRUE)){
    install.packages("dplyr",dependencies=TRUE)
  }

install.packages(HDL.package.name, dependencies=TRUE, repos = NULL)