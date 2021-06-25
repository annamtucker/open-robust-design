
pkgs = c("tidyverse", "cowplot", "lubridate", 
         "tidybayes", "nimble", "broom")

check <- sapply(pkgs, library,
                warn.conflicts = TRUE,
                character.only = TRUE)

# if(any(!check)){
#   pkgs.missing <- pkgs[!check]
#   install.packages(pkgs.missing)
#   check <- sapply(pkgs.missing,require,warn.conflicts = TRUE,character.only = TRUE)
# }

source("scripts/00_functions.R")

theme_set(theme_cowplot())