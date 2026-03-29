# Suppress "no visible binding for global variable" NOTEs from R CMD check.
# These variables are column names in data frames created within plot functions
# and referenced inside ggplot2::aes(), which uses non-standard evaluation.
utils::globalVariables(c("RCP", "Method", "Approach", "Scale", "N_label"))
