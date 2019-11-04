#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
CLI_Template.R (-h | --help | --version)
CLI_Template.R DIR
Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.
Arguments:
DIR    Provide directory where scyttools.args.Rdata file is located
" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

cat("\nLoading arguments from", ARGS_DIR, "\n")

load(paste(ARGS_DIR, "scyttools.args.Rdata", sep = ""))

RESULTS_DIR <- args$OUT

source("scyttools_functions.R")

##########################################################################
############################ R code goes here ############################
##########################################################################

# your R code!

##########################################################################
############################     End code     ############################
##########################################################################