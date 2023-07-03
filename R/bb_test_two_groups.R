#!/usr/bin/env Rscript

if(!require(countdata)) {install.packages("countdata", dependencies = TRUE, quiet = TRUE, repos="https://cloud.r-project.org/")}
if(!require(stats)) {install.packages("stats", dependencies = TRUE, quiet = TRUE, repos="https://cloud.r-project.org/")}
if(!require(optparse)) {install.packages("optparse", dependencies = TRUE, quiet = TRUE, repos="https://cloud.r-project.org/")}

# Parse Argeuments
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, help="Input-CSV-File (comma seperated)", metavar="character"),
    make_option(c("-o", "--out"), type="character", default=NULL, help="output-CSV-File (will be comma seperated) with addtional columns at the end", metavar="character"),
    make_option(c("-A", "--groupA"), type="character", default="GroupA1,GroupA2", help="Columns which belong to group A. Entries should be comma-seperated", metavar="character"),
    make_option(c("-B", "--groupB"), type="character", default="GroupB1,GroupB2", help="Columns which belong to group A. Entries should be comma-seperated", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Open csv_file
inf <- read.csv(opt$file, header=TRUE, check.names = FALSE)

# Set groups to list
gA <- unlist(strsplit(opt$groupA, ","))
gB <- unlist(strsplit(opt$groupB, ","))

# Get the subset dataframes
gAdf <- subset(inf, select=gA)
gBdf <- subset(inf, select=gB)

# Set Groups and whole dataframe
groups <- c(rep("GroupA", length(gA)), rep("GroupB", length(gB))   )
gABdf <- cbind(gAdf, gBdf)
tx_sum <- colSums(gABdf)

# Execute bb test
bb_test_out <- bb.test(gABdf, tx_sum, groups)

# Export results to table
write.table(cbind(inf, bb_test_out$p.value), file = opt$out, sep = ",", row.names = FALSE)
