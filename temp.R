## Merge features with ABC table

# save.image("merge_e2g_candidates.rda")
# stop()

# required packages
suppressPackageStartupMessages({
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(argparse)
library(dplyr)
})

# Function to verify that the required columns are in the data
verify_columns <- function(data, required_columns) {
  missing_cols <- setdiff(required_columns, colnames(data))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing: ", paste(missing_cols, collapse = ", "))
  }
}

## Define functions --------------------------------------------------------------------------------

# function to merge features to e2g_candidates table
merge_feature_to_e2g_candidates <- function(e2g_candidates, features, feature_score_cols, new_feature_names, agg_fun, fill_value) {
  
  # convert fill value to numeric
  fill_value <- suppressWarnings(as.numeric(fill_value))
  
  # Confirm desired features are available
  available_features <- feature_score_cols[feature_score_cols %in% names(features)]
  missing_features <- setdiff(feature_score_cols, available_features)
  if (length(missing_features) > 0) {
    warning("These features from config are missing in input: ", paste(missing_features, collapse = ", "))
  }

  # fill blanks or NA in new_feature_names with the original names
  missing_mask <- is.na(new_feature_names) | new_feature_names == ""
  new_feature_names[missing_mask] <- feature_score_cols[missing_mask]

  # Only retain relevant columns from feature
  mapping <- setNames(new_feature_names, feature_score_cols)

  # Rename cols in features table
  cols_to_rename <- intersect(names(mapping), names(features))
  features_renamed <- rename_with(features, .cols = all_of(cols_to_rename), .fn = function(old) {
    mapped <- mapping[old]
    replace(old, !is.na(mapped), unname(mapped))
  })

  # only retain relevant columns from feature
  available_features_renamed <- new_feature_names[new_feature_names %in% names(features_renamed)]
  feature <- select(features_renamed, ElementChr, ElementStart, ElementEnd, GeneEnsemblID, all_of(available_features_renamed))

  # create GRanges for e2g_candidates  and feature E-G pairs
  e2g_candidates_gr <- with(e2g_candidates, GRanges(seqnames = paste0(ElementChr, ":", GeneEnsemblID),
                              ranges = IRanges(ElementStart, ElementEnd)))
  
  feat_gr <- with(feature, GRanges(seqnames = paste0(ElementChr, ":", GeneEnsemblID),
                                   ranges = IRanges(ElementStart, ElementEnd)))
  
  # set same seqlevels for both GRanges objects to avoid warnings
  seqlevels_all_pairs <- as.character(unique(c(seqnames(e2g_candidates_gr), seqnames(feat_gr))))
  seqlevels(e2g_candidates_gr) <- seqlevels_all_pairs
  seqlevels(feat_gr) <- seqlevels_all_pairs
  
  # find overlaps between e2g_candidates and features
  ovl <- findOverlaps(e2g_candidates_gr, feat_gr)
  
  # merge e2g_candidates table with features
  merged <- cbind(e2g_candidates[queryHits(ovl)], feature[subjectHits(ovl), available_features_renamed, with = FALSE])
  
  # aggregate feature scores, where multiple features overlapped one e2g_candidates enhancer
  agg_cols <- setdiff(colnames(merged), available_features_renamed)
  merged <- aggregate_features(merged, feature_cols = available_features_renamed,
                               agg_cols = agg_cols, agg_fun = agg_fun)
  
  # report how many e2g_candidates pairs overlap >= 1 features
  perc_overlap <- nrow(merged) / nrow(e2g_candidates)
  if (perc_overlap == 0) {
    stop("No overlapping E2G pairs were found: features can not be merged")
  }
  message("Features overlapping E2G candidates: ", round(perc_overlap * 100, digits = 2), "%, merged ", nrow(merged), " pairs")
  
  # get pairs from e2g_candidates table that are missing from features 
  missing <- e2g_candidates[setdiff(seq_len(nrow(e2g_candidates)), queryHits(ovl)), ]
  
  # # # fill in missing value for features in missing pairs
  # new_cols <- as.data.frame(
  # lapply(available_features_renamed, function(x) rep(fill_value, nrow(missing))),
  # stringsAsFactors = FALSE
  # )
  # names(new_cols) <- available_features_renamed
  # # create/overwrite columns in missing
  # missing[available_features_renamed] <- new_cols

  if (nrow(missing)>0) {
    for (i in available_features_renamed) {
      missing[[i]] <- fill_value
    }
  }
  
  # combine merged and missing pairs to create output
  output <- rbindlist(list(merged, missing), fill = TRUE)
  
  # sort output by cre position
  output <- output[order(ElementChr, ElementStart, ElementEnd, GeneEnsemblID), ]
  
  return(output)
  
}

# function to aggregate multiple overlaps
aggregate_features <- function(merged, feature_cols, agg_cols, agg_fun) {
  
  # aggregate feature columns as specified by agg_fun
  agg_list <- mapply(FUN = function(feat_col, agg_func, agg_cols) {
    agg_func <- get(agg_func)
    merged[, setNames(.(agg_func(get(feat_col))), feat_col), by = agg_cols]
  }, feat_col = feature_cols, agg_func = agg_fun, MoreArgs = list(agg_cols = agg_cols),
  SIMPLIFY = FALSE)
  
  # merge all the aggregates together to make collapsed data.frame
  output <- as.data.table(Reduce(function(df1, df2) merge(df1, df2, by = agg_cols), agg_list))
  
  return(output)
  
}

## Process features --------------------------------------------------------------------------------

message("Loading data...")

parser <- ArgumentParser()

parser$add_argument(
  "-e",
  "--e2g-candidates",
  type = "character",
  default = NULL,
  help = "Path to E2G candidate universe",
  metavar = "character"
)
parser$add_argument(
  "-f",
  "--features",
  type = "character",
  default = NULL,
  help = "Path to E2G feature table",
  metavar = "character"
)
parser$add_argument(
  "-c",
  "--config",
  type = "character",
  metavar = "character",
  default = NULL,
  help = "Path to feature merge configuration table"
)
parser$add_argument(
  "-o",
  "--output",
  type = "character",
  default = NULL,
  help = "Name of merged feature data output file",
  metavar = "character"
)

# Get command line data
args <- parser$parse_args()

# load input files
e2g_candidates <- fread(args$e2g_candidates)
features <- fread(args$features, header = TRUE, sep = "\t")
config <- fread(args$config)

# Validate the required columns
required_e2g_cols <- c("ElementChr", "ElementStart", "ElementEnd", "GeneEnsemblID")
verify_columns(e2g_candidates, required_e2g_cols)
verify_columns(features, required_e2g_cols)
verify_columns(config, c("feature_col", "feature_rename", "aggregation_function"))

# Merge features with e2g_candidates
output <- merge_feature_to_e2g_candidates(
  e2g_candidates,
  features,
  feature_score_cols = config$feature_col,
  new_feature_names = config$feature_rename,
  agg_fun = config$aggregation_function,
  fill_value = NA_real_
)

# Write output to file
fwrite(output, file = args$output, sep = "\t")