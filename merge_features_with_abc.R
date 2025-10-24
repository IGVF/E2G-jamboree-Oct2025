## Merge features with ABC table

# save.image("merge_e2g_candidates.rda")
# stop()

# required packages
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(argparse)
library(dyplr)

## Define functions --------------------------------------------------------------------------------

# function to merge features to e2g_candidates table
merge_feature_to_e2g_candidates <- function(e2g_candidates, feature, feature_score_cols, new_feature_names, agg_fun, fill_value) {
  
  # convert fill value to numeric
  fill_value <- suppressWarnings(as.numeric(fill_value))

  # confirm desired features are available
  wanted <- feature_score_cols
  available_features <- wanted[wanted %in% names(features)]
  missing <- setdiff(wanted, available_features)
  if (length(missing) > 0) warning("These features from config are missing in input: ", paste(missing, collapse = ", "))

  # E2G region columns
  e2g_region_cols <- c("ElementChr", "ElementStart", "ElementEnd", "GeneEnsemblID")

  # corresponding new names for the present columns
  available_new_feature_names <- new_feature_names[match(present_score_cols, feature_score_cols)]

  # select in config order then rename
  feature_renamed <- feature %>%
  select(all_of(e2g_region_cols), all_of(available_features)) %>%
  rename_with(.cols = all_of(available_features),
              .fn = function(old) {
                # map each old name to its new name using config order
                available_new_feature_names[match(old, available_features)]
              })
  feature_renamed <- feature %>%
  select(ElementChr, ElementStart, ElementEnd, GeneEnsemblID, all_of(feature_score_cols)) %>%
  rename_with(.cols = everything(),
              .fn = function(old) {
                # map old names to new names using config order
                new_feature_names[match(old, feature_score_cols)]
              })
  if (any(duplicated(names(df_out)))) warning("Duplicate target names in feature_rename: ", paste(unique(names(feature_renamed)[duplicated(names(feature_renamed))]), collapse = ", "))
  
  # only retain relevant columns from feature
  feature <- select(feature_renamed, ElementChr, ElementStart, ElementEnd, GeneEnsemblID, all_of(feature_score_cols))

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
  merged <- cbind(e2g_candidates[queryHits(ovl)], feature[subjectHits(ovl), feature_score_cols, with = FALSE])
  
  # aggregate feature scores, where multiple features overlapped one e2g_candidates enhancer
  agg_cols <- setdiff(colnames(merged), feature_score_cols)
  merged <- aggregate_features(merged, feature_score_cols = feature_score_cols,
                               agg_cols = agg_cols, agg_fun = agg_fun)
  
  # report how many e2g_candidates pairs overlap >= 1 features
  perc_overlap <- nrow(merged) / nrow(e2g_candidates)
  message("Features overlapping E2G candidate universe: ", round(perc_overlap * 100, digits = 2), "%")
  
  # get pairs from e2g_candidates table that are missing from features 
  missing <- e2g_candidates[setdiff(seq_len(nrow(e2g_candidates)), queryHits(ovl)), ]
  
  # fill in missing value for features in missing pairs
  for (i in feature_score_cols) {
    missing[[i]] <- fill_value
  }
  
  # combine merged and missing pairs to create output
  output <- rbind(merged, missing)
  
  # sort output by cre position
  output <- output[order(ElementChr, ElementStart, ElementEnd, GeneEnsemblID), ]
  
  return(output)
  
}

# function to aggregate multiple overlaps
aggregate_features <- function(merged, feature_score_cols, agg_cols, agg_fun) {
  
  # aggregate feature columns as specified by agg_fun
  agg_list <- mapply(FUN = function(feat_col, agg_func, agg_cols) {
    agg_func <- get(agg_func)
    merged[, setNames(.(agg_func(get(feat_col))), feat_col), by = agg_cols]
  }, feat_col = feature_score_cols, agg_func = agg_fun, MoreArgs = list(agg_cols = agg_cols),
  SIMPLIFY = FALSE)
  
  # merge all the aggregates together to make collapsed data.frame
  output <- as.data.table(Reduce(function(df1, df2) merge(df1, df2, by = agg_cols), agg_list))
  
  return(output)
  
}

## Process features --------------------------------------------------------------------------------

message("Loading data")

parser <- ArgumentParser()

parser$add_argument(
    "-e",
    "--e2g-candidates",
    type = "character",
    default = NULL,
    help = "path to E2G candidate universe",
    metavar = "character"
)
parser$add_argument(
    "-f",
    "--features",
    type = "character",
    default = NULL,
    help = "path to E2G feature table",
    metavar = "character"
)
parser$add_argument(
    "-c",
    "--config",
    type = "character",
    metavar = "character",
    default = NULL,
    help = "path to feature merge configuration table",
)
parser$add_argument(
    "-o",
    "--output",
    type = "character",
    default = NULL,
    help = "name of merged feature data",
    metavar = "character"
)

# Get command line data
args <- parser$parse_args()

# load ABC output
e2g_candidates <- fread(args$e2g_candidates)

# load feature file
features <- fread(args$features, header = TRUE, sep = "\t")

# load feature config file
config <- fread(args$config)
# Check for required columns
stopifnot("feature_col" %in% names(config), "feature_rename" %in% names(config), "aggregate_function" %in% names(config))

# only select relevant columns from e2g_candidates
e2g_candidates <- select(e2g_candidates, ElementChr, ElementStart, ElementEnd, ElementName, GeneSymbol, GeneEnsemblID)

# merge features with e2g_candidates
output <- merge_feature_to_e2g_candidates(e2g_candidates, feature = features, feature_score_cols = config$feature_col,
                               new_feature_names = config$feature_rename,
                               agg_fun = config$aggregate_function, fill_value = NA_real_)

# select columns for final output
output <- select(output, ElementChr, ElementStart, ElementEnd, ElementName, GeneSymbol, GeneEnsemblID, all_of(config$feature_rename))

# write output to file
fwrite(output, file = args$output, sep = "\t")
