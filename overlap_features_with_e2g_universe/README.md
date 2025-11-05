# Overlap feature(s) into a different E2G candidate universe
Adds specified feature columns from the feature table to the E2G candidate universe.  

Element-gene pairs are overlapped across the E2G universe and your feature table.
* Where there is a 1:1 match, the value from the feature column will be added as data for the E2G pair from the candidate universe.
* Where multiple E2G candidates in the feature table overlap one E2G candidate in the E2G universe, the config-specified aggregation function will be used to aggregate the duplicate feature values.
* Where one E2G candidate from the E2G universe overlaps multiple E2G candidates from the feature table, teh feature value will be copied across matching rows.

All columns from the E2G canddiate universe will the maintained in the final output.  If you have feature columns in the E2G candidate universe file you will create a feature table with features from multiple sources.  You can also use this script as a fancy way to left join multiple features that were calculated in the same E2G candidate universe.

The script will report the percentage of E2G candidates in the E2G universe which were overlapped by an E2G candidate in the feature table.  This is the percentage of rows in the output that will have new feature data.

## How to use:

### 1: Chose an E2G candidate universe
Your E2G candidate universe must have the columns ElementChr, ElementStart, ElementEnd, and GeneEnsemblID and be in tsv.gz format.  Any additional columns will be retained in the output.  

It is recommended to use an [E2G candidate universe from Synapse](https://www.synapse.org/Synapse:syn68258924).

### 2. Select a feature file to merge
For the purpose of overlapping, your feature table must have the columns ElementChr, ElementStart, ElementEnd, and GeneEnsemblID and be in tsv.gz format.  Specify the columns you want to add to the E2G candidate universe in the configuration table (see below).  **If your feature(s) of interest have the same name as a column in the E2G pair universe, you must specify a new name for those feature(s) in the configuration table.

### 3. Output
Name of the output file (ending in tsv.gz)

### 4. Configuration Table
Example:
| feature_col | feature_rename | aggregation_function |
|---|---|---|
| Score | pgBoost_Score | mean |
| Percentile | | mean |

**feature_col** is where you add the names of columns of interest in the feature table

**feature_rename** is the what you want the feature to be called in your output table

**aggregation_function** is the name of an R function which can be used to when multiple values of the feature map to the same E2G pair in teh E2G universe.

### 5. Run the script

Activate the eg_classifier_features virtual environment.
```bash
mamba env create -f envs/eg_classifier_features.yml
mamba activate eg_classifier_features
```

Overlap your E2G pairs:
```bash
Rscript merge_features_with_e2g_universe.R \
--e2g-candidates path/to/e2g/candidate/universe.tsv.gz \
--features path/to/feature/table.tsv.gz \
--output path/to/output.tsv.gz \
--config path/to/feature_config.tsv
```

Or for short:
```bash
Rscript merge_features_with_e2g_universe.R \
-e path/to/e2g/candidate/universe.tsv.gz \
-f path/to/feature/table.tsv.gz \
-o path/to/output.tsv.gz \
-c path/to/feature_config.tsv
```





