import logging
import re
from pathlib import Path
import pandas as pd
import sys
import gzip
import argparse
from typing import List, Dict, Set, Tuple
import numpy as np
import os
import difflib


# Data Validation
DTYPES = {
    'ElementChr': 'string',
    'ElementStart': 'Int64',  # Use nullable integer to catch NaNs manually
    'ElementEnd': 'Int64',    # Use nullable integer to catch NaNs manually
    'ElementName': 'string',
    'ElementClass': 'string',
    'GeneSymbol': 'string',
    'GeneEnsemblID': 'string',
    'GeneTSS': 'Int64',
    'SampleSummaryShort': 'string',
    'Score': 'float'          # Float is inherently nullable (uses np.nan)
}
REQUIRED_COLS = list(DTYPES.keys())

# Matches chr + (1-22 or X/Y/M)
CHR_REGEX = re.compile(r"^chr((?:[1-9]|1\d|2[0-2])|[XYM])$")

# Matches a valid Ensembl ID format, e.g., 'ENSG00000139618' or 'ENST00000384233.3'
ENSEMBL_ID_REGEX = re.compile(r"^ENS[A-Z]{1,5}\d{11}(?:\.\d+)?$")

# Allowed values for ElementClass (case-sensitive)
ALLOWED_ELEMENT_CLASSES = {'promoter', 'genic', 'intergenic'}

# Metadata Check
REQUIRED_META_KEYS = {
    "Source", "Version", "GenomeReference", "URL", "Assays",
    "SampleAgnostic", "SampleTermName", "SampleTermID", "SampleSummaryShort", "ScoreType"
}
OPTIONAL_META_KEYS = {"ScoreThreshold", "Metadata"}

# Pre-compiled regexes for performance and clarity.
# META_LINE_REGEX = re.compile(r"^#\s*([^:]+):\s*(.*)$")
GENOME_REF_ACCESSION_REGEX = re.compile(r"^IGVF[A-Z]{2}[0-9]{4}[A-Z]{4}$")
# GENOME_REF_ALIAS_REGEX = re.compile(r"^\S+:\S+$") # Matches 'string:string'

META_DESCRIPTIONS: Dict[str, str] = {
    "Source": "the predictive model used",
    "Version": "version of predictive model",
    "GenomeReference": "Accession ID to the IGVF genome reference object (ex: IGVFDS0280IQAI)",
    "URL": "link to code repository, documentation or source of E-G links",
    "Assays": "molecular assay(s) used to define candidate elements",
    "SampleAgnostic": "boolean value 'True' if prediction is applicable to all biosamples, 'False' if not",
    "SampleTermName": "biosample/cell type name (required if not SampleAgnostic)",
    "SampleTermID": "UBERON, CL for cellType or biosample; not needed if SampleAgnostic (https://www.ebi.ac.uk/efo/)",
    "SampleSummaryShort": "brief description of the sample, including treatments",
    "ScoreType": "one of [positive_score, negative_score, p_value, adj_p_value, divergent, boolean]",
    "ScoreThreshold": "(thresholded predictions only) - cutoff or threholding strategy used to select predictions)",
    "Metadata": "IGVF data portal accession ID for full metadata",
}

def _validate_chr_column(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates the 'ElementChr' column format."""
    if 'ElementChr' in reported:
        return [], []
    
    bad_rows = chunk[~chunk['ElementChr'].str.match(CHR_REGEX, na=False)]
    if not bad_rows.empty:
        row = bad_rows.iloc[0]
        line = start_row + row.name
        error = (f"Invalid Format [L{line}]: 'ElementChr' value '{row['ElementChr']}' "
                 "is not in Gencode/UCSC notation.  Possible values: [chr1, chr2, ..., chr22, chrX chrY, chrM]")
        reported.add('ElementChr')
        return [error], []
    return [], []

def _validate_coordinate_columns(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates 'ElementStart' and 'ElementEnd' columns."""
    errors = []
    for col in ['ElementStart', 'ElementEnd']:
        if col not in chunk.columns: continue

        # Check for missing values
        if f'{col}_nan' not in reported and chunk[col].isnull().any():
            row = chunk[chunk[col].isnull()].iloc[0]
            line = start_row + row.name
            errors.append(f"Missing Value [L{line}]: '{col}' coordinate is required and cannot be blank.")
            reported.add(f'{col}_nan')

        # Check for non-positive values
        if f'{col}_value' not in reported and (chunk[col].dropna() < 0).any():
            row = chunk[(chunk[col].notna()) & (chunk[col] < 0)].iloc[0]
            line = start_row + row.name
            errors.append(
                f"Invalid Value [L{line}]: '{col}' coordinate must be a nonnegative integer, but found '{int(row[col])}'."
            )
            reported.add(f'{col}_value')
    return errors, []

def _validate_non_empty_strings(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates that certain string columns are not empty."""
    errors = []
    for col in ['ElementName', 'GeneSymbol', 'SampleSummaryShort']:
        if col not in chunk.columns or col in reported: continue
        
        invalid = chunk[chunk[col].isnull() | (chunk[col].str.strip().str.len() == 0)]
        if not invalid.empty:
            row = invalid.iloc[0]
            line = start_row + row.name
            errors.append(f"Invalid Value [L{line}]: Column '{col}' cannot be blank or empty.")
            reported.add(col)
    return errors, []

def _validate_element_class(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates 'ElementClass' against a list of allowed values (Warning only)."""
    if 'ElementClass' not in chunk.columns or 'ElementClass' in reported:
        return [], []
        
    invalid = chunk[~chunk['ElementClass'].str.lower().isin(ALLOWED_ELEMENT_CLASSES) & chunk['ElementClass'].notna()]
    if not invalid.empty:
        row = invalid.iloc[0]
        line = start_row + row.name
        warning = (f"Data Warning [L{line}]: 'ElementClass' has an unrecognized value "
                   f"'{row['ElementClass']}'. Allowed values are {ALLOWED_ELEMENT_CLASSES}.")
        reported.add('ElementClass')
        return [], [warning]
    return [], []

def _validate_ensembl_id(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates the 'GeneEnsemblID' format."""
    if 'GeneEnsemblID' not in chunk.columns or 'GeneEnsemblID' in reported:
        return [], []

    invalid = chunk['GeneEnsemblID'].dropna()[~chunk['GeneEnsemblID'].dropna().str.match(ENSEMBL_ID_REGEX)]
    if not invalid.empty:
        row = invalid.to_frame().iloc[0]
        line = start_row + row.name
        error = (f"Invalid Format [L{line}]: 'GeneEnsemblID' value '{row['GeneEnsemblID']}' "
                 "is not a valid Ensembl ID format (e.g., ENSG00000136997).")
        reported.add('GeneEnsemblID')
        return [error], []
    return [], []

def _validate_nonnegative_TSS(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates 'GeneTSS' is a non-negative integer"""
    if 'GeneTSS' not in chunk.columns or 'GeneTSS_negative' in reported:
        return [], []
    
    invalid = chunk[chunk['GeneTSS'].dropna() < 0]
    if not invalid.empty:
        row = invalid.iloc[0]
        line_num = start_row + row.name
        error = f"Invalid Value [L{line_num}]: 'GeneTSS' must be a nonnegative coordinate, not {int(row['GeneTSS'])}."
        reported.add('GeneTSS_negative')
        return [error], []
    return [], []

def _validate_not_blank_TSS(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates 'GeneTSS' conatins no blanks (Warning only)"""
    if 'GeneTSS' not in chunk.columns or 'GeneTSS_blank_warn' in reported or 'GeneTSS_negative' in reported:
        return [], []
    
    invalid = chunk[chunk['GeneTSS'].isnull()]
    if not invalid.empty:
        row = invalid.iloc[0]
        line_num = start_row + row.name
        warning = f"Data Warning [L{line_num}]: 'GeneTSS' should not contain blank/NaN values."
        reported.add('GeneTSS_blank_warn')
        return [], [warning]
    return [], []

def _validate_score(chunk: pd.DataFrame, start_row: int, reported: Set[str]) -> Tuple[List[str], List[str]]:
    """Validates 'Score' is numeric and contains no NaN."""
    if 'Score' not in chunk.columns or 'Score' in reported:
        return [], []

    invalid_rows = chunk[chunk['Score'].isnull()]
    if not invalid_rows.empty:
        first_bad_row = invalid_rows.iloc[0]
        line_num = start_row + first_bad_row.name
        error = f"Invalid Value [L{line_num}]: 'Score' column cannot contain blank or NaN values."
        reported.add('Score')
        return [error], []
    return [], []


def _check_data_rows(file_path: Path, parsed_meta: Dict[str, str], check_all_rows: bool = False) -> Tuple[List[str], List[str]]:
    """
    Orchestrates data row validation, providing specific feedback on missing columns
    and validating data only in the columns that are correctly named.
    """
    data_errors: List[str] = []
    data_warnings: List[str] = []

    # --- 1. Header Analysis ---
    try:
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            for i, line in enumerate(f, 1):
                if not line.strip().startswith('#'):
                    actual_columns = [col.strip() for col in line.split('\t')]
                    header_line_num = i
                    break
    except Exception as e:
        data_errors.append(f"Fatal Parsing Error: Could not read file header. Error: {e}")
        return data_errors, data_warnings

    expected_set = set(REQUIRED_COLS)
    actual_set = set(actual_columns)
    
    missing_columns = expected_set - actual_set
    found_columns = expected_set & actual_set
    extra_columns = actual_set - expected_set # For typo suggestions

    # --- 2. Report Missing Columns ---
    if missing_columns:
        for col in sorted(list(missing_columns)):
            suggestion = ""
            close_matches = difflib.get_close_matches(col, list(extra_columns), n=1, cutoff=0.7)
            if close_matches:
                suggestion = f" Closest matching column found: '{close_matches[0]}'."
            data_errors.append(f"Missing Column: The required column '{col}' was not found.{suggestion}")

    # If no valid columns were found, we can't proceed.
    if not found_columns:
        return data_errors, data_warnings

    # --- 3. Data Validation on Found Columns ---
    reported_error_categories: Set[str] = set()
    chunksize = 100000
    
    # Filter dtypes to only those columns we are actually loading
    dtypes_to_use = {k: v for k, v in DTYPES.items() if k in found_columns}

    try:
        chunk_iterator = pd.read_csv(
            file_path,
            sep='\t',
            comment='#',
            compression='gzip',
            skiprows=header_line_num, # Skip metadata AND the header we already read
            header=None, # We've already processed the header
            names=actual_columns, # Use actual names from file
            usecols=list(found_columns), # IMPORTANT: Only load correctly named columns
            dtype=dtypes_to_use,
            chunksize=chunksize,
            low_memory=False
        )

        for i, chunk in enumerate(chunk_iterator):
            chunk_start_row = header_line_num + (i * chunksize) + 1

            # Call each validation helper
            errors, warnings = _validate_chr_column(chunk, chunk_start_row, reported_error_categories)
            data_errors.extend(errors)
            data_warnings.extend(warnings)

            errors, _ = _validate_coordinate_columns(chunk, chunk_start_row, reported_error_categories)
            data_errors.extend(errors)
            # This helper doesn't produce warnings

            errors, _ = _validate_non_empty_strings(chunk, chunk_start_row, reported_error_categories)
            data_errors.extend(errors)
            # This helper doesn't produce warnings

            _, warnings = _validate_element_class(chunk, chunk_start_row, reported_error_categories)
            # This helper doesn't produce errors
            data_warnings.extend(warnings)
            
            errors, _ = _validate_ensembl_id(chunk, chunk_start_row, reported_error_categories)
            data_errors.extend(errors)
            # This helper doesn't produce warnings
            
            errors, _ = _validate_nonnegative_TSS(chunk, chunk_start_row, reported_error_categories)
            data_errors.extend(errors)
            # This helper doesn't produce warnings

            _, warnings = _validate_not_blank_TSS(chunk, chunk_start_row, reported_error_categories)
            # This helper doesn't produce errors
            data_warnings.extend(warnings)

            errors, _ = _validate_score(chunk, chunk_start_row, reported_error_categories)
            data_errors.extend(errors)
            # This helper doesn't produce warnings

            if not check_all_rows and (data_errors or data_warnings):
                break

    except (ValueError, TypeError) as e:
        msg = (f"Fatal Parsing Error: Data in a row could not be converted to its expected type "
               f"(e.g., text in a number column). Please check data format. Pandas error: {e}")
        data_errors.append(msg)
    except Exception as e:
        data_errors.append(f"An unexpected error occurred during data processing: {e}")

    return data_errors, data_warnings


def _check_metadata_header(file_path: Path) -> Tuple[List[str], Dict[str, str]]:
    """
    Reads and parses metadata from a gzipped file, returning advisories.

    This function reads comment lines from the top of a file, parses them as
    key-value metadata, and checks for common issues. It does not cause a
    hard failure but returns a list of advisory warnings.

    Args:
        file_path: The path to the gzipped prediction file.

    Returns:
        A tuple containing:
        - A list of string warnings about the metadata.
        - A dictionary of the parsed metadata.
    """
    warnings: List[str] = []
    found_meta: Dict[str, str] = {}
    key_locations: Dict[str, int] = {}

    try:
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                stripped = line.strip()

                if not stripped:
                    continue  # Skip blank lines

                if stripped.startswith('#'):
                    # Strip leading '#' and try to split into key-value
                    parts = stripped.lstrip('#').split(':', 1)
                    if len(parts) == 2:
                        key, value = parts
                        key = key.strip()
                        
                        # Rule: If the key is not known, ignore this line as a simple comment.
                        if key not in META_DESCRIPTIONS:
                            continue

                        value = value.strip()
                        if key in found_meta:
                            # Rule: Duplicate keys are a simple warning, not a critical error.
                            first_occurrence = key_locations.get(key, 'N/A')
                            warnings.append(
                                f"Duplicate metadata key '{key}' on line {line_num}. "
                                f"The first value from line {first_occurrence} will be used."
                            )
                            continue  # Ignore this duplicate key, keeping the first one
                        
                        found_meta[key] = value
                        key_locations[key] = line_num
                    # else: it's a comment without a ':', so we ignore it.
                else:
                    # We've reached the first non-comment, non-blank line. Metadata is over.
                    break
    except (gzip.BadGzipFile, EOFError) as e:
        msg = f"ERROR: File '{file_path.name}' is not a valid gzip file or is empty. Cannot check metadata."
        logging.warning(msg) # Downgraded from error to warning
        warnings.append(msg)
        return warnings, {}
    except Exception as e:
        msg = f"ERROR: Could not read metadata from '{file_path.name}': {e}"
        logging.warning(msg) # Downgraded from error to warning
        warnings.append(msg)
        return warnings, {}

    # --- Post-parsing validation logic ---

    is_sample_agnostic = found_meta.get("SampleAgnostic", "").lower() == 'true'

    # Check for missing REQUIRED keys
    missing_keys = REQUIRED_META_KEYS - set(found_meta.keys())
    for key in sorted(list(missing_keys)):
        if key == "SampleTermName" and is_sample_agnostic:
            continue
        description = META_DESCRIPTIONS.get(key, "No description available.")
        warnings.append(f"Missing required field '{key}'. This field should contain: {description}")

    # Check for issues in existing keys
    for key, value in found_meta.items():
        if not value:
            if key == "SampleTermName" and is_sample_agnostic:
                continue
            if key in META_DESCRIPTIONS:
                description = META_DESCRIPTIONS.get(key, "No description available.")
                warnings.append(f"Field '{key}' is present but empty. It should be filled with: {description}")
            continue

        if key == "GenomeReference" and not GENOME_REF_ACCESSION_REGEX.match(value):
            warnings.append(f"Field 'GenomeReference' has a malformed value '{value}'. It should be a proper accession ID like 'IGVFFI0000GXML'.")

        if key == "SampleAgnostic" and value.lower() not in ['true', 'false']:
            warnings.append(f"Field 'SampleAgnostic' must be 'True' or 'False' (case-insensitive), but found '{value}'.")

    return warnings, found_meta

# --- Validator Function ---

def validate_prediction_file(file_path: Path) -> Tuple[bool, str]:
    """
    Validates a prediction file against standard format rules.

    This function orchestrates metadata and data validation, then compiles
    all feedback into a single, user-friendly string.

    Args:
        file_path: The path to the prediction file to validate.

    Returns:
        A tuple containing:
        - A boolean, True if the file has no data errors, False otherwise.
        - A formatted string detailing all errors and warnings found.
    """
    feedback_parts = []

    # Step 1: Check metadata header
    metadata_warnings, parsed_meta = _check_metadata_header(file_path)

    # Step 2: Check data rows
    data_errors, data_warnings = _check_data_rows(file_path, parsed_meta)

    # Step 3: Compile feedback into the desired format
    if data_errors:
        feedback_parts.append("\nData Errors:")
        for error in data_errors:
            feedback_parts.append(f">> {error}")

    if data_warnings:
        feedback_parts.append("\nData Warnings:")
        for warning in data_warnings:
            feedback_parts.append(f">> {warning}")
            
    if metadata_warnings:
        feedback_parts.append("\nMetadata Warnings:")
        for warning in metadata_warnings:
            feedback_parts.append(f">> {warning}")

    # Step 4: Determine final status and add summary
    overall_is_valid = (len(data_errors) == 0)
    
    if overall_is_valid and not (data_warnings or metadata_warnings):
        # Case: Perfect file
        feedback_parts.append(f"\nFILE IS CORRECTLY FORMATED: {file_path.name}")
        feedback_parts.append("Great job!")
    elif overall_is_valid:
        # Case: Valid, but with warnings
        feedback_parts.append(f"\nFILE IS IN IGVF STANDARD FORMAT: {file_path.name}")
        feedback_parts.append("Note: Please review the warnings above to improve file quality.")
    else:
        # Case: Invalid due to data errors
        feedback_parts.append(f"\nFILE IS *NOT* IN IGVF STANDARD FORMAT: {file_path.name}")
        feedback_parts.append("Please address the data errors to ensure your E2G predictions are in standard format.")

    return overall_is_valid, "\n".join(feedback_parts)


def is_valid_readable_file(path_str: str) -> Path:
    """
    Custom argparse type for a file path that must exist and be readable.
    Converts the string path to a pathlib.Path object.
    """
    path = Path(path_str)

    # 1. Check if the path exists
    if not path.exists():
        raise argparse.ArgumentTypeError(f"The file '{path}' does not exist.")

    # 2. Check if it's a file (not a directory)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"The path '{path}' is a directory, not a file.")

    # 3. Check for read permissions
    if not os.access(path, os.R_OK):
        raise argparse.ArgumentTypeError(f"The file '{path}' is not readable.")

    return path

# --- Command-Line Interface ---

def main(prediction_file: Path):
    """A command-line tool to validate IGVF E2G prediction file formats."""

    # Run the validation
    print(f"Checking formatting of {prediction_file}...")
    is_valid, feedback_message = validate_prediction_file(prediction_file)

    # Print the comprehensive feedback
    print(feedback_message)
    
    # Exit with a status code indicating validity
    if not is_valid:
        sys.exit(1)
    
    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check E2G prediction file against IGVF standard format guidelines"
    )
    parser.add_argument(
        '-f', '--prediction-file',
        dest='prediction_file',
        # Use custom function to validate and convert the input
        type=is_valid_readable_file,
        required=True,
        help='Path to the E2G prediction file to be validated.'
    )
    args = parser.parse_args()
    main(args.prediction_file)