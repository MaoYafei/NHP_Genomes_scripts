"""
Run subsequence validations. Per reference region, extract the number of bases from the aligned contig or assembly and quantify their length.
"""

import numpy as np
import os, sys
import pandas as pd

from Bio import SeqIO

SDIR = os.path.realpath(os.path.dirname(srcdir("Snakefile")))

sys.path.append(SDIR)


import subseqlib.seq
import subseqlib.stats
import subseqlib.pipeline
import subseqlib.validate


#
# Definitions
#

PWD = os.getcwd()


if config.get("subseq_definitions_path", None):
    suffix = config.get("subseq_definitions_path")
else:
    suffix = "definitions.snakefile"

include: os.path.join(PWD, suffix)




#
# Validation tables
#

# subseq_val
#
# Make a validation table
rule subseq_val:
    input:
        tsv='tables/subseq/{sample}/{source}_{caller}/sv_{svtype}/{alnsample}_{alnsource}.tsv.gz'
    output:
        tsv='tables/validation/{sample}/{source}_{caller}/sv_{svtype}/{alnsample}_{alnsource}/{strategy}.tsv.gz'
    run:
        
        # Check alignsource generator function
        if ALNSOURCE_PLOIDY_DICT[wildcards.alnsource] != subseqlib.stats.align_summary_diploid:
            raise RuntimeError('Validation with size50 requires results from align_summary_diploid')
        
        # Read
        df = pd.read_csv(input.tsv, sep='\t')
        
        # Validate
        df = subseqlib.validate.validate_summary(df, strategy=wildcards.strategy)
        
        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

#
# Subsequence tables
#

# subseq_merge_setdef
#
# Merge subseq results across set definitions (one table for variant callset vs one alignment source).
rule subseq_merge_setdef:
    input:
        tsv=expand('temp/tables/sample/{{sample}}/{{source}}_{{caller}}/{set_def}/sv_{{svtype}}/{{alnsample}}_{{alnsource}}.tsv.gz', set_def=SET_DEF.keys())
    output:
        tsv='tables/subseq/{sample}/{source}_{caller}/sv_{svtype}/{alnsample}_{alnsource}.tsv.gz'
    run:
        return pd.concat(
            [pd.read_csv(tsv_file_name, sep='\t') for tsv_file_name in input.tsv],
            axis=0
        ).to_csv(
            output.tsv, sep='\t', index=False, compression='gzip'
        )

# subseq_tab_window_single
#
# Get stats for a single window.
rule subseq_tab_window_single:
    input:
        bed='tables/subseq/{sample}/{source}_{caller}/{set_def}/setdef_sv_{svtype}.bed.gz',
        aln=lambda wildcards: subseqlib.pipeline.get_aln_source(wildcards, ALNSOURCE_PATTERN_DICT)
    output:
        tsv=temp('temp/tables/sample/{sample}/{source}_{caller}/{set_def}/sv_{svtype}/{alnsample}_{alnsource}.tsv.gz')
    run:
        
        # Read BED
        
        # Read region BED
        region_bed = pd.read_csv(input.bed, sep='\t')
        
        region_bed['ALNSAMPLE'] = wildcards.sample
        region_bed['ALNSOURCE'] = wildcards.alnsource

        region_bed = region_bed.loc[
            :, ['ID', 'SVTYPE', 'SVLEN', 'WINDOW', 'WINDOW_SIZE', 'SAMPLE', 'CALLER', 'ALNSAMPLE', 'ALNSOURCE']
        ].reset_index(drop=True)
        
        # Get alignment source
        aln_input_exists = os.path.isfile(input.aln)
        
        # Get summary function
        summary_func = ALNSOURCE_PLOIDY_DICT[wildcards.alnsource]
        
        # Build a dict of stat records
        
        # Collect stats for each region
        if aln_input_exists:
            # Alignment file exists, get stats
            
            stat_list = [summary_func(subseqlib.seq.get_len_list(window, input.aln, SUBSEQ_EXE)) for window in region_bed['WINDOW']]
                            
        else:
            # Alignment file does not exist, get empty lists
            
            stat_list = [summary_func(list())] * region_bed.shape[0]
        
        # Merge stats
        if len(stat_list) > 0:
            # Merge summary records
            df = pd.DataFrame(pd.concat(stat_list, axis=1)).T.reset_index(drop=True)
            
        else:
            # Create an empty DataFrame with the right headers if there were no variants in the input
            df = pd.DataFrame([], columns=summary_func(list()).index)
        
        df['HAS_ALN'] = aln_input_exists

        # Convert frame types
        dtype_dict = {field: dtype for field, dtype in subseqlib.stats.ALIGN_SUMMARY_FIELD_DTYPE.items() if field in df.columns}

        df = df.astype(dtype_dict)

        # Prepend region information
        df = pd.concat([region_bed, df], axis=1)

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# subseq_table_sv_bed
#
# Make SV validation table.
#
# Wildcards:
# * sample: Biological sample
# * set_def: Entry in SET_DEF (without the sv prefix)
# * caller: Variant caller
# * svtype: ins or del
rule subseq_table_setdef_bed:
    input:
        bed=lambda wildcards: subseqlib.pipeline.get_variant_input(wildcards, VARIANT_BED_PATTERN)
    output:
        bed='tables/subseq/{sample}/{source}_{caller}/{set_def}/setdef_sv_{svtype}.bed.gz'
    run:
        
        # Get parameters
        if wildcards.set_def not in SET_DEF:
            raise RuntimeError('Set definition not found in SET_DEF: ' + wildcards.set_def)

        svlen_min, svlen_max, win_size = SET_DEF[wildcards.set_def]
        
        # Process variants
        if 'bed' in input.keys() and bool(input.bed):  # Function returns an empty list if no input, and Snakemake omits those from the input object
            
            df = pd.read_csv(input.bed, sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN'))

            df_fai = subseqlib.seq.get_ref_fai(REF_FA + '.fai')
            
            # Filter chromosomes that were not aligned to (alt and decoys for short-read callers)
            df = df.loc[df['#CHROM'].apply(lambda val: val in df_fai.index)]

            # Subset
            if svlen_min is not None:
                df = df.loc[df['SVLEN'] >= svlen_min]

            if svlen_max is not None:
                df = df.loc[df['SVLEN'] < svlen_max]

            # Add sample and caller
            df['SAMPLE'] = wildcards.sample
            df['CALLER'] = wildcards.caller

            # Annotate window (replace POS and END with the window coordinates)
            df['WIN_FLANK'] = win_size

            # Define windows, replace POS and END
            df['ORG_POS'] = df['POS']
            df['ORG_END'] = df['END']

            if df.shape[0] > 0:
                df['POS'] = df.apply(lambda row: np.max([row['POS'] - row['WIN_FLANK'], 0]), axis=1)
                df['END'] = df.apply(lambda row: np.min([row['END'] + row['WIN_FLANK'], df_fai[row['#CHROM']]]), axis=1)

            df['WIN_L'] = df['ORG_POS'] - df['POS']
            df['WIN_R'] = df['END'] - df['ORG_END']

            df['WINDOW_SIZE'] = df['WIN_L'] + df['WIN_R']

            del(df['ORG_POS'], df['ORG_END'])

            if df.shape[0] > 0:
                df['WINDOW'] = df.apply(lambda row: '{#CHROM}:{POS}-{END}'.format(**row), axis=1)
            else:
                df['WINDOW'] = []
        else:
            # No variant calls
            
            df = pd.DataFrame([], columns=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SAMPLE', 'CALLER', 'WIN_FLANK', 'WIN_L', 'WIN_R', 'WINDOW_SIZE', 'WINDOW'])

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
