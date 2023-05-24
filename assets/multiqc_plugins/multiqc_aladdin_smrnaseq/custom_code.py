#!/usr/bin/env python
"""
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_aladdin_smrnaseq_version = get_distribution("multiqc_aladdin_smrnaseq").version

# Add default config options that can be overriden by user config
def plugin_before_config():
    
    # Use the zymo template by default
    config.template = 'aladdin'
    
# Add additional config options
def plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsed.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', False) is True:
        return None

    log.info("Running smrnaseq MultiQC Plugins v{}".format(config.multiqc_aladdin_smrnaseq_version))

    # Add to the search patterns used by modules
    if 'plot_sample_distance/heatmap' not in config.sp:
        config.update_dict( config.sp, { 'plot_sample_distance/heatmap': [ {'fn' : '*sample_distance_matrix*'}, {'fn' : '*sample_similarity_matrix*'} ] } )
    if 'plot_sample_distance/pca' not in config.sp:
        config.update_dict( config.sp, { 'plot_sample_distance/pca': [ {'fn' : '*sample_pca_plot.tsv'}, {'fn' : '*sample_PCA_plot.tsv'}, {'fn' : '*sample_MDS_plot.tsv'}, {'fn' : '*sample_mds_plot.tsv'} ] } )
    if 'plot_gene_heatmap' not in config.sp:
        config.update_dict( config.sp, { 'plot_gene_heatmap': [ { 'fn' : '*gene_heatmap.tsv' }, { 'fn' : '*transcript_heatmap.tsv' } ] } )
    if 'download_data' not in config.sp:
        config.update_dict( config.sp, { 'download_data' : { 'fn' : '*download_links.json', 'shared': True } } )
    if 'isomiRs' not in config.sp:
        config.update_dict( config.sp, { 'isomiRs' : { 'fn' : '*isomiRs_results.tsv' } } ) 
    if 'deseq2_rnatypes/data' not in config.sp:
        config.update_dict( config.sp, {'deseq2_rnatypes/data' : { 'fn' : '*deseq2_results.tsv'} } )
    if 'Trim_Galore/logs' not in config.sp:
        config.update_dict( config.sp, { 'Trim_Galore/logs' : { 'contents': 'This is cutadapt', 'shared': True } } )
    if 'Trim_Galore/trimmedfastq' not in config.sp:
        config.update_dict( config.sp, {'Trim_Galore/trimmedfastq' : { 'fn' : '*trimmed_fastqc.zip', 'shared' : True } } )
    
    # Some additional filename cleaning
    config.fn_clean_exts.extend([
        '_R1',
        '_R2',
        '_isomiRs_results',
        '_deseq2_results'
    ])
