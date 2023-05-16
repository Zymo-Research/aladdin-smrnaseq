#!/usr/bin/env python
"""
Setup code for Zymo Research MultiQC plugin.
For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '1.0.0'

setup(
    name = 'multiqc_zymo',
    version = version,
    description = "MultiQC plugins for smrnaseq pipeline",
    packages = find_packages(),
    include_package_data = True,
    install_requires = ['multiqc==1.9'],
    entry_points = {
        'multiqc.templates.v1': [
            'zymo = multiqc_zymo.templates.zymo'
        ],
        'multiqc.modules.v1': [
            'plot_sample_distance = multiqc_zymo.modules.plot_sample_distance:MultiqcModule',
            'plot_gene_heatmap = multiqc_zymo.modules.plot_gene_heatmap:MultiqcModule',
            'isomiRs = multiqc_zymo.modules.isomiRs:MultiqcModule',
            'deseq2_rnatypes = multiqc_zymo.modules.deseq2_rnatypes:MultiqcModule',
            'Trim_Galore = multiqc_zymo.modules.Trim_Galore:MultiqcModule'
        ],
        'multiqc.hooks.v1': [
            'before_config = multiqc_zymo.custom_code:plugin_before_config',
            'execution_start = multiqc_zymo.custom_code:plugin_execution_start'
        ]
    }
)