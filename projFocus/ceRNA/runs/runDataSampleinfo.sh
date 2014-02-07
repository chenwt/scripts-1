#!/bin/bash
#By: J.He
#TODO: 

awk '{print substr($2,0,19)}' brca_meth_level3_sampleNameSummary_11102013.txt >brca_meth_level3_sampleNameSummary_11102013.txt_sample
grep -wf brca_wgs_bam_summary_02042014.tsv_TumorSample brca_meth_level3_sampleNameSummary_11102013.txt_sample > overlap_wgs_meth_sample

