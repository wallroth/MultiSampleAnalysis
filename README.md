# MultiSampleAnalysis
R functions for advanced stratified analyses of multi-sample data

Author: Raphael Wallroth

Maintainer: Raphael Wallroth <raphaelwallroth@gmail.com>

Imports: plyr, utils, data.table

Description: This package provides various tools to preprocess, transform and analyse multi-sample data. Many functions in R and its various packages will simply consider every sample of your data as a new, independent observation. However, in cases where multiple measurements per stimulus event are collected this assumption does not hold true and dependencies of within-trial samples have to be considered. This package was specifically designed for the analysis of electroencephalogram (EEG) data but should be easily applicable to any data with more than one measurement per trial (e.g. eye-tracking, kinematics, MEG, etc.). 

ToDo: extensive documentation to be able to compile a library

for now, use with: source("decode.R") 

Functions are loaded into their own environments so they won't appear in your workspace.