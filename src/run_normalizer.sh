#! /bin/bash
EMAIL=$1
# You will probably need to install BiocManager, preprocessCore, and rhdf5 manually (I've had no real success in
# automating this process) before attempting to  run this script. Run MODULE_CMD before installing any of the R
# packages - some of them depend on these modules.
#
# Be aware that paths to zlib in the cluster are not the same as when you run on a workstation. I found that I needed
# to modify the configure and Makevars.in for Rhdf5lib to ensure that the zlib path was read. Set ZLIB_HOME to the path
# to zlib's "lib" directory in configure, and make sure that the Makevars.in also references this path when building the 
# final package. You'll need to do this for Rhdf5lib and rhdf5 (only in Makevars).
# preprocessCore should install easily enough, once the gcc/6.5.0 module is loaded.

# Make sure R is loaded.
# HDF needs zlib
# gfortran is needed to compile preprocessCore and it is in gcc module
# let's load the HDF module as well.
MODULE_CMD="module load R/3.6.0 && module load zlib/1.2.8 && module load gcc/6.5.0 && module load hdf/1.8.14"
R_CMD="$MODULE_CMD && Rscript normalize_expression_data_hpc.R ./human_matrix.rda ./normalized_human_matrix.hdf"
# Submit the job.
# If you want to try using less than 185GB RAM then you don't need to use the ultra_high_mem queue.
qsub -q ultra_high_mem -o out.log -e err.log -l h_vmem=200G -N normalize-expression-values -cwd -b y -M $EMAIL "$R_CMD"
