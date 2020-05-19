This repository contains the PIC-seq algorithm, and all the needed code and metadata files to analyze and generate figures for Giladi A. & Cohen M. et al. Nature Biotechnology 2020
In order to run the scripts, download processed data from GSE135382 to the folder output/umi.tab
Unzip and change file names by running these shell commands (in output/umi.tab):
gzip -d *
ls -1 | awk -F'_' '{print $0,$2}' | xargs -n 2 mv

To start analysis, run from the root directory: Rscript scripts/run_pic_seq.r

Please send questions to Amir Giladi: amir.goldberg@weizmann.ac.il

Package versions:
-----------------
tglkmeans_0.2.0 (install.packages('tglkmeans', repos=c(getOption('repos'), 'https://tanaylab.github.io/repo')))
Hmisc_4.2-0
ggplot2_3.1.0
Formula_1.2-3
survival_2.43-3
lattice_0.20-38
glmnet_2.0-16
foreach_1.4.4
Matrix_1.2-18
compositions_1.40-2
bayesm_3.1-0.1
energy_1.7-5
robustbase_0.93-3
tensorA_0.36.1
gplots_3.0.1.1
plotrix_3.7-4
plyr_1.8.4
RANN_2.6.1
reshape2_1.4.3
KernSmooth_2.23-15
metacell_0.3.41
tgstat_2.3.5
misha_4.0.6
gtools_3.8.1
