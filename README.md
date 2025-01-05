# TDbasedUFEHiC
TDbasedUFEHiC: R function that can select optimal bin sizes in Hi-C data

Novel AI-powered computational method using tensor decomposition can discover the common optimal bin sizes when integrating multiple Hi-C datasets
Y-h Taguchi, Turki Turki
bioRxiv 2024.09.29.615651; doi: https://doi.org/10.1101/2024.09.29.615651

The following is how to execute  TDbasedUFEHiC

% wget --user-agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"  https://pmc.ncbi.nlm.nih.gov/articles/instance/7410828/bin/41586_2020_2493_MOESM12_ESM.txt 

% wget https://raw.githubusercontent.com/tagtag/TDbasedUFEHiC/main/sample2_1000000.matrix.gz

% wget https://raw.githubusercontent.com/tagtag/TDbasedUFEHiC/main/sample2_1000000_abs.bed

% gunzip sample2_1000000.matrix.gz

% wget https://raw.githubusercontent.com/tagtag/TDbasedUFEHiC/main/TDbasedUFE_HiC.R

% R

R> source("TDbasedUFE_HiC.R")

R> result <-  TDbasedUFE_HiC("sample2_1000000.matrix","sample2_1000000_abs.bed",sel=T)

The following is the explanaiton of contents of result obtained in the above.
result has three components: COR, SVD, MEANSD

COR includes correlation abalysis shown in Tables 2 to 4 in the paper (Peason and Spearman).

COR$CTCF -> Table 2

COR$PLS -> Table 3

COR$pELS -> Table 4

COR$dELS -> Table 5

However, it is upto the first singular value vectors. If you need results for more singular value vectors, change l_list in the argument of TDbasedUFE_HiC.

SVD includes SVD result. TDbasedUFEHiC generates various pdf files

hist_*.pdf: This is the historgam of 1-P value attributed to regions. Usually, user do not mind it directly.

plot_*.pdf: This is the graph to compute optimal sigma_l. Usually, user do not mind it directly.

image_*.pdf: Cluster structure of the selected region. It corresponds to upper region of Fig. 7 and Fig.8. Please note that in the present case we consider only one sample (profile). Thus the apearance differs.

SVD_*.pdf: Visual representation of selected region. It corresponds to the lower region of Fig. 7

index_*.ccv: The list of selected regions. Acutal genomic loci  must be retrived from  sample2_1000000_abs.bed.

P_*: R object that includes P-values attributed to regions.

*: corresponds to the contents of l_list.

That's all.


