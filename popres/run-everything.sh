#!/bin/bash

echo "Running POPRES_SNPdata_recode12.R to recode the .tped file into a numeric matrix."
echo "Creates coded_data_chr1.txt"
Rscript -e "source('POPRES_SNPdata_recode12.R')"

echo "Running POPRES_cov.R to find the covariance matrix for the entire chromosome."
echo "Creates cov_data_for_chr1.txt"
Rscript -e "source('POPRES_cov.R')"

echo "Running POPRES_PCA_win100.R to create a matrix which has one column for each window,"
echo "and whose rows are eigen vectors/values in the order (vec1, value1, vec2, value2)."
echo "Creates fluctuation_PCA_win_100_chr1.txt"
Rscript -e "source('POPRES_PCA_win100.R')"

echo "Running POPRES_distance.R to compute the distance matrix between each pair of windows."
echo "Creates quick_method_pairwise_distance_between_win_100_chr1"
Rscript -e "source('POPRES_distance.R')"

echo "Running POPRES_MDS.R to make MDS plots of the distance matrix in quick_method_pairwise_distance_between_win_100_chr1"
echo "Creates MDS_2D_win100_chr1.pdf and MDS_1D_win100_chr1.pdf"
Rscript -e "source('POPRES_MDS.R')"
