# COMS4761

A ReadMe file that tells exactly what is each file, how to compile (if relevant) and test-run the project on sample inputs to get the sample output, 
how to really run the project on large files, if relevant, what are the parameters, what are the system requirements (e.g. are you using a particular
version of MatLab/Python/R/Ruby/Java/C++/Cobol? Do you run on a particularly powerful machine/cloud instance? Which standard or add-on libraries would
you need to have been previously installed?)

Files:
The baseline_v1.2 annotations can be downloaded from https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/GRCh38/baseline_v1.2.tgz.

Data for clumping using plink can be found at https://www.cog-genomics.org/plink/2.0/resources. These are large data files. As such, once I downloaded them, I ran the file get_eur.sh to separate out only the European data, which is a much smaller group and thus speeds up computation.

GWAS and PRS data are saved in the diseases folder, which is too big for GitHub but hopefully not for Courseworks!

get_eur.sh - This file keeps only European IDs for clumping.

clumping.sh - You need to change the author when running this on a new dataset. This file removes duplicate rsIDs, clumps and superclumps the data, saving them into keep.ss.{chromosome #} and superclump.{author}.ss.{chromosome #}.

annotate.R - This file takes in the converted annotations (GrCh38 --> GrCh37) and the GWAS clumped data generated in clumping.sh, appending annotations to the clumped GWAS data.

score_and_cluster.py - You need to change the author name and PGS file when running this on a new dataset. This file runs PCA and GMM to calculate the risky and nonrisky SNPs as identified by the GMM. This file calculates functional annotation scores for all three methods. This file outputs text files saving all of the data.

plotting_diseases.R - This file generates plots - you need to change the autohr name to plot for different diseases.

pca_gmm.py, pca_PGS.py, pca_super.py - These files perform PCA to group the 8 diseases. These take in a file for each disease and return the PCA data.

plotting_PCA.R - I ran this three times, changing the file that is being read, depending on the method. The output of this file is a PCA plot.

super_comparisons.R, PGS_comparisons.R, gmm_comparisons.R - These files generate plots showing the difference in functional annotation for all 8 disesases. 

heatmap.R - This script plots a dendrogram for the 8 different diseases. You need to change the input files to generate a dendrogram for the different methods. 

Additional PRS and GWAS data can be downloaded from PGS Catalog and GWAS Catalog.


Test Input:
To simplify the process, you can download the file keepss, which has already performed the first two steps (the get_eur and clumping). Additionally, you can download converted_annot.txt, which has already run annotate.R. All of this uses the Mahajan, type 2 diabetes data. You can then run score_and_cluster.py. To plot the comparisons, you will need to download the other files that end with .to_plot.gz. You can then run plotting_diseases.R, pca_gmm.py, pca_PGS.py, pca_super.py, plotting_PCA.R, super_comparisons.R, PGS_comparisons.R, gmm_comparisons.R, and heatmap.R.


Outputs
The outputs of these files will be a series of plots and CSV files.
