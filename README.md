# Cardio project
**This is a basic instruction for those working with plink program** 
## First part
We had initial vcf file with genotypes frequencies from patients with pulmonary arterial hypertension(PAH). We wanted to compare them with genotypes frequencies of healthy donors from 1000G project. For this purpose we chose populations of Non-Finnish Europeans and Americans from Utah State, which we believed have close ethnical background to Slavic population, which most our patients belonged to.  
1. We recoded genotypes with reference/alternative allele to 0/1 format accordingly
2. Converted file vcf format to ped 
3. recoded gnomad .fam file as controls
4. Merged case(PAH) and controls (NFE and Utah)
5. Performed PCA
6. Performed ascocciation study based on chi-square analysis
## Second part
1. Based on the result of ascocciation study we performed manhattan plot via qqman package in R. However the originak font size was not suitable for us, thus I changed it manualy in the function.
2. Based on the resulted PCA from plink we build the plot with R
