# POCKET
PriOritizing the Candidate genes by incorporating  information of Knowledge-based gene sets, Effects of variants, GWAS and TWAS

<p style="text-align:justify"> Multi-omics datasets were used to prioritize the candidate gene in GWAS QTL regions, and a comprehensive scoring system was established. First, based on variation effect annotated with SNP effector (Cingolani et al., 2012) and GWAS P value of the variation, we used a Min-Max scaled score to evaluate the effects of variation in gene region. Second, based on TWAS P value and cis-eQTL results, we evaluated the gene expression effect in candidate region. Third, we used variations in gene region and upstream region of gene to categorize the gene into different haplotypes and calculated the haplotype-based association score. Fourth, we predicted the potential probability of whether the genes related to the phenotype or not, depends on 8,283 features which we collected from four datasets: (i) GO category, (ii) InterPro protein classification, (iii) gene expression datasets from Lu et al. (Lu et al., 2018), (iv) DEGs of known SOC related mutant or OE lines, (v) ICA ccomponents identified from population transcriptome in 20 and 40 DAF, then we used SVM to predict the gene function. Lastly, we summarized the scores from four processes to determine which genes were more likely effects the phenotype.</p>

Note: we tested the POCKET in Rapeseed and Rice. And we assume that POCKET can be used for other species.

# Requirements
We tested the code on linux platform. 
Requirements are:

GCC C and C++ compiler (gcc, g++)

Python >= 3.6

pip for installing python packages

Numpy

Scipy

scikit-learn

limix

pandas

Joblib

<a href = 'https://www.cog-genomics.org/plink2/'> plink software</a>

# Installation

If all the requirements are met, you can install the POCKET library with the command:

python3 setup.py build install --user

or using

pip install gene_pocket

# Examples

You can follow <a href='./examples/POCKET_example.ipynb'> the example </a> to use the POCKET.

# Citation
Tang, S., Zhao, H., Lu, S., Yu, L., Zhang, G., Zhang, Y., . . . Guo, L. Genome- and transcriptome-wide association studies provide insights into the genetic basis of natural variation of seed oil content in Brassica napus. Molecular Plant. doi:10.1016/j.molp.2020.12.003
