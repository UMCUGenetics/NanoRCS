# 06_fragment_length_NMF_TF

Non-negative matrix factorization is a non-supervised method where a matrix V is factorized into two non-negative 
matrices W and H, both matrices are smaller than the original matrix V. One of these matrices H, the signature matrix, 
has as many columns as the original matrix and represents the preference of observing each fragment length for each 
cfDNA source. The other matrix, the weight matrix W, has as many rows as the input matrix and represents the 
contributions of each cfDNA source to each sample. The number of cfDNA sources is a hyperparameter that needs to be 
set in advance. Renaud et al., 2022, utilized NMF for determining the contribution of different cfDNA 
sources to fragment length signatures in 86 prostate cancer samples. We adapted the 2 signatures extracted f
rom this analysis, and selected the region between 30 to 220 bp, and used them as signatures to decompose the cfDNA 
source in our sample sets. The fragments between 30-220 bp are selected and normalized to 1. 

NMF with fixed signatures with function non_negative_factorization from sklearn.decomposition (v1.1.2) 
is applied to obtain cfDNA source contribution for each sample. All values are capped at 1.0. Implementation 
for calculating fragment lengths and applying NMF are depicted in:
- `10_figures/Fig4A_Fig5C_SupplFig6A_Suppl10C.py`
- `10_figures/Fig4BCD_SupplFig5_SupplFig6C.py` 

For comparison with AFM imaging, analysis script of converting AFM measurement from nm to bp length and compared 
to sequencing results is depicted in:
- `10_figures/Fig4F_SupplFig7L.py`

### Input:
- `data/Renaud_et_al_Fig2B_WGS_Sigs_data.tsv`
- `data/Fig4_AFM_length_nm.csv`
- output from `01_preprocess_nanorcs`: `/path/to/NanoRCS/output/processed_data/01_preprocess_nanorcs/consensus_filtered/stats/`
- - output from `02_preprocess_novaseq`: `/path/to/NanoRCS/output/processed_output/02_preprocessing_novaseq/results/stats_ecco/`