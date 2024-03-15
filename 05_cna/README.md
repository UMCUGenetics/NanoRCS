# 05_cna

ichorCNA is a widely adopted tool for cell-free DNA copy number alteration detection and tumor fraction estimation tool. 
Here you can find `ichorCNA` cloned as a submodule with small adaptation in the snakemake pipeline.

The adaptations are listed below:
- Allow absolute output path location instead of relative path originally
- Include conda environment automatically install by snakemake pipeline in envs folder
- Resourses udpate 

ichorCNA is included as a submodule in this repository. The code difference can be found with git diff to the original 
commit: `git diff 5bfc03e`.

We derived CNA with two mode, when derived TF is smaller than 0.03, we perform MRD mode which does not infer 
subclonal status. For NanoRCS CNA inference, we use customized PoN. For NovaSeq we use PoN provided by the authors.

- (1) high-ploidy and high-tumor fraction mode config file: `config_high_CNA_high_ploidy_cusPoN.yaml`

- (2) MRD mode: `config_MRD_cusPoN.yaml`

More documentation of ichorCNA please be referred to: https://github.com/broadinstitute/ichorCNA

## Downstream analysis: IchorCNA provide multiple solutions especially with allowing high-ploidy solutions 
and subclonal solutions. Sometimes manual curation is important to get the most likely true copy number profiles. 
For samples with known CNA and ploidy in the tumor tissue, we manually select the highest scoring solution with fitting 
ploidy. Exporting the output of the selected solution is done via R studio.  

### Output: 

Downstream plotting scripts are located in `10_figures`. Scripts refer to ichorCNA output path relative 
`/path/to/NanoRCS/output/processed_data/05_cna/ichorCNA_curated_solution/`.