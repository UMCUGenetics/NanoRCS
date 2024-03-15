# After vcf2bed.sh, perform combining of HC variants
bedops --merge HC01-N_snv_only_coor.bed HC02-N_snv_only_coor.bed HC03-N_snv_only_coor.bed > HC_snv.bed