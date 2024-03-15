#!/usr/bin/python
import argparse
import subprocess
import os

def parse_arguments():
	parser = argparse.ArgumentParser(description="Process paths and list of names.")
	parser.add_argument("--bedtools_path", type=str, help="Path to bedtools.")
	parser.add_argument("--temp_path", type=str, help="Path to temporary directory.")
	parser.add_argument("--output_path", type=str, help="Path to output directory.")
	parser.add_argument("--mutect2_prefix", type=str, help="Prefix for Mutect2 files.")
	parser.add_argument("--strelka_prefix", type=str, help="Prefix for Strelka files.")
	parser.add_argument("--haplotype_caller_prefix", type=str, help="Prefix for Haplotype Caller.")
	parser.add_argument("--samples_names", nargs='+', type=str, help="List of sample names.")
	parser.add_argument("--samples_suffix_tumor", type=str, default= '-T', help="Suffix for each tumor sample by sarek.")
	parser.add_argument("--samples_suffix_normal", type=str, default= '-N', help="Suffix for each normal sample by sarek.")

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_arguments()
	# Defile output
	outpath = args.output_path
	outfile = '{}counts.txt'.format(args.output_path)
	midpath = args.temp_path

	# Tool path
	bedtools= args.bedtools_path

	# Forming file lists from given prefix:
	mutfiles = []
	strfiles = []
	germfiles = []
	samples = args.samples_names
	for idx, sample in enumerate(samples):
		mutfiles.append(args.mutect2_prefix+f'/{sample}{args.samples_suffix_tumor}_vs_{sample}{args.samples_suffix_normal}/{sample}{args.samples_suffix_tumor}_vs_{sample}{args.samples_suffix_normal}.mutect2.somatic_snvs_snpEff.ann.vcf')
		strfiles.append(args.mutect2_prefix+f'/{sample}{args.samples_suffix_tumor}_vs_{sample}{args.samples_suffix_normal}/{sample}{args.samples_suffix_tumor}_vs_{sample}{args.samples_suffix_normal}.strelka.somatic_snvs_snpEff.ann.vcf')
		germfiles.append(args.mutect2_prefix+f'/{sample}{args.samples_suffix_normal}/{sample}{args.samples_suffix_normal}.haplotypecaller.filtered_snpEff.ann.vcf.gz')

	# Main: Loop through samples, perform filtering and writing
	for idx,sample in enumerate(samples):
		print("\n-----\nPROCESSING: ",sample)

	# Get the names of the input files
		mut_filename = mutfiles[idx].split('input/')[1]
		mut_filename = mut_filename.split('.vcf')[0]
		str_filename = strfiles[idx].split('input/')[1]
		str_filename = str_filename.split('.vcf')[0]

	# Select autosomes for mutect and strelka vcfs
		c_mut_autosomes = "awk '$1 ~ /^#/ || $1 !~/[[:alpha:]]/' {} > {}{}_CHR.vcf".format(mutfiles[idx],midpath,mut_filename)
		c_str_autosomes = "awk '$1 ~ /^#/ || $1 !~/[[:alpha:]]/' {} > {}{}_CHR.vcf".format(strfiles[idx],midpath,str_filename)
		os.system(c_mut_autosomes)
		os.system(c_str_autosomes)

	# Remove germline variants
		c_mut_germline = "{} intersect -v -wa -header -a {}{}_CHR.vcf -b {} > {}{}_CHR_SOM.vcf".format(bedtools,midpath,mut_filename,germfiles[idx],midpath,mut_filename)
		c_str_germline = "{} intersect -v -wa -header -a {}{}_CHR.vcf -b {} > {}{}_CHR_SOM.vcf".format(bedtools,midpath,str_filename,germfiles[idx],midpath,str_filename)
		os.system(c_mut_germline)
		os.system(c_str_germline)

	# Remove indels (Mutect2 only)
		c_mut_indels_1 = "awk '$1 ~ /^#/' {}{}_CHR_SOM.vcf > {}{}_CHR_SOM_SNV.vcf".format(midpath,mut_filename,midpath,mut_filename)
		c_mut_indels_2 = "awk 'length($4) == 1 && length($5) == 1' {}{}_CHR_SOM.vcf >> {}{}_CHR_SOM_SNV.vcf".format(midpath,mut_filename,midpath,mut_filename)
		os.system(c_mut_indels_1)
		os.system(c_mut_indels_2)

	# FILTER is PASS (Strelka only)
		c_str_filter = "awk '$1 ~ /^#/ || $7 == \"PASS\"' {}{}_CHR_SOM.vcf > {}{}_CHR_SOM_PASS.vcf".format(midpath,str_filename,midpath,str_filename)
		os.system(c_str_filter)

	# FILTER is NOT 'germline' OR 'panel_of_normals' (Mutect2 only)
		c_mut_filter_1 = "awk '$1 ~ /^#/' {}{}_CHR_SOM_SNV.vcf > {}{}_CHR_SOM_SNV_FILTER.vcf".format(midpath,mut_filename,midpath,mut_filename)
		c_mut_filter_2 = "awk '$7 !~ /panel_of_normals/ && $7 !~ /germline/' {}{}_CHR_SOM_SNV.vcf >> {}{}_CHR_SOM_SNV_FILTER.vcf".format(midpath,mut_filename,midpath,mut_filename)
		os.system(c_mut_filter_1)
		os.system(c_mut_filter_2)

	# Intersect the two files and keep mutect calls
		c_intersect = "{} intersect -wa -header -a {}{}_CHR_SOM_SNV_FILTER.vcf -b {}{}_CHR_SOM_PASS.vcf > {}{}_filtered_somatic.vcf".format(bedtools,midpath,mut_filename,midpath,str_filename,outpath,sample)
		os.system(c_intersect)




	# Count in outputfile
	print("\n-----\nCOUNTING\n")
	out = open(outfile,'w+')

	# Print header
	header = "Sample\tCaller\tUnfiltered\tAutosomal\tNon-germline\tSNV\tFilter\tFinal"
	out.write(header)

	# For each sample
	for idx,sample in enumerate(samples):

	# Get the names of the input files
		mut_filename = mutfiles[idx].split('input/')[1]
		mut_filename = mut_filename.split('.vcf')[0]
		str_filename = strfiles[idx].split('input/')[1]
		str_filename = str_filename.split('.vcf')[0]

	# Count all Mutect variants
		mut1 = int(subprocess.check_output("awk '$1 !~ /^#/' {} | wc -l".format(mutfiles[idx]),shell= True))
		mut2 = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_CHR.vcf | wc -l".format(midpath,mut_filename), shell=True))
		mut3 = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_CHR_SOM.vcf | wc -l".format(midpath,mut_filename), shell=True))
		mut4 = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_CHR_SOM_SNV.vcf | wc -l".format(midpath,mut_filename), shell=True))
		mut5 = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_CHR_SOM_SNV_FILTER.vcf | wc -l".format(midpath,mut_filename), shell=True))

	# Count all Strelka variants
		str1 = int(subprocess.check_output("awk '$1 !~ /^#/' {} | wc -l".format(strfiles[idx]),shell= True))
		str2 = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_CHR.vcf | wc -l".format(midpath,str_filename), shell=True))
		str3 = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_CHR_SOM.vcf | wc -l".format(midpath,str_filename), shell=True))
		str4 = "NA"
		str5 = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_CHR_SOM_PASS.vcf | wc -l".format(midpath,str_filename), shell=True))

	# Count final variants
		final = int(subprocess.check_output("awk '$1 !~ /^#/' {}{}_filtered_somatic.vcf | wc -l".format(outpath,sample), shell = True))

	# Print
		result = "\n"
		result = "{}{}\tmutect2\t{}\t{}\t{}\t{}\t{}\t{}".format(result,sample,mut1,mut2,mut3,mut4,mut5,final)
		result = "{}\n{}\tstrelka\t{}\t{}\t{}\t{}\t{}\tNA".format(result,sample,str1,str2,str3,str4,str5)
		out.write(result)

	# Close outputfile
	out.close()




	# Remove intermediate files
	print("\n-----\nREMOVING INTERMEDIATE FILES = ON\n")

	remfiles = os.listdir('{}/'.format(midpath))

	for remfile in remfiles:
		os.remove("{}/{}".format(midpath,remfile))
