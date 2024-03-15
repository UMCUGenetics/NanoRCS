# %%
import pysam
import pandas as pd
import os
import time

# %%

header = '\t'.join(["SAMPLE","VAF","COV","GENE","VARIANT"])+'\n'

# List of variant types
type_coding = {'coding_sequence_variant', 
               'conservating_inframe_insertion', 'conservative_inframe_deletion', 
               'disruptive_inframe_deletion', 'disruptive_inframe_insertion', 
               'exon_loss', 'exon_loss_variant', 'exon_lost_variant', 'exon_variant', 
               'frameshift_truncation', 'frameshift_variant', 
               'incomplete_terminal_codon_variant', 'initiator_codon_variant',
               'inframe_deletion', 'inframe_insertion',
               'missense_variant', 'rare_amino_acid_variant', 
               'start_lost', 'start_retained', 'start_retained_variant', 
               'stop_gained', 'stop_lost', 'stop_retained_variant', 
               'synonymous_variant'}
type_splice = {'splice_acceptor_variant', 'splice_donor_5th_base_variant',
               'splice_donor_region_variant', 'splice_donor_variant',
               'splice_polypyrimidine_tract_variant', 'splice_region_variant'}
type_intron = {'intron_variant','conserved_intron_variant'}
type_stream = {'3_prime_UTR_truncation', '3_prime_UTR_variant', 
               '5_prime_UTR_premature_start_codon_gain_variant', 
               '5_prime_UTR_truncation', '5_prime_UTR_variant', 
               'TFBS_ablation', 'TFBS_amplification', 
               'TF_binding_site', 'TF_binding_site_variant', 
               'conserved_intergenic_variant', 'downstream_gene_variant', 
               'intergenic_region', 'intergenic_variant', 
               'internal_feature_elongation', 
               'mature_miRNA_variant', 'miRNA', 
               'non_coding_transcript_variant', 
               'regulatory_region', 'regulatory_region_ablation', 
               'regulatory_region_amplification', 'regulatory_region_variant', 
               'UTR_variant', 'upstream_gene_variant'}


all_types_real = set()
all_types_miss = set()
def process_vcf(vcf_in, sample_id, file_out):
    global all_types_real
    global all_types_miss
    start = time.time()

    # Find the tumor sample in the VCF, probably sample 1
    sample_t = -1
    for i,sample in enumerate(vcf_in.header.samples):
        if str(sample).endswith('T'):
            sample_t = i
            break
    assert(i>=0), "No sample ending with T found (no tumor samples in VCF?)"

    # vcf_in = pysam.VariantFile(path_vcf)
    all_types = set()

    # Only process autosomal chromosomes
    for chrom in range(1,23):
        # Progress report
        print(chrom)

        for rec in vcf_in.fetch(str(chrom)):
            # Only keep PASS and SNV (lengths of 1)
            if ('PASS' in rec.filter) and len(rec.ref)==1 and len('-'.join(rec.alts))==1:
                # Placeholder variables
                gene = '.'
                type = '.'
                types = set()

                # GENE
                # not quite clear what we're supposed to grab here, if ANN has info then look for SEW and then what?
                if 'ANN' in rec.info:
                    ann = rec.info['ANN']
                    if 'SEW' in rec.info:
                        sew = rec.info['SEW']
                        gene, types = sew[0], set(sew[2].split('&'))
                    elif 'intragenic_variant' not in '|'.join(ann).split('|'):
                        print("WARNING: NO GENE NAME FOUND IN SEW, GENE IS NA",ann)
                    # else:                    
                        # Join and split with '|' to combine list of strings and split strings in elements again
                        # assert('intragenic_variant' in '|'.join(ann).split('|')), "ERROR: NO GENE NAME FOUND IN SEW, GENE IS NA"
                        

                # TYPE
                # CODING VARIANTS
                if types & type_coding:
                    type = 'exon'
                # SPLICE VARIANTS
                elif types & type_splice:
                    type = 'splicesite'
                # INTRON VARIANTS
                elif types & type_intron:
                    type = 'intron'
                # Upstream and downstream variants
                elif types & type_stream:
                    type = '.'
                    gene = '.'
                # Anything we didn't catch so far
                else:
                    # Clearly we missed something here
                    all_types_miss |= types
                    # continue
                # Track all possible types we've seen
                all_types_real |= types

                # DRIVERS
                if rec.info['REPORTED']:
                    type += '.DRIVER'

                # Handle missing PURPLE_AF cases
                vaf = 'NA'
                if rec.info['PURPLE_AF']:
                    vaf = min(rec.info['PURPLE_AF'], 1) 
                cov = rec.samples[sample_t]['DP'] 

                gene_type = 'NA'
                if gene != '.':
                    # print(gene,type)
                    gene_type = gene+'_'+type

                # print('\t'.join([str(x) for x in [sample_id,vaf,cov,gene_type,'_'.join([str(chrom),str(rec.pos),rec.ref,'.'.join(rec.alts)])]]))
                file_out.write('\t'.join([str(x) for x in [sample_id,vaf,cov,gene_type,'_'.join([str(chrom),str(rec.pos),rec.ref,'.'.join(rec.alts)])]])+'\n')

            # break
    # Report types we didn't catch before
    print("Time to process:", time.time() - start)

# %%
# Test one file processing
# with open("./data/debug_sample.txt", "w") as file_out:
#     file_out.write(header)
#     process_vcf(pysam.VariantFile("./data/DO218428T.purple.somatic.postprocessed.vcf.gz"), 'DO218428T', file_out)

# %%
df_selected = pd.read_csv('./data/selected_samples_table.txt',sep='\t')
uniq_types = pd.unique(df_selected['cancer_type'])
path_to_samples = './data/'
path_to_output = './output/'
sample_suffix = 'T.purple.somatic.postprocessed.vcf.gz'

for cancertype in uniq_types:
    with open(path_to_output+cancertype.replace(' ','-')+".csv", "w") as file_out:
        file_out.write(header)
        print(cancertype)
        samples = df_selected[df_selected['cancer_type'] == cancertype]['sample_id']
        for sample in samples:
            print('Process:', path_to_samples+sample+sample_suffix)
            if os.path.exists(path_to_samples+sample+sample_suffix):
                process_vcf(pysam.VariantFile(path_to_samples+sample+sample_suffix),sample,file_out)
            else:
                print("Unable to find:", path_to_samples+sample+sample_suffix)
                # pass

# Check if we caught all type annotations properly:
print('All seen variant types:', all_types_real)
print('Possibly skipped types:', all_types_real - type_coding - type_splice - type_intron - type_stream)
print('Uncaught variant types:', all_types_miss)
# %%
