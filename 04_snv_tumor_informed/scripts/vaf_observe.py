import pandas as pd
import argparse
# python vaf_observe.py --query_map_qual {params.query_map_qual} --name {params.sample_repeated} \
#         --pickle {input.pickle} --tissue_tf {input.vcf_tf} \
#         --bed {output.bed} --vaf_observe {output.txt}


def write_observe(pickle, bedfile, vaf_observe, vaf_before_adjust, sample):
    a = pd.read_pickle(pickle)
    a.drop_duplicates(inplace = True)
    a = a[a['query_map_qual'] >= args.query_map_qual]
    a = a[a['CHROM'] != 'X']
    a['Sample'] = sample
    a= a.reset_index(drop = True)
    codes = {'REF':0, 'ALT':1, 'ERR':2}
    a['snv_type'] = a['snv_overlap_type'].map(codes)
    # Disable write header
    #with open(a_path, "w") as f:
    #         f.write(f"#num_measurements={df_combined[df_combined[sample]>=0].shape[0]}\n")
    df_vaf = a[['ALT_freq', 'snv_type']].sort_values('ALT_freq')
    # Write to VAF txt file without adjustment of VAF
    df_vaf.to_csv(args.vaf_before_adjust,
                    header=False,
                    index = False,
                    sep = "\t",
                    mode='w')

    # Adjust VAF based on tissue tumor purity. If above 1.0, adjust to 1.0
    df_vaf['ALT_freq'] = df_vaf['ALT_freq'].apply(lambda x: 1.0 if x / args.tissue_tf > 1.0 else x / args.tissue_tf  )
    df_vaf.to_csv(args.vaf_observe,
                    header=False,
                    index = False,
                    sep = "\t",
                    mode='w')

    # Write the variants ALT discovered in targeted sample as BED file. (or vcf file)
    ALT_write = a[a['snv_overlap_type'] == 'ALT'][['CHROM','start','end','REF','ALT','QUAL','FILTER','ALT_freq']]
    ALT_write.to_csv(bedfile, sep = '\t',index=False)





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-q",
        "--query_map_qual",
        type=int,
        help="Minimal mapping quality of the base of interest in cfDNA",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        help="Sample name",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--pickle",
        type=str,
        help="Path to pickle file",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--tissue_tf",
        type=float,
        default=1.0,
        help="tissue tumor purity for VAF adjustment",
    )
    parser.add_argument(
        "-a",
        "--bed",
        type=str,
        help="Output bed file of observed variants in cfDNA.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--vaf_observe",
        type=str,
        help="Output VAF per observed variants after adjustment of tissue purity.",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--vaf_before_adjust",
        type=str,
        help="Output VAF per observed variants before adjustment of tissue purity.",
        required=True,
    )


    args = parser.parse_args()
    # perform write
    write_observe(args.pickle, args.bed, args.vaf_observe, args.vaf_before_adjust, args.name)
