import pandas as pd
import os
import argparse

def get_vaf_from_single_vaf_table(vaf_table_txt_path):
    vaf = pd.read_csv(vaf_table_txt_path, sep="\t")
    return vaf

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--vaf-path",
        type=str,
        help="Path to yaml parameters",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-path",
        type=str,
        help="Path to yaml parameters",
        required=True,
    )

    args = parser.parse_args()
    for a_file in os.listdir(args.vaf_path):
        print(a_file)
        if a_file.endswith('.csv'):
            cancer_type = a_file.split('.')[0]
            # print(cancer_type)
            vaf_table = get_vaf_from_single_vaf_table(os.path.join(args.vaf_path,a_file))
            os.makedirs(os.path.join(args.output_path, cancer_type))
            for sample in vaf_table["SAMPLE"].unique():
                #print(sample)
                vaf_list = list(vaf_table[vaf_table['SAMPLE'] == sample]['VAF'])
                vaf_subset = vaf_table[vaf_table['SAMPLE'] == sample]
                #pd.to_pickle(vaf_subset, os.path.join(args.output_path, cancer_type ,f'{cancer_type}_{sample}.pickle.gz') )
                with open(os.path.join(args.output_path, f'{cancer_type}_variant_count.txt'), 'a+') as f:
                    f.write(f'{cancer_type}\t{sample}\t{vaf_subset.shape[0]}\n')
    print(f'all files written to {args.output_path}')
