import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from math import floor
import numpy as np
from matplotlib.patches import Rectangle
import numpy as np
import argparse
import pandas as pd
from scipy.stats import binom
import itertools
import scipy.stats as stats
import matplotlib.pyplot as plt
import os
import cyvcf2
from scipy.signal import resample
import time
from tqdm import tqdm
#import logging
import sys
from numba import jit



@jit(nopython=True, parallel=False)
# Input x: vafs* tfs
# Output x: an array of 1,0 based on the vafs* tfs
def draw_allele(x):
    for i in range(x.shape[0]):
        x[i] = np.random.binomial(n = 1, p = x[i])
    return x


def get_tf(n_alt, df_all, df_mrd):
    print("Number of alt allele", n_alt)
    if n_alt <= 50:
        slice = df_mrd[n_alt]
        print('Mode: low ALT COV, MRD')
        mode = 'mrd'
    elif n_alt > 50:
        n_alt_bin = floor(n_alt / 10) * 10
        slice = df_all[n_alt_bin]
        print('Mode: high ALT COV')
        mode = 'high_tf'
    return slice, mode


def cumulative_sum_percentiles(input_array):
    cum_sum = np.cumsum(input_array)
    cum_sum_before_and_last = np.insert(cum_sum, 0, 0)
    if cum_sum.max() > 0:
        percentiles_array = cum_sum_before_and_last / cum_sum.max()
    else:
        percentiles_array =  np.zeros(len(input_array))
    return percentiles_array


def find_percentile_indices(input_slice, low_percentile, high_percentile):
    """
    To test if the selecting of index is correct:

    :param input_slice: An np.array([[tfs],[]])
    :param low_percentile:
    :param high_percentile:
    :return:
    """
    input_slice = input_slice.reset_index()
    input_slice.columns = ["TF","count"]
    input_array = np.array(input_slice['count'])
    if np.all(input_array == 0):
        percentile_5 = -1
        percentile_95 = -1
        estimated_tf = -1
        index_low = -1
        index_high = -1
        print("Warning: No TF bin contains possible value. Simulate larger range of TF...is a bug if you simulate from 0-1.")


    else:
        estimated_tf_id = input_slice["count"].idxmax()

        estimated_tf = input_slice['TF'][estimated_tf_id]
        ## Start with double tail statistics:
        # Percentile array always have 0 at the first position, because it represent the percentile "before" the first bin. It contains 1 more bin than input_array.
        percentiles_array = cumulative_sum_percentiles(input_array)
        # At which tumor fraction, the percentile is smaller than 2.5% on each tail
        # index of the percentile array
        index_low = np.where((percentiles_array ) < low_percentile)[0][-1] + 1 -1   # +1 for The bin right to the bin with less than tail. Then -1 for the actual bin in TFs. (percentile include left edge as an extra bin).
        index_high = np.where(percentiles_array > high_percentile )[0][0]  - 1 -1  # -1 Take bin left of the bin with tail percentile. Then -1 for the actual bin in TFs. (percentile include left edge as an extra bin).
        if index_low == 1:
            print("Warning: lowest TF estimate limit is at simulation limit. Simulate lower range. (unless lowest TF is at 0).")
        if index_high == len(input_slice):
            print("Warning: Highest TF estimate limit is at simulation limit. Simulate higher range.  (unless highest TF is already at 1.0).")
        if estimated_tf_id > index_low:
            percentile_5 = input_slice['TF'][index_low]
            percentile_95 = input_slice['TF'][index_high]
        elif estimated_tf_id <= index_low:
            print("Max bin is smaller than left tail limit. Use single-tailed statistics instead.")
            percentile_5 = input_slice['TF'][0]
            single_tail_high_percentile = high_percentile - low_percentile
            index_high_single_tail = np.where(percentiles_array > single_tail_high_percentile)[0][0] -1 -1
            if index_high_single_tail < 0:
                # An index to an array cannot be negative. In case the percentiles_array is only 2 values (input_slice only 1 value), then enter this exception.
                index_high_single_tail = 0
            percentile_95 = input_slice['TF'][index_high_single_tail]
            index_low = 0
            index_high = index_high_single_tail

        print("TF max probability:", estimated_tf)
        print("\nTF 5% Percentile:", percentile_5)
        print("TF 95% Percentile:", percentile_95)



    return percentile_5, percentile_95,  estimated_tf, index_low, index_high


def plot_monte_carlo_tf_estimate(df_all, df_mrd, n_alt, index_low, index_high, estimated_tf, mode, outpath_fig):
    # plot all
    if float(estimated_tf) < 0:
        plt.subplots(figsize=(15, 6))
        plt.savefig(outpath_fig)
    else:
        if mode == 'high_tf':
            try:
                right_most_column_with_value = int(df_all.columns[df_all.sum() > 0].max())
                df2 = df_all.loc[:, :right_most_column_with_value + 1]
            except Exception as e:
                print(e, "right_most_column_with_value")
                df2 = df_all.copy()
            plt.subplots(figsize=(15, 6))
            ## TODO: x = 0 is empty. Shouldn't be.
            g = sns.heatmap(df2, cmap='Greys')
            ## Get x axis for n_alt or n_alt bin
            number_of_observed_alt_bin = floor(n_alt / 10) * 10
            slice_high_tf = df_all[number_of_observed_alt_bin]
            ## Get x, y axis for Estimated TF on heatmap
            x_index_all = df_all.columns.get_loc(number_of_observed_alt_bin)
            y_index_all = df_all.index.get_loc(estimated_tf)
            ## Plot red box around the highest estimated TF
            g.add_patch(Rectangle((x_index_all, index_low), 1, index_high-index_low+1, fill=False, edgecolor='orange', lw=1))
            g.add_patch(Rectangle((0, y_index_all), x_index_all, 1, fill=False, edgecolor='blue', lw=1))
            g.add_patch(Rectangle((x_index_all, y_index_all), 1, 1, fill=False, edgecolor='red', lw=1))
            plt.savefig(outpath_fig)
            # plt.show()
        # plot mrd
        if mode == 'mrd':
            plt.subplots(figsize=(15, 6))
            df_mrd = df_mrd.loc[(df_mrd.sum(axis=1) > 0),]

            g2 = sns.heatmap(df_mrd, cmap='Greys')
            # Add box around

            ## Get x axis for n_alt, ## Get y axis for Estimated TF
            x_index = df_mrd.columns.get_loc(n_alt)
            y_index = df_mrd.index.get_loc(estimated_tf)
            ## Plot Red box around the highest probable TF.
            g2.add_patch(Rectangle((x_index, y_index), 1, 1, fill=False, edgecolor='red', lw=3))
            # Orange block at 5%-95% confidence
            g2.add_patch(Rectangle((x_index, index_low), 1, index_high-index_low+1, fill=False, edgecolor='orange', lw=1))
            # blue lines to align to Y axis
            g2.add_patch(Rectangle((0, y_index), x_index, 1, fill=False, edgecolor='blue', lw=1))
            plt.savefig(outpath_fig)


def check_path(output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)


def store_df(df_detection_per_id, output_path, df_types= ['pickle', 'csv',]):
    if 'pickle' in df_types:
        df_detection_per_id.to_pickle(f"{output_path}.pickle.gz")
    if 'csv' in df_types:
        df_detection_per_id.to_csv(f"{output_path}.csv", sep = '\t')


def simulate_VAF_dist(VAF_mean, VAF_sigma, trials):
    lower, upper = 0.0, 1.0
    mu, sigma = VAF_mean, VAF_sigma
    X = stats.truncnorm(
        (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    # N = stats.norm(loc=mu, scale=sigma)
    VAF = X.rvs(trials)
    # fig, ax = plt.subplots(2, sharex=True)
    # ax[0].hist(X.rvs(10000), )
    # ax[1].hist(N.rvs(10000), )
    return VAF


def initialize_df(rejecting_threshold_percentile=[0, 0.680, 0.900, 0.950, 0.997]):
    names = []
    for percentile in rejecting_threshold_percentile:
        names += [f"rejecting_count-t{percentile}",
                  f"TP_trial_ratio-t{percentile}",
                  f"FN_trial_ratio-t{percentile}",
                  f"FP_trial_ratio-t{percentile}",
                  f"TN_trial_ratio-t{percentile}",
                  f"detect_rate-t{percentile}",
                  f"false_discovery_rate-t{percentile}",
                  ]
    # Initiate data collection dataframe. Each row stores one input condition (combination of parameters).
    df_detection_per_id = pd.DataFrame(columns=["sample",
                                                "SNV_site",
                                                "SNV_site_covered",
                                                "sequencing_depth",
                                                "tumor_fraction",
                                                "sequencing_error_rate",
                                                "trials_count",
                                                "TP_snv_count_mean",
                                                "FN_snv_count_mean",
                                                "FP_snv_count_mean",
                                                "TN_snv_count_mean",
                                                ]+ names)

    return df_detection_per_id


def determine_confusion_matrix(snv_positive, snv_false_negative, snv_false_positive, rejecting_threshold=0):
    """
    snv_positive =       np.array([9,1,2,0,0,2])
    snv_false_negative = np.array([0,1,0,0,0,2])
    snv_false_positive = np.array([0,0,1,1,0,2])

    # tp: true observation
    # fn: flip back observation
    # fp: error observation, independent of tp observation
    rejecting_threshold = 0
    observe_tp = np.array([1,0,1,0,0,1])  # 3
    observe_fn = np.array([0,1,0,0,0,0])  # 1
    observe_fp = np.array([0,0,0,1,0,0])  # 1
    observe_tn = np.array([0,0,0,0,1,0])  # 1
    observe_p =  [1,0,1,1,0,1]
    observe_n =  [0,1,0,0,1,0]
    ## Test 2
    rejecting_threshold > 2
    snv_positive =       np.array([9,1,2,0,0,2])
    snv_false_negative = np.array([0,1,0,0,0,2])
    snv_false_positive = np.array([0,0,1,1,0,2])

    observe_tp = np.array([1,0,1,0,0,0])  # 2
    observe_fn = np.array([0,1,0,0,0,1])  # 2
    observe_fp = np.array([0,0,0,0,0,0])  # 0
    observe_tn = np.array([0,0,0,1,1,0])  # 1

    ## Definition of true positive is: truth is positive (at least 1 SNV), measured positive
    ## Definition of false negative is: truth is positive, measured negative.
    ## Definition of false positive is: truth is negative, measured positive
    ## Definition of true negative is: truth is negative, measured negative.

    # Measure positive or negative: (TP + FP - FN)  > rejecting_threshold
    # Truth is positive or negative; TP > 0

    """
    measurement = snv_positive + snv_false_positive - snv_false_negative
    truth = snv_positive
    tp_trial_count = (truth > 0 ) * (measurement > rejecting_threshold)
    fn_trial_count = (truth > 0 ) * (measurement <= rejecting_threshold)
    fp_trial_count = (truth == 0 ) * (measurement > rejecting_threshold)
    tn_trial_count = (truth == 0 ) * (measurement <= rejecting_threshold)

    return tp_trial_count, fn_trial_count, fp_trial_count, tn_trial_count


def simulate(params_combinations,
             SNV_site,
             df_detection_per_id,
             VAF,
             sample=None,
             df_variant_count=None,
             df_variant_count_mrd= None,
             tf_estimation=False,
             rejecting_threshold_percentile=[0.680, 0.900, 0.950, 0.997]
             ):
    pbar = tqdm(total=len(params_combinations),
                desc='collection simulation')
    # for each combination, perform simulation
    existing_rows = df_detection_per_id.shape[0]
    detection100_min_tf_per_GE_error_rate = {}
    for row_count, comb in enumerate(params_combinations):
        pbar.update(1)
        # unzip parameters from 'params_combinations'
        GE, tf, error_rate = comb
        # Without resampling (each SNV can only be sequenced once)
        sequencing_depth = GE

        # How many unique SNV are covered in each trial? (an approximate average, not a distribution)
        # Fixed number for each sequencing depth.
        all_unique_site_count_per_trial = int(SNV_site * sequencing_depth)
        # Array
        if sequencing_depth < 1:
            all_sites_counts = np.random.binomial(SNV_site, sequencing_depth, size=trials)
        # in case of tf estimate:
        elif sequencing_depth == 1:
            all_sites_counts = [SNV_site]* trials
        else:
            # Sequencing depth >1:
            # Coverage per site is determined as a normal distribution around the coverage depth
            cov_per_site = np.random.normal(loc=sequencing_depth, scale=sequencing_depth/2, size=(SNV_count, trials))
            cov = np.sum(cov_per_site, axis = 0)
            all_sites_counts = np.rint(cov).astype(np.int32)

        # Determine false positive and false negative error rate.
        fp_seq_rate = error_rate / 3
        fn_seq_rate = error_rate
        #### note: Each of the following parameters are all a array of 100_000 values. ####
        ## Calculate how many positive sites do I cover in a given parameter set
        snv_positive = np.zeros((trials,), dtype=int)
        # TODO: Try removing this for loop to make the program much quicker. Draw for 10000 trials instead of per trial

        for trial in range(trials):
            if len(VAF) > all_sites_counts[trial]:
                vaf_subset = np.random.choice(VAF, size = all_sites_counts[trial], replace=False)
                assert np.all(vaf_subset * tf <= 1), (vaf_subset*tf)[vaf_subset*tf > 0 ]
            # In case of tf_estimate
            elif len(VAF) == all_sites_counts[trial]:
                vaf_subset = VAF
            else:
                # Choose VAF from vaf list based on amount of sites.
                vaf_subset = np.random.choice(VAF, size = all_sites_counts[trial], replace=True)
            # This line below is the time taking step.
            assert len(vaf_subset) > 0 , "Check input VAF file if it is empty!"
            snv_positive[trial] = draw_allele(vaf_subset * tf).sum()
            # The rest of the measured SNV are reference allele (snv_negative)
        snv_negative = all_sites_counts - snv_positive

        # Based on the sequencing error rate, calculate out of all positive SNV, how many would be observed as "negative" (reference allele or any other alleles.)
        fn = np.random.binomial(snv_positive, fn_seq_rate, trials)
        # The remaining are real true positive, measured positive (ALT allele).
        # Based on the sequencing error rate, calculate out of all SNV positions where REF allele is measured, how many would be observed as "ALT" (ALT allele only)
        #assert np.random.binomial(snv_negative, fp_seq_rate, trials), f"GE:{GE}, SNV_site:{SNV_site}, tf:{tf}, error_rate:{error_rate}, sample:{sample}"
        try:
            fp_count_per_trial = np.random.binomial(snv_negative, fp_seq_rate, trials)
        except ValueError as e:
            print("fp got value error", e )
            print( f"snv_negative has negative values: {snv_negative[snv_negative < 0]} , GE:{GE}, SNV_site:{len(VAF)}, tf:{tf}, error_rate:{error_rate}, sample:{sample}")
        # Remaining are real true negative
        tp_trial_count, fn_trial_count, fp_trial_count, tn_trial_count = determine_confusion_matrix(snv_positive, fn, fp_count_per_trial)


        # In a trial, any SNV detected is defined as "Detection". Regardless if it is detection of true SNV (TP) or false SNV (FP).
        detection = (tp_trial_count + fp_count_per_trial) > 0
        detection = detection.astype(int)

        detect_at_least_1_fp = fp_count_per_trial > 0
        detect_at_least_1_fp = detect_at_least_1_fp.astype(int)
        # false detection is a binary true/false array.
        # tp>0 is also a binary true/false array.
        # false discovery is also a boolean. True and False.
        false_discovery = detect_at_least_1_fp - (snv_positive>0) > 0
        # fn:  This is NOT correct, did not take into account of False positive, these would then be false positive trials, instead of false negative.
        # of all detection instance, how many of these detection is actually false detection? (no tp, only fp).
        false_discovery_rate = np.sum(false_discovery) / np.sum(detection)
        # TODO: locations within it where the condition is False will remain uninitialized. # check where tp+ fp = 0 if false_discovery = Nan

        # False negative
        confusion_matrix = []
        percent_reject_null_count= []
        for percentile in rejecting_threshold_percentile:
            rejecting_threashold_X_percentile = binom.ppf(percentile, all_unique_site_count_per_trial, error_rate)

            tp_r, fn_r, fp_r, tn_r= determine_confusion_matrix(snv_positive,
                                       fn,
                                       fp_count_per_trial,
                                       rejecting_threshold=rejecting_threashold_X_percentile)
            percentage_reject_null_X = ( np.sum(tp_r) + np.sum(fp_r) ) / trials
            percent_reject_null_count.append(percentage_reject_null_X)
            # nan often
            if (np.sum(tp_r) + np.sum(fp_r)) == 0:
                # If there is no discovery/detection, there is no false discovery rate. Set it to max (=1)
                fdr_threshold_X = 1
            else:
                fdr_threshold_X = np.sum(fp_r) / (np.sum(tp_r) + np.sum(fp_r))

            confusion_matrix += [rejecting_threashold_X_percentile,
                                     np.sum(tp_r) / trials,
                                     np.sum(fn_r) / trials,
                                     np.sum(fp_r) / trials,
                                     np.sum(tn_r) / trials,
                                     percentage_reject_null_X,
                                     fdr_threshold_X,
                                     ]

        df_detection_per_id.loc[row_count+existing_rows] = [sample,
                                                            SNV_site,
                                                            all_unique_site_count_per_trial,
                                                            sequencing_depth,
                                                            tf,
                                                            error_rate,
                                                            trials,
                                                            np.mean(snv_positive),
                                                            np.mean(fn),
                                                            np.mean(fp_count_per_trial),
                                                            0,
                                                            0,
                                                            np.sum(tp_trial_count) / trials,
                                                            np.sum(fn_trial_count) / trials,
                                                            np.sum(fp_trial_count) / trials,
                                                            np.sum(tn_trial_count) / trials,
                                                            (np.sum(tp_trial_count) + np.sum(fp_trial_count)) / trials,
                                                            false_discovery_rate,
                                                            ] + confusion_matrix
        epsilon = 0.001  # A small value
        if ((np.sum(tp_trial_count) + np.sum(fp_trial_count)) / trials == 1) & (false_discovery_rate == 0):
            if (error_rate, GE) not in detection100_min_tf_per_GE_error_rate.keys():
                detection100_min_tf_per_GE_error_rate[(error_rate, GE)] = tf
                # print(detection100_min_tf_per_GE_error_rate)
        # 0 will become negative, therefore fall into the 1st bin (and only 0) is in this bin.
        # 10000 will become 9999.999 which fall into the last bin (from 9990-10000)
        alt_detection_per_trial_size_trial = (snv_positive + fp_count_per_trial - fn) - epsilon
        if tf_estimation == True:
            ## TO store monte carlo simulation per parameter combination per sample (disable for large sample size).
            #Exclude the last bin at n = 50 because it contains everything above. We add 1 bin (n=51) to the end, and remove it when
            # Storing into the dataframe.
            alt_count_minimal, header_mrd = np.histogram(
                alt_detection_per_trial_size_trial,
                bins=52,
                range=(-1, 51)
            )

            if df_variant_count_mrd is None:
                # Exclude the last bin (n = 51)
                df_variant_count_mrd = pd.DataFrame(columns=[int(x) for x in header_mrd[1:-1]])
                df_variant_count_mrd.loc[str(round(tf, 6))] = list(alt_count_minimal[:-1])
            else:
                df_variant_count_mrd.loc[str(round(tf, 6))] = list(alt_count_minimal[:-1])
            alt_count, header_all = \
                np.histogram(
                alt_detection_per_trial_size_trial,
                bins=1001,
                range=(-10, 10000)
            )
            if df_variant_count is None:
                df_variant_count = pd.DataFrame(columns=[int(x) for x in header_all[1:]])
                df_variant_count.loc[str(round(tf, 6))] = list(alt_count)
            else:
                df_variant_count.loc[str(round(tf, 6))] = list(alt_count)
        else:
            df_variant_count = None
            df_variant_count_mrd = None
    # print('Error rate & GEs and the min TF with DR=1, FDR=0')
    # print(detection100_min_tf_per_GE_error_rate)
    return df_detection_per_id, df_variant_count, df_variant_count_mrd


def check_vcf_format(input_vcf):
    if input_vcf.endswith("vcf"):
        print("VCF file is not indexed")
        print(
            f"bgzip -c {input_vcf} > {input_vcf}.gz",
        )
        print(f"tabix -p vcf {input_vcf}.gz")
        exit(1)


def get_vaf_from_vcf(somatic_vcf, output_file=None):
    somatic_vcf = cyvcf2.VCF(somatic_vcf)
    vaf = []
    vaf_list = []
    # var_type = collections.Counter()
    # var_subtype = collections.Counter()
    # aaf_collect = collections.Counter()
    for ssnv_raw in somatic_vcf:
        # Check indel length
        if ssnv_raw.is_snp:
            vaf.append(str(ssnv_raw.gt_alt_freqs[0]))
            vaf_list.append(ssnv_raw.gt_alt_freqs[0])
            # var_type[ssnv_raw.var_type] += 1
            # var_subtype[ssnv_raw.var_subtype] += 1
            # aaf_collect[ssnv_raw.aaf] += 1
            ## If there are >2 ALT alleles which contains one of them none snp, will not be included in this list (would be unknown)
    if output_file is not None:
        with open(output_file, "w") as f:
            f.write("\n".join(vaf))
    # vaf is a list of vaf for a certain vcf file.
    return vaf_list


def get_vaf_from_vaf_tables(vaf_table_path):
    directory = vaf_table_path
    vafs = {}
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if filename.endswith("txt") and filename.startswith("PCAWG"):
            vaf = pd.read_csv(f, sep="\t")
            name = filename.split("_")[1].strip(".txt")
            print(name)
            vafs[name] = vaf['VAF']
    return vafs


def get_vaf_from_single_vaf_table(vaf_table_txt_path):
    vaf = pd.read_csv(vaf_table_txt_path, sep="\t")
    return vaf

class Files:
    def __init__(self, params):
        self.params = params
        self.suffix = None
        self.id = None
        self.filename = None
        self.output_folder = None
        self.cancer_type = self.params['cancer_type']
        self.VAF_type = self.params['VAF_type']
        self.get_cancer_type()
        self.filename_handling()

    def get_cancer_type(self):
        if self.cancer_type == 'None' and self.VAF_type == 'combined':
            self.cancer_type = self.params['vaf_table_path'].split('/')[-1].split('.')[0]
        elif self.cancer_type == 'None' and self.VAF_type == 'VCF':
            self.cancer_type = self.params['somatic_vcf_path'].split('/')[-1].split('.')[0]
        elif self.cancer_type == 'None' and self.VAF_type == 'simulated':
            self.cancer_type = 'simulated'
        else:
            self.cancer_type = self.params['cancer_type']

    def filename_handling(self):
        self.id = self.params['id'].replace(" ", "_")
        self.suffix = self.params['VAF_type']
        self.output_folder = self.params['output_folder']
        self.filename = f"{self.cancer_type}_{self.id}_{self.suffix}"
        check_path(self.output_folder)


def determine_VAF_tf_estimate(tumor_fractions, error_rate, variant_detection_table_path, name):
    variant_detection_table = variant_detection_table_path
    # print(f"VAF_type mode: estiamte_tf. This mode estimate the tumor fraction of input sequencing run with a variant detection table.")
    VAF = pd.read_csv(variant_detection_table, sep = '\t', names=['vaf', 'obs'], comment='#')
    VAF = np.array(VAF.vaf)
    SNV_count = len(VAF)
    VAF_list = [(VAF, SNV_count, name)]
    list_of_params = [GEs, tumor_fractions, [error_rate]]
    params_combinations = list(itertools.product(*list_of_params))
    # Cap VAF between 0 and 1.
    VAF_list_updated = []
    for (VAF, SNV_sites, sample) in VAF_list:
        VAF[VAF < 0] = 0.0000000001
        VAF[VAF > 1] = 1.0
        VAF_list_updated.append((VAF, SNV_sites, sample))

    return params_combinations, VAF_list_updated


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e",
        "--error_rate",
        type=float,
        help="A sequencing error rate",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--type",
        type=str,
        help="define range of TF to monte carlo simulate and step type",
        default="linear",
    )

    parser.add_argument(
        "-i",
        "--variant_detection_table",
        type=str,
        help="Path to variant detected in sample with their VAF",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--outpath_figure",
        type=str,
        help="Path to write figure.",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--output_summary",
        type=str,
        help="Path to write output summary TF estimate",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_folder",
        type=str,
        help="Path to write output summary TF estimate",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        help="sample name",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--tumor_fraction_max",
        type=float,
        help="max tumor fraction to simulate",
        default=1.0,
    )
    args = parser.parse_args()
    cancer_type = args.name
    id = args.name
    output_folder = args.output_folder
    filename = args.name
    error_rate = args.error_rate

    # Check if vcf format is compatible. must be gzip and tabix.
    # Amount of trials:
    trials = 100_00
    genome_size = 2_806_377_226
    GEs = [1]
    tumor_fraction_min= 0.0
    tumor_fraction_max= args.tumor_fraction_max
    tumor_fraction_step= 100
    step_type = args.type

    if step_type == 'log':
        if tumor_fraction_min == 0:
            tumor_fraction_min = 1 / 10**tumor_fraction_step # /20
        tumor_fractions = np.logspace(np.log10(tumor_fraction_min), np.log10(tumor_fraction_max), num=tumor_fraction_step)
        tumor_fractions = np.array([0] + list(tumor_fractions) )
    elif step_type == 'linear':
        # np.linspace --> include start and end bin.
        tumor_fractions = np.linspace(tumor_fraction_min, tumor_fraction_max , num=tumor_fraction_step+1)
        tumor_fractions[tumor_fractions > 1] = 1.0
    else:
        print("For params step_type, please select from 'log' or 'range'. ")
        exit(1)

    # print(f"Tumor fractions selected: {tumor_fractions}")
    # fig, ax = plt.subplots()
    # ax.set_yticks(np.arange(0, 0.55, 0.05))
    # ax.yaxis.grid(color='gray', linestyle='dashed')
    # ax.stem(tumor_fractions, bottom=0.01)
    # ax.set_ylabel('tumor fraction')
    # ax.set_xlabel('simulation')
    # plt.savefig(f"{args.outpath_figure}-tf.pdf",dpi=250)
    # fig, ax = plt.subplots()
    # ax.yaxis.grid(color='gray', linestyle='dashed')
    # ax.stem(tumor_fractions[tumor_fractions<=0.01], bottom=0.001)
    # ax.set_ylabel('tumor fraction')
    # ax.set_xlabel('simulation')
    # plt.savefig(f"{args.outpath_figure}-tf-less-than-1percent.pdf",dpi=250)

    # initiate dataframe, collect data. For combined, collect data multiple times.
    df_detection_per_id = initialize_df(rejecting_threshold_percentile=[0, 0.680, 0.900, 0.950, 0.997])
    tf_estimation = True
    # Initialize df_variant_count
    df_variant_count_mrd = None
    df_variant_count = None
    # Handling different input VAF type.
    params_combinations, VAF_list = determine_VAF_tf_estimate(tumor_fractions, error_rate, args.variant_detection_table, args.name)
    # print(f"{len(VAF_list)} set of VAFs (samples) &")
    # print(f"{len(params_combinations)} combination parameters to simulate per VAF (sample).")
    for i, (VAF, SNV_count, sample) in enumerate(VAF_list):
        timea = time.time()
        # Simulation started
        df_detection_per_id, df_variant_count, df_variant_count_mrd = simulate(params_combinations,
                                                                               SNV_count,
                                                                               df_detection_per_id,
                                                                               VAF,
                                                                               sample,
                                                                               df_variant_count=df_variant_count,
                                                                               df_variant_count_mrd=df_variant_count_mrd,
                                                                               tf_estimation=tf_estimation,
                                                                               rejecting_threshold_percentile=[0.680, 0.900, 0.950, 0.997])
        print(f"\nsample:{sample}, SNV_count:{SNV_count}, {i+1}/{len(VAF_list)}. Took {time.time() - timea} seconds.")
    store_df(df_variant_count, f"{output_folder}/{filename}_{str(tumor_fraction_max)}_monte_carlo", df_types=['pickle', 'csv'])
    store_df(df_variant_count_mrd,  f"{output_folder}/{filename}_{str(tumor_fraction_max)}_monte_carlo_mrd", df_types=['pickle', 'csv'])
    # Store summary:
    store_df(df_detection_per_id, f"{output_folder}/{filename}_{str(tumor_fraction_max)}_all_simulations", df_types=['pickle', 'csv'])
    # print(f"Output stored at:{output_folder}/{filename}.pickle.gz, {output_folder}/{filename}.csv")

    df_variant_count.columns = df_variant_count.columns.astype(int)
    df_variant_count_mrd.columns = df_variant_count_mrd.columns.astype(int)

    x = pd.read_csv(args.variant_detection_table, sep = '\t', names=['vaf', 'obs'], comment='#')
    ## Obs = 2: error, = 1: ALT, = 0: REF
    n_err = x[x['obs'] == 2].shape[0]
    n_alt = x[x['obs'] == 1].shape[0]
    n_ref = x[x['obs'] == 0].shape[0]
    n_all = x.shape[0]
    # Set error allele to ref (not detecting ALT).
    x[x['obs'] == 2] = 0
    assert all(element in [0, 1] for element in x.obs)
    slice, mode = get_tf(n_alt, df_variant_count, df_variant_count_mrd)
    percentile_5, percentile_95, estimated_tf, index_low, index_high = find_percentile_indices(slice, 0.025, 0.975)
    plot_monte_carlo_tf_estimate(df_variant_count, df_variant_count_mrd, n_alt, index_low, index_high, estimated_tf, mode, args.outpath_figure)

    # Store tf_estimate:
    print("\nSummary: \nSample:", cancer_type, "\nMeasured variants:", n_all, "\nALT measurement:", n_alt, "\nTF:", estimated_tf)
    out_df = pd.DataFrame(columns=['sample',
                                   'all_count',
                                   'ALT_count',
                                   'REF_count',
                                   'ERR_count_regarded_as_REF',
                                   'SNV_derived_TF',
                                   'SNV_derived_TF_range_p0025',
                                   'SNV_derived_TF_range_p0975',
                                   'simulate_tf_min',
                                   'simulate_tf_max',
                                   'simulate_tf_steps']
                          )
    out_df.loc[1] = [cancer_type, n_all, n_alt, n_ref, n_err, estimated_tf, percentile_5, percentile_95,tumor_fraction_min, tumor_fraction_max , tumor_fraction_step]
    out_df.to_csv(args.output_summary, sep = "\t")
