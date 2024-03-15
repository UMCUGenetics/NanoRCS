import numpy as np
import argparse
import pandas as pd
from scipy.stats import binom
import itertools
import scipy.stats as stats
import os
import time
import sys
from numba import jit

@jit(nopython=True, parallel=False)
# Input x: vafs* tfs
# Output x: an array of 1,0 based on the vafs* tfs
def draw_allele(x):
    for i in range(x.shape[0]):
        x[i] = np.random.binomial(n = 1, p = x[i])
    return x

# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# from cfdetect.draw_allele import draw_allele


# Fixed parameters
GENOME_SIZE = 3_117_275_501
GENOME_SIZE_xMB = GENOME_SIZE / 1000_000

def check_path(output_path):
    if not os.path.exists(output_path):
        print("python line 30: path did not exist")
        os.makedirs(output_path)
        print("path created.")
    else:
        print("path already existed, skip creating.")


def store_df_no_filename(df_detection_per_id, output_path, df_types= ['pickle', 'csv',]):
    if 'pickle' in df_types:
        df_detection_per_id.to_pickle(output_path)
    if 'csv' in df_types:
        df_detection_per_id.to_csv(output_path, index=False, header=False, sep = '\t')


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
    existing_rows = df_detection_per_id.shape[0]
    detection100_min_tf_per_GE_error_rate = {}
    for row_count, comb in enumerate(params_combinations):
        # unzip parameters from 'params_combinations'
        GE, tf, error_rate = comb
        # Without resampling (each SNV can only be sequenced once)
        sequencing_depth = GE

        # How many unique SNV are covered in each trial? (an approximate average, not a distribution)
        # Fixed number for each sequencing depth. This is used to simulate the background FP based on error rate and use this to determine detection.
        all_unique_site_count_per_trial = int(SNV_site * sequencing_depth)


        # How many sites are covered in a run. A site can be covered twice (2 observation) in case of COV>1x, they are
        # regarded as independent observation because they could be either derived from normal cells or cancer cells.
        if sequencing_depth < 1:
            all_sites_counts = np.random.binomial(SNV_site, sequencing_depth, size=trials)
        else:
            # Sequencing depth >1:
            # Coverage per site is determined as a normal distribution around the coverage depth,
            # TODO: replace scale with STD of runs
            cov_per_site = np.random.normal(loc=sequencing_depth, scale=sequencing_depth/2, size=(SNV_count, trials))
            cov = np.sum(cov_per_site, axis = 0)
            all_sites_counts = np.rint(cov).astype(np.int32)

        # Determine false positive and false negative error rate.
        # If truth is A>T. All combination of changes are {A,C,G,T} > {A,C,G,T}
        # FP rate is A,C,G > T
        # FN rate is T > A,C,G
        fp_seq_rate = error_rate / 3
        fn_seq_rate = error_rate
        #### note: Each of the following parameters are all a array of 100_000 values. ####
        ## Calculate how many positive sites do I cover in a given parameter set
        snv_positive = np.zeros((trials,), dtype=int)
        # TODO: Try removing this for loop to make the program much quicker. Draw for 10000 trials instead of per trial

        # In all trials, I decide a subset of variants (with certain VAF)
        for trial in range(trials):
            if len(VAF) > all_sites_counts[trial]:
                # To speed up: Try randomize the VAF, and take the first X.
                # np.randint very fast (indexes) and then take from the VAF. (this is with replacement)
                vaf_subset = np.random.choice(VAF, size = all_sites_counts[trial], replace=False)
                assert np.all(vaf_subset * tf <= 1), (vaf_subset*tf)[vaf_subset*tf > 0 ]
            else:
                # Choose VAF from vaf list based on amount of sites.
                vaf_subset = np.random.choice(VAF, size = all_sites_counts[trial], replace=True)
            # This line below is the time taking step.
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
        # Truth= 5 TP, 495 TN.  What I get is 1 FN, 1 FP --> I get 5TP, 495TN.
        # Per trial: We detect cancer or not?
        # Make a confusion matrix per patient simulation ( out of 10_000 trials)
        tp_trial_count, fn_trial_count, fp_trial_count, tn_trial_count = determine_confusion_matrix(snv_positive, fn, fp_count_per_trial)


        # In a trial, any SNV detected is defined as "Detection". Regardless if it is detection of true SNV (TP) or false SNV (FP).
        detection = (tp_trial_count + fp_count_per_trial) > 0
        detection = detection.astype(int)

        detect_at_least_1_fp = fp_count_per_trial > 0
        detect_at_least_1_fp = detect_at_least_1_fp.astype(int)
        # false detection is a binary true/false array.
        # tp>0 is also a binary true/false array.
        # false discovery is also a boolean. True and False. FDR = (FP / (TP+FP))
        false_discovery = detect_at_least_1_fp - (snv_positive>0) > 0
        # fn:  This is NOT correct, did not take into account of False positive, these would then be false positive trials, instead of false negative.
        # of all detection instance, how many of these detection is actually false detection? (no tp, only fp).
        false_discovery_rate = np.sum(false_discovery) / np.sum(detection)
        # TODO: locations within it where the condition is False will remain uninitialized. # check where tp+ fp = 0 if false_discovery = Nan

        # False negative
        confusion_matrix = []
        percent_reject_null_count= []
        # This should be outside of the loop, each technique per threshold percentile  I get 1 cutoff
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
            ## TODO: report false negative rate per patient ( (1-tp_r) / (tp_r + fn_r)) --> How many patients we missed with this threshold.


            confusion_matrix += [rejecting_threashold_X_percentile,
                                     np.sum(tp_r) / trials,
                                     np.sum(fn_r) / trials,
                                     np.sum(fp_r) / trials,
                                     np.sum(tn_r) / trials,
                                     percentage_reject_null_X,
                                     fdr_threshold_X,
                                     ]
        ## TODO: speed up
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
    print('Error rate & GEs and the min TF with DR=1, FDR=0')
    print(detection100_min_tf_per_GE_error_rate)
    return df_detection_per_id, df_variant_count, df_variant_count_mrd



def check_vcf_format(input_vcf):
    if input_vcf.endswith("vcf"):
        print("VCF file is not indexed")
        print(
            f"bgzip -c {input_vcf} > {input_vcf}.gz",
        )
        print(f"tabix -p vcf {input_vcf}.gz")
        exit(1)


def get_vaf_from_single_vaf_table_pickle(vaf_table_txt_path):
    vaf = pd.read_pickle(vaf_table_txt_path)
    return vaf



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


def determine_VAF_smk(VAF_type, args, GEs, tumor_fractions, error_rates):
    if VAF_type == 'combined':
        print(f"You have supplied VAF per sample.")
        vaf_table = get_vaf_from_single_vaf_table_pickle(args.variant_detection_table)
        VAF_list = []
        for sample in vaf_table["SAMPLE"].unique():
            vaf_list = list(vaf_table[vaf_table['SAMPLE'] == sample]['VAF'])
            # resample VAF to count of trials
            SNV_count= len(vaf_list)
            VAF_list.append((np.array(vaf_list), SNV_count, sample))
    else:
        print("VAF_type should be one of these: simulated, combined, VCF, estimate_tf.")
        exit(1)

    list_of_params = [GEs, tumor_fractions, error_rates]
    params_combinations = list(itertools.product(*list_of_params))

    # Cap VAF between 0 and 1.
    VAF_list_updated = []
    for (VAF, SNV_sites, sample) in VAF_list:
        VAF[VAF < 0] = 0.0000000001
        VAF[VAF > 1] = 1.0
        VAF_list_updated.append((VAF, SNV_sites, sample))

    return params_combinations, VAF_list_updated

def get_tfs(tumor_fraction_min, tumor_fraction_max, tumor_fraction_step, step_type):
    if step_type == 'log':
        # TODO: The range of TF estimate is not correct for type log.
        if tumor_fraction_min == 0:
            tumor_fraction_min = 1 / 10**tumor_fraction_step
        tumor_fractions = np.logspace(np.log10(tumor_fraction_min), np.log10(tumor_fraction_max), num=tumor_fraction_step)
        tumor_fractions = np.array([0] + list(tumor_fractions) )
    elif step_type == 'linear':
        # np.linspace --> include start and end bin.
        tumor_fractions = np.linspace(tumor_fraction_min, tumor_fraction_max , num=tumor_fraction_step+1)
        tumor_fractions[tumor_fractions > 1] = 1.0
    else:
        print("For params step_type, please select from 'log' or 'range'. ")
        exit(1)
    return tumor_fractions

if __name__ == "__main__":
    TRIALS = 100_00
    genome_size = 3_117_275_501

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--variant_detection_table",
        type=str,
        help="folder containing vaf table per sample",
        required=True,
    )
    parser.add_argument(
        "-op",
        "--output_pickle_path",
        type=str,
        help="folder to write per sample output",
        required=True,
    )
    parser.add_argument(
        "-oc",
        "--output_csv_path",
        type=str,
        help="folder to write per sample output",
        required=False,
    )

    parser.add_argument(
        "-g",
        "--ge",
        type=float,
        help="genome equivalent",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--error-rate",
        type=float,
        help="error rate",
        required=True,
    )
    parser.add_argument(
        "-tfmin",
        "--tumor_fraction_min",
        type=float,
        help="tumor fraction min",
        required=True,
    )
    parser.add_argument(
        "-tfmax",
        "--tumor_fraction_max",
        type=float,
        help="tumor fraction max",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--tumor_fraction_step",
        type=int,
        help="tumor fraction steps",
        required=True,
    )
    parser.add_argument(
        "-tft",
        "--step-type",
        type=str,
        help="tumor fraction type",
        default='log',
        required=True,
    )
    args = parser.parse_args()
    prefix = args.variant_detection_table.split('/')[-1].split('.')[0]
    filename = f'{prefix}_{args.ge}_{args.error_rate}'
    print(filename)
    # YAML
    trials = TRIALS
    # Simulation parameters: sequencing related
    # Parallel by SMK. Each python runtime: only 1 error rate and 1 ge
    error_rates = [args.error_rate]
    GEs = [args.ge]

    # For 1 million mapped bam reads, We get 0.055x coverage (5.5% genome covered) on chr1-22.
    ## Stats is calculated from GCT01,GCT02.
    # GEs resembles the coverage. 0.05 per million reads.
    tumor_fractions = get_tfs(args.tumor_fraction_min, args.tumor_fraction_max, args.tumor_fraction_step, args.step_type)

    # initiate dataframe, collect data. For combined, collect data multiple times.
    df_detection_per_id = initialize_df(rejecting_threshold_percentile=[0, 0.680, 0.900, 0.950, 0.997])
    tf_estimation = False
    # Initialize df_variant_count
    df_variant_count_mrd = None
    df_variant_count = None
    # Handling different input VAF type.
    params_combinations, VAF_list = determine_VAF_smk('combined', args, GEs, tumor_fractions, error_rates)
    print(f"{len(VAF_list)} set of VAFs (samples) &")
    ### TODO: parameter combination --> VAF list and Amount of SNV should not be multiplied.
    print(f"{len(params_combinations)} combination parameters to simulate per VAF (sample).")

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

    # Store summary:
    pickle_path = args.output_pickle_path.rsplit("/",1)[0]
    #csv_path = args.output_csv_path.rsplit("/",1)[0]
    # print("create dir:", pickle_path)
    #check_path(csv_path)
    check_path(pickle_path)
    store_df_no_filename(df_detection_per_id, args.output_pickle_path, df_types=['pickle'])
    print(f"pickle file written to {args.output_pickle_path}")
    #store_df_no_filename(df_detection_per_id, args.output_csv_path, df_types=['csv'])
    # print(f"Output stored at:{args.output_folder}/{filename}.pickle.gz, {args.output_folder}/{filename}.csv")
    #logging.info(f"Output stored at:{f.output_folder}/{f.filename}.pickle.gz, {f.output_folder}/{f.filename}.csv")
