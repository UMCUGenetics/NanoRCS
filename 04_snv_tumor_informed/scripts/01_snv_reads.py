import argparse
import os
import re
import shlex
import subprocess
import time
import warnings
from pathlib import Path
import math
import cyvcf2
import numpy as np
import pandas as pd
import pysam
from cyvcf2 import VCF, Writer

with warnings.catch_warnings():
    warnings.simplefilter("ignore")


def prepare_vcf(vcf_path, input_vcf):
    # input_vcf.header
    input_vcf = VCF(vcf_path)
    input_vcf.add_info_to_header(
        {
            "ID": "gene",
            "Description": "overlapping gene",
            "Type": "Character",
            "Number": "1",
        }
    )
    w = Writer(vcf_path, input_vcf)
    return w


def write_vcf_record(w, v):
    w.write_record(v)


def prob_to_phred(prob: float):
    """
    Convert probability of base call being correct into phred score
    Values are clipped to stay within 0 to 60 phred range
    Args:
        prob  (float): probability of base call being correct
    Returns:
        phred_score (byte)
    """
    return np.rint(-10 * np.log10(np.clip(1 - prob, 1 - 0.999999, 0.999999))).astype(
        "B"
    )


def phred_to_prob(phred):
    """Convert a phred score (ASCII) or integer to a numeric probability
    Args:
        phred (str/int) : score to convert
    returns:
        probability(float)
    """

    try:
        if isinstance(phred, int):
            return math.pow(10, -(phred) / 10)
        return math.pow(10, -(ord(phred) - 33) / 10)
    except ValueError:
        return 1


def examine_path(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return 0


def check_vcf_format(in_vcf):
    if not os.path.exists(in_vcf):
        print(f"VCF file does not exist: {in_vcf}")
        exit(1)
    if in_vcf.endswith("vcf"):
        print("VCF file is not gzip and indexed. Execute the following")
        print(
            f"bgzip -c {in_vcf} > {in_vcf}.gz;",
        )
        print(f"tabix -p vcf {in_vcf}.gz")
        exit(1)
    if not os.path.exists(in_vcf+'.tbi'):
        print("VCF file is not indexed. Execute the following:")
        print(f"tabix -p vcf {in_vcf}.gz")
        exit(1)


def pileup_cleanup(pileup_record):
    """
    exceptions includes:
    ^ begin of a sequence
    $ after end of a sequence
    X-NNN (N amount depends on the following deletion)
    * deletion
    >>>pileup_cleanup("^ATT")
    "ATT"
    >>>pileup_cleanup("A-1N*T")
    "A*T"
    >>>pileup_cleanup("A-2N**")
    "A"
    >>>pileup_cleanup("*T-1N*")
    "T"
    >>>pileup_cleanup("T$")
    "T"
    >>>pileup_cleanup("G-3NNN***")
    "G"
    """
    pileup_record = re.sub("[0-9]", "", pileup_record)
    pileup_record = pileup_record.replace("-", "")
    pileup_record = pileup_record.replace("N", "")
    pileup_record = pileup_record.replace("^", "")
    pileup_record = pileup_record.replace("$", "")
    pileup_record = pileup_record.replace("*", "")

    return pileup_record


# Get only the reads that overlap with somatic sSNV. Can skip.
def get_overlapping_vcf_reads_bedtools(
    somatic_vcf, input_bam_path, output_overlap_path
):
    pass
    subprocess.call(f"touch {output_overlap_path}.touchtest")
    cmd = shlex.split(
        f"bedtools intersect -wa -a {input_bam_path}  -b {somatic_vcf} > {output_overlap_path}"
    )
    subprocess.call(cmd)


class Alignments:
    def __init__(self, bam, out_ref_bam, out_alt_bam, out_exception_bam):
        self.alignments = []
        # Use this to assign alignments to out bam is faster
        self.ordered_alignments_tags = []
        self.phased_tags = []  # -1, 0, 1
        self.template_bam = bam
        self.MUT_alignments = []
        self.WT_alignments = []
        self.exception_alignments = []
        self.not_phased_ref_alignments = []
        self.ref_out_bam = out_ref_bam
        self.alt_out_bam = out_alt_bam
        self.exception_bam = out_exception_bam
        self.associated_gsnv = []
        self.associated_gsnv_labels = []

    def __len__(self):
        return len(self.alignments)

    # DA: allele, DS: siteCoordinate, DB: allele for gSNV , DC: siteCoordinate for gSNV
    def write_records(self):
        # Write wild type alignments by looking at ordered alignment tags
        array = np.array(self.ordered_alignments_tags)
        indices = list(np.where(array == 0)[0])
        self.WT_alignments = map(self.alignments.__getitem__, indices)
        # Get MUT
        array = np.array(self.ordered_alignments_tags)
        indices = list(np.where(array > 0)[0])
        self.MUT_alignments = map(self.alignments.__getitem__, indices)

        array = np.array(self.ordered_alignments_tags)
        indices = list(np.where(array == -2)[0])
        self.exception_alignments = map(self.alignments.__getitem__, indices)

        for read in self.WT_alignments:
            self.ref_out_bam.write(read)
        for read in self.MUT_alignments:
            self.alt_out_bam.write(read)
        for read in self.exception_alignments:
            self.exception_bam.write(read)


class SomaticNucleotideVariants:
    def __init__(self, snv, criteria=None, phasing=False, cfdna_length_range=350):
        # inherit gsnv object from cyvcf2
        if criteria is None:
            criteria = {"FILTER": None, "QUAL": 30, "LEN": 0}
        self.snv = snv
        self.CHROM = snv.CHROM
        self.start = snv.start
        self.end = snv.end
        self.ALT = snv.ALT
        self.REF = snv.REF
        # Boolean. If we want to phase it or not.
        self.phasing = phasing
        # add customized features
        self.ssnv_length = snv.end - snv.start
        self.snv_type = None
        # A list containing only alignments that contains gsnv. Updated over time.
        self.pileup_alignments = []
        self.query_sequences = []
        self.query_qualities = []
        self.ssnv_ref_alleles = []
        self.ssnv_alt_alleles = []
        self.criteria = criteria
        self.pass_filter = False
        self.phasing_info = None  # Get information if it should be ALT-REF or ALT-ALT for the tumor allele. From vcf?
        self.gsnv_ref_alleles = []
        self.gsnv_alt_alleles = []
        self.remaining_positions = None
        self.align_to_ref = []
        self.align_to_alt = {}
        self.related_gsnv = []
        self.filter()
        self.cfdna_length_range = cfdna_length_range
        self.related_gsnv = []

    def get_gsnvs(self):
        self.related_gsnv = []

    def snv_type(self):
        if self.ssnv_length > 1:
            self.snv_type = "indel"
        else:
            self.snv_type = "gsnv"

    def filter(self):
        self.pass_filter = True
        # # LEN is used as a special criteria that QUAL is not used.
        # if self.criteria["LEN"] == 0:
        #     # PASS is rendered as None. Anything else are not None.
        #     if self.snv.FILTER == self.criteria["FILTER"]:
        #         self.pass_filter = True
        #
        # elif self.snv.FILTER == self.criteria["FILTER"]:
        #     if self.snv.QUAL is None:
        #         self.pass_filter = False
        #
        # elif self.snv.FILTER == self.criteria["FILTER"]:
        #     if self.snv.QUAL >= self.criteria["QUAL"]:
        #         self.pass_filter = True
        # else:
        #     self.pass_filter = False


# For snv
def get_alignments_and_query_per_snv(
    snv, bam, out_alignment_records, snv_db, tumor_index, min_qual=40,
):

    pileup_columns = bam.pileup(
        snv.CHROM, snv.start, snv.end, min_mapping_quality=0, min_base_quality=min_qual
    )
    # all columns of reads that overlap with this particular position.

    snv.remaining_positions = snv.ssnv_length
    # if snv.ssnv_length > 1:
    #     pass
    # is_sv, is_mnp, is_indel, is_deletion, is_transition(subtype ts, not transition (trasversion = subtype tv)
    # if snv.is_snp:
    #     pass

    # If no read pileup at this location, the for loop would finish.
    for pileup_column in pileup_columns:
        if pileup_column.reference_pos == snv.start:
            snv.pileup_alignments = [
                pileup_column.pileups[i].alignment
                for i in range(len(pileup_column.pileups))
            ]

            query_sequences = pileup_column.get_query_sequences(add_indels=True)
            snv.query_qualities = pileup_column.get_query_qualities()
            snv.query_sequences = [v.upper() for v in query_sequences]
            snv.remaining_positions -= 1
            if snv.remaining_positions == 0:
                break
                # to collect more positions
        elif (pileup_column.reference_pos < snv.end) & (
            pileup_column.reference_pos > snv.start
        ):
            # If there is a deletion at this position, the length of the query_sequences would be shorter.
            # It is not possible to know which base should stitch to which. (Same problem as Inez faced)
            query_sequences_next = pileup_column.get_query_sequences(add_indels=True)
            query_sequences_next = [v.upper() for v in query_sequences_next]

            # Edge case where the alignment does not cover the first base of the VCF record.
            if len(snv.query_sequences) == 0:
                # First position is not covered. The covered part cannot be used to determine which allele it is.
                return 0
            # if all reads covering the indel positions are not prematurely short, then process normally
            elif len(query_sequences_next) == len(snv.query_sequences):
                snv.query_sequences = [
                    a + b for (a, b) in zip(snv.query_sequences, query_sequences_next)
                ]

            # Edge cases where one of the bam record is much shorter, and did not cover the complete indel length
            else:
                # kick out the alignments that are shorter than the gsnv/indel length
                # because we cannot determine the allele there (enter group "reads_not_aligned")
                alignments = [
                    pileup_column.pileups[i].alignment
                    for i in range(len(pileup_column.pileups))
                ]
                # Find the alignment that is not in the first position. Any alignment that either start later,
                # or ends prematurely will not be in the final list
                alignments_to_kick_out = set(snv.pileup_alignments) - set(alignments)
                for alignment in alignments_to_kick_out:
                    index_to_kick_out = snv.pileup_alignments.index(alignment)
                    snv.query_sequences.pop(index_to_kick_out)
                    snv.pileup_alignments.pop(index_to_kick_out)
                # Now the query_sequences only contains the bases of the remaining alignments
                snv.query_sequences = [
                    a + b for (a, b) in zip(snv.query_sequences, query_sequences_next)
                ]
            snv.remaining_positions -= 1
            if snv.remaining_positions == 0:
                break
        else:
            pass
    # If bam overlaps with SNV:
    if snv.remaining_positions == 0:
        ref_alt_id = np.full((len(snv.query_sequences)), -1)
        for i, (query_sequence, query_sequence_quality) in enumerate(
            zip(snv.query_sequences, snv.query_qualities)
        ):
            query_sequence = pileup_cleanup(query_sequence)
            # print(snv.REF, snv.ALT, query_sequence)
            if (
                (query_sequence not in snv.ALT)
                and (query_sequence != snv.REF)
                and (query_sequence != "N")
                and (query_sequence != "")
                and (len(query_sequence) == 1)
            ):
                # Dep ref_alt_id
                ref_alt_id[i] = -2
                ## Error allele
                snv_db.snv_db_list.append(
                    [
                        f"{snv.snv.CHROM}_{snv.snv.start}_{snv.snv.end}",
                        snv.snv.CHROM,
                        snv.snv.start,
                        snv.snv.end,
                        snv.snv.REF,
                        snv.snv.ALT[0],
                        snv.snv.QUAL,
                        snv.snv.FILTER,
                        snv.snv.gt_alt_freqs[tumor_index],
                        query_sequence,
                        query_sequence_quality,
                        snv.pileup_alignments[i].mapq,
                        "ERR",
                    ]
                )

                # Does not align to any known allele
                # snv.pileup_alignments[i].set_tag(
                #     "DP", f"{snv.CHROM}_{snv.start}_{snv.end}"
                # )
                # snv.pileup_alignments[i].set_tag("DS", None)
                # snv.pileup_alignments[i].set_tag("DA", -1)

                # Updated alignment is added to a new list.
                if snv.pileup_alignments[i] not in out_alignment_records.alignments:
                    out_alignment_records.alignments.append(snv.pileup_alignments[i])
                    out_alignment_records.ordered_alignments_tags.append(-2)

                # For debugging:
                # print("NOT MATCHED\n", f"{snv.CHROM}:{snv.start}-{snv.end}\n",  snv.pileup_alignments[i])
                # Write all unmapped sequence to file. Very valuable.

            elif query_sequence == snv.REF:
                # Dep ref_alt_id
                ref_alt_id[i] = 0
                snv_db.snv_db_list.append(
                    [
                        f"{snv.snv.CHROM}_{snv.snv.start}_{snv.snv.end}",
                        snv.snv.CHROM,
                        snv.snv.start,
                        snv.snv.end,
                        snv.snv.REF,
                        snv.snv.ALT[0],
                        snv.snv.QUAL,
                        snv.snv.FILTER,
                        snv.snv.gt_alt_freqs[tumor_index],
                        query_sequence,
                        query_sequence_quality,
                        snv.pileup_alignments[i].mapq,
                        "REF",
                    ]
                )

                # Add to ordered list

                # Update tag to a bam  alignmentRead
                # snv.pileup_alignments[i].set_tag(
                #     "DP", f"{snv.CHROM}_{snv.start}_{snv.end}"
                # )
                # snv.pileup_alignments[i].set_tag("DS", f"{snv.REF}")
                # snv.pileup_alignments[i].set_tag("DA", 0)
                # Updated alignment is added to a new list.
                if snv.pileup_alignments[i] not in out_alignment_records.alignments:
                    out_alignment_records.ordered_alignments_tags.append(0)
                    out_alignment_records.alignments.append(snv.pileup_alignments[i])
            for j in range(len(snv.ALT)):
                allele_num = j + 1  # start counting by 1
                if query_sequence == snv.ALT[j]:
                    ref_alt_id[i] = allele_num
                    snv_db.snv_db_list.append(
                        [
                            f"{snv.snv.CHROM}_{snv.snv.start}_{snv.snv.end}",
                            snv.snv.CHROM,
                            snv.snv.start,
                            snv.snv.end,
                            snv.snv.REF,
                            snv.snv.ALT[0],
                            snv.snv.QUAL,
                            snv.snv.FILTER,
                            snv.snv.gt_alt_freqs[tumor_index],
                            query_sequence,
                            query_sequence_quality,
                            snv.pileup_alignments[i].mapq,
                            "ALT",
                        ]
                    )
                    # snv.pileup_alignments[i].set_tag(
                    #     "DP", f"{snv.CHROM}_{snv.start}_{snv.end}"
                    # )
                    # snv.pileup_alignments[i].set_tag("DS", f"{snv.ALT[j]}")
                    # snv.pileup_alignments[i].set_tag("DA", allele_num)
                    if snv.pileup_alignments[i] not in out_alignment_records.alignments:
                        out_alignment_records.alignments.append(
                            snv.pileup_alignments[i]
                        )
                        out_alignment_records.ordered_alignments_tags.append(allele_num)

                    else:
                        pass
            else:
                pass
                # print("This SNV has no pileup bam read.")
                # or the position is an "N".

                # For debugging:
                # print(f"ALT: {snv.CHROM}:{snv.start}-{snv.end}\n", snv.pileup_alignments[i])

        # print("Allele info:", ref_alt_id)

    else:
        pass
        # There are remaining positions


class MatchingSnvDB:
    def __init__(self, out_path):
        self.out_path = out_path
        self.snv_db = None
        self.snv_db_list = []
        # self.snv = None
        # self.query_sequence = None
        # self.snv_overlap_type = None

    def finalize_db(self):
        self.snv_db = pd.DataFrame(
            self.snv_db_list,
            columns=[
                "variant_id",
                "CHROM",
                "start",
                "end",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "ALT_freq",
                "query_base",
                "query_base_qual",
                "query_map_qual",
                "snv_overlap_type",
            ],
        )

    # def perge_variant(self):
    #     self.snv = None
    #     self.query_sequence = None
    #     self.snv_overlap_type = None

    def store_db(self):
        self.snv_db.to_pickle(self.out_path)


def find_ssnv(
    bam_path,
    somatic_vcf_path,
    germline_vcf_path,
    out_ref_bam_path,
    out_alt_bam_path,
    out_ERR_bam_path,
    out_variants,
    out_all_variants,
    phasing=False,
    min_qual=40,
    record_all_variants=False,
):
    check_vcf_format(somatic_vcf_path)
    somatic_vcf = cyvcf2.VCF(somatic_vcf_path)
    if phasing is True:
        check_vcf_format(germline_vcf_path)
        germline_vcf = cyvcf2.VCF(germline_vcf_path)
    bam = pysam.AlignmentFile(bam_path)
    pysam.index(bam_path)
    out_ref_bam = pysam.AlignmentFile(out_ref_bam_path, "wb", template=bam)
    out_alt_bam = pysam.AlignmentFile(out_alt_bam_path, "wb", template=bam)
    out_ERR_bam = pysam.AlignmentFile(out_ERR_bam_path, "wb", template=bam)
    count_filtered_ssnv = 0
    count_bam_read = 0
    allele_ref_counts = 0
    allele_alt_counts = 0
    allele_not_matched = 0
    allele_debug_should_be_zero = 0

    # Create DB for writing, store info to write at ssnv object, write after for loop.
    snv_db = MatchingSnvDB(out_variants)
    if record_all_variants:
        snv_db_all = MatchingSnvDB(out_all_variants)
    # Check if it is a multi sample vcf.If it is, find the index of tumor sample
    if len(somatic_vcf.samples) != 1:
        print("The vcf contains multiple samples. Check index of the samples of interest. Here is a list of samples:")
        print(somatic_vcf.samples)
        for ind, sample in enumerate(somatic_vcf.samples):
            if ('-T' in sample ) or ('tumor' in sample):
                print(f"Automatic detecting...'-T' or 'tumor' found. "
                      f"Is '{sample}' the tumor sample? \n"
                      f"Continue with this assumption. Please kill the program if this is not correct.")
                tumor_index = ind
            else:
                pass
        if tumor_index is None:
            print("Could not detect which sample is tumor sample. Exit the program.")
            exit(1)

    else:
        print("The vcf contains one samples. Check if the sample is expected:")
        print(f"Continue with this assumption. Please kill the program if this is not correct.")
        print(somatic_vcf.samples)
        tumor_index = 0
    count = 0
    for ssnv_raw in somatic_vcf:
        count += 1
        if count % 100000 == 0:
            print(count/1000, 'k')
        # ssnv.snv.QUAL is empty for GCT02

        ssnv = SomaticNucleotideVariants(
            ssnv_raw,
            criteria={"FILTER": None, "QUAL": 30, "LEN": 0},
            phasing=phasing,
        )
        if ssnv_raw.format('DP')[tumor_index][0] <= 30:
            continue
        if ssnv_raw.format('AD')[tumor_index][1] <= 2:
            continue
        # If a SNV has dbSNP membership == True, do not use this variant.
        if ssnv.snv.INFO.get("DB") is not None:
            continue
        out_alignment_records = Alignments(bam, out_ref_bam, out_alt_bam, out_ERR_bam)
        # if variant pass filter, action per variant. Record all variants?
        if ssnv.pass_filter and ssnv.snv.is_snp and len(ssnv.snv.ALT) == 1:
            get_alignments_and_query_per_snv(
                ssnv, bam, out_alignment_records, snv_db, tumor_index=tumor_index, min_qual=min_qual,
            )
            count_filtered_ssnv += 1
            if record_all_variants:
                snv_db_all.snv_db_list.append(
                    [
                        f"{ssnv.snv.CHROM}_{ssnv.snv.start}_{ssnv.snv.end}",
                        ssnv.snv.CHROM,
                        ssnv.snv.start,
                        ssnv.snv.end,
                        ssnv.snv.REF,
                        ssnv.snv.ALT[0],
                        ssnv.snv.QUAL,
                        ssnv.snv.FILTER,
                        ssnv.snv.gt_alt_freqs[tumor_index],
                        None,
                        None,
                        None,
                        None,
                    ]
                )

            # Per snv positions, write records
            if len(out_alignment_records) > 0:
                out_alignment_records.write_records()

                count_bam_read += len(out_alignment_records)
                arr = np.array(out_alignment_records.ordered_alignments_tags)
                allele_ref_counts += np.count_nonzero(arr == 0)
                allele_alt_counts += np.count_nonzero(arr > 0)
                allele_not_matched += np.count_nonzero(arr == -2)
                allele_debug_should_be_zero += np.count_nonzero(arr == -1)

        else:
            # The SNV did not pass filtering criteria
            continue

    # Write all snvs to db
    if record_all_variants:
        snv_db_all.finalize_db()
        snv_db_all.store_db()
    snv_db.finalize_db()
    snv_db.store_db()
    # Check each gsnv alignment record if they are in snv list, if so, add info to snv read.
    print(
        "filtered SNV:",
        count_filtered_ssnv,
        "overlapping reads:",
        count_bam_read,
        "Ref allele:",
        allele_ref_counts,
        "Alt allele:",
        allele_alt_counts,
    )
    print("allele not matched either", allele_not_matched)
    print(
        "ALT allele ratio in matching molecules (related to Tumor fraction) =Alt/(Ref+Alt) ratio:",
        round(allele_alt_counts / (allele_alt_counts + allele_ref_counts), 5),
    )
    print(
        "Error rate =Err/(Ref+Alt+err) ratio:",
        round(
            allele_not_matched
            / (allele_alt_counts + allele_ref_counts + allele_not_matched),
            5,
        ),
    )
    # print("Phred", 1 - prob_to_phred(allele_not_matched / (allele_alt_counts + allele_ref_counts + allele_not_matched)))
    bam.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b",
        "--input_bam_path",
        type=str,
        help="Input bam file of interest for molecule to tag.",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--somatic_vcf_path",
        type=str,
        help="Input path of somatic vcf file (gz + tabix).",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--germline_vcf_path",
        type=str,
        help="Input germline vcf file. If nothing is given, no phasing is done.",
        default=None,
    )
    parser.add_argument(
        "-q",
        "--min_qual",
        type=int,
        default=40,
        help="Minimal  base quality of which should be considered as detection.",
    )

    parser.add_argument(
        "-alt",
        "--out_bam_alt",
        type=str,
        help="Path ALT",
    )
    parser.add_argument(
        "-ref",
        "--out_bam_ref",
        type=str,
        help="Path REF",
    )
    parser.add_argument(
        "-other",
        "--out_bam_other",
        type=str,
        help="Path Error",
    )
    parser.add_argument(
        "-all",
        "--all_variants",
        type=Path,
        help="Path all variants",
    )
    parser.add_argument(
        "-pickle",
        "--out_pickle",
        type=str,
        help="Path Error",
    )

    parser.add_argument(
        "-o",
        "--out_bam_path",
        type=str,
        help="Path REF",
    )

    parser.add_argument(
        "-p",
        "--phasing",
        type=bool,
        default=False,
        help="If phasing should be applied.",
    )
    args = parser.parse_args()
    examine_path(args.out_bam_ref.rsplit("/",1)[0])
    # TODO: remove suffix
    out_ref_bam_path = args.out_bam_ref
    # print("out_ref_bam_path", out_ref_bam_path)
    out_variants = args.out_pickle
    out_alt_bam_path = args.out_bam_alt
    out_ERR_bam_path = args.out_bam_other
    out_all_variants = args.all_variants
    # Not implemented
    # filename = os.path.basename(args.out_bam_ref).split(".")[0]
    # out_cfdna_vcf = f"{args.out_bam_path}/{filename}_cfDNA_minQUAL{args.min_qual}.vcf"
    # TODO: remove unsorted files.
    # print("input bam:", args.input_bam_path)
    # print("input somatic vcf:", args.somatic_vcf_path)


    start_time = time.time()
    find_ssnv(
        args.input_bam_path,
        args.somatic_vcf_path,
        args.germline_vcf_path,
        out_ref_bam_path,
        out_alt_bam_path,
        out_ERR_bam_path,
        out_variants,
        out_all_variants,
        phasing=False,
        min_qual=args.min_qual,
        record_all_variants=True,
    )
    end_time = time.time()
    # Execute
    pysam.sort('-o', out_ref_bam_path+'.sorted.bam', out_ref_bam_path)
    pysam.index(out_ref_bam_path+'.sorted.bam')
    print(
        f"Program finishes successfully. Time spent: {round((end_time - start_time)/60,2)} mins."
    )
