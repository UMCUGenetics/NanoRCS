import argparse
import pysam
from collections import defaultdict, Counter
 
def main(input_bam_file, output_bam_file):
    # Open the input BAM file for reading and the output BAM file for writing
    with pysam.AlignmentFile(input_bam_file, "rb") as input_bam:
        # Create a dictionary to count read occurrences
        read_counts = Counter()
        read_dict = defaultdict(list)
        # Iterate through the reads in the input BAM file
        for read in input_bam:
            read_name = read.query_name
            # Count the occurrences of each read
            read_counts[read_name] += 1
            read_dict[read_name].append(read)
    with open(output_bam_file+'_read_count.dict', 'w+') as f:
        f.write(str(read_counts))
    assessed_read_names = []
    with pysam.AlignmentFile(input_bam_file, "rb") as input_bam,  \
            pysam.AlignmentFile(output_bam_file, "wb", header=input_bam.header) as output_bam:
        for read in input_bam:
            read_name = read.query_name
            if read_name in assessed_read_names:
                continue
            else:
                assessed_read_names.append(read_name)
                if read_counts[read_name] >= 3:
                    read_list = read_dict[read_name]
                    clip_distances = []
                    for target_read in read_list:
                    # For example, check if the CIGAR string starts with 'H' or 'S'
                        cigar_string = target_read.cigarstring
                        clipping_distance = cigar_string.split('H')[0].split('S')[0]
                        try:
                            clipping_distance = int(clipping_distance)
                        except Exception as e:
                            print(cigar_string)
                            clipping_distance = 100000
                        clip_distances.append(clipping_distance)
                    location_in_list_first_alignment = clip_distances.index(min(clip_distances))
                    first_alignment = read_list[location_in_list_first_alignment]
                    # Write the read to the output BAM file
                    output_bam.write(first_alignment)

if __name__ == "__main__":
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Process BAM file with specified criteria")
    # Add arguments for input and output BAM file paths
    parser.add_argument("-i", "--input_bam_file", help="Path to the input BAM file")
    parser.add_argument("-o", "--output_bam_file", help="Path to the output BAM file")
    # Parse the command-line arguments
    args = parser.parse_args()
    # Call the main function with the provided input and output BAM file paths
    main(args.input_bam_file, args.output_bam_file)
