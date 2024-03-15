#/bin/bash
# This file contains the logic of steps to perform. Actual runs were parallelized during hpc computing per chromosome.

# Input BAM file
input_bam=$1
final_output_bam=$2
out_folder=$3
# Index input bam
samtools index "$input_bam"
# Use samtools to list all contigs in the input BAM file
contig_list=$(samtools view -H "$input_bam" | grep -E '^@SQ' | cut -f 2 -d ':' | cut -f 2 -d '@')
# Select target contigs:
contig_list=$("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
# Loop through the contigs and process each one separately
for contig in $contig_list; do
    # Create a temporary BAM file for the current contig
    temp_bam="${out_folder}/temp_${contig}.bam"
    # Extract reads for the current contig using samtools view
    samtools view -b "$input_bam" "$contig" > "$temp_bam"
    # Execute the Python script on the temporary BAM file
    python first_align_nanorcs.py -i "$temp_bam" -o "${out_folder}/output_${contig}.bam"
    # Parallel example with slurm
    #sbatch -J ${contig}_first_alignment -o ./${contig}_first_alignment.slurm.out -e ./${contig}_first_alignment.slurm.err --time=24:00:00 --mem=128G 01_get_first_alignment_wrapper.sh ${temp_bam} ${out_folder}/output_${contig}.bam
done


# Merge all the output BAM files into a single final BAM file
ls -1 ${out_folder}/output_*.bam > ${input_bam}_merge_list.txt
samtools merge -@ 32 -O BAM -b ${input_bam}_merge_list.txt -o $final_output_bam
samtools sort -@ 32 ${final_output_bam} > ${final_output_bam}.sorted.bam
samtools index ${final_output_bam}.sorted.bam

# Optionally, clean up temporary files after processing
# rm "${out_folder}/$temp_bam"
# rm "${out_folder}/output_*.bam"
# rm "${input_bam}_merge_list.txt"
