import os

output_file = snakemake.output.summary
input_files = snakemake.input.stats_files

with open(output_file, "w") as out:
    # Write the header
    header = "Sample\tTotal_Mapped_fragments\tPeak_Count\tReads_In_Peaks\tFRiP_Score\n"
    out.write(header)
    
    # Iterate through all sample-specific QC files
    for filepath in input_files:
        if os.path.exists(filepath):
            with open(filepath, "r") as infile:
                # Simply copy the line from the individual file
                out.write(infile.read())
