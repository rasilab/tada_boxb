#!/usr/bin/awk -f

# Usage: awk -v barcode_len=10 -v output_dir="output_dir" -v read2="read2.fastq" -f split_by_barcode_with_filter.awk barcodes.csv read1.fastq

BEGIN {
    FS = ",";  # Set the field separator to comma
}

NR == FNR && NR > 1 {  # Skip the header line in the barcode list
    barcode_filter[$2] = 1;  # Mark the barcode as valid
    # initialize these empty files
    file1 = output_dir "/" barcode "_R1.fastq";  # File for Read 1
    file2 = output_dir "/" barcode "_R2.fastq";  # File for Read 2
    print "" > file1;
    print "" > file2;
    next
}

NR > FNR {
    # Process Read 1 (first FASTQ file)
    if (FNR % 4 == 1) { 
        header1 = $0;  # Read 1 sequence identifier
    } 
    else if (FNR % 4 == 2) { 
        seq1 = $0;  # Read 1 sequence
        barcode = substr(seq1, 1, barcode_len);  # Extract barcode
    } 
    else if (FNR % 4 == 3) { 
        plus1 = $0;  # Read 1 plus line
    } 
    else if (FNR % 4 == 0) { 
        qual1 = $0;  # Read 1 quality scores

        # Process corresponding Read 2 lines
        getline header2 < read2;  # Read 2 sequence identifier
        getline seq2 < read2;     # Read 2 sequence
        getline plus2 < read2;    # Read 2 plus line
        getline qual2 < read2;    # Read 2 quality scores

        # Check if the barcode is in the barcode filter
        if (barcode in barcode_filter) {
            # Store the entire FASTQ entry for Read 1 and Read 2 in memory by barcode
            read1_buffer[barcode] = read1_buffer[barcode] header1 "\n" seq1 "\n" plus1 "\n" qual1 "\n"
            read2_buffer[barcode] = read2_buffer[barcode] header2 "\n" seq2 "\n" plus2 "\n" qual2 "\n"
        }
    }
    if (NR % 400000 == 0) {
        for (barcode in barcode_filter) {
            # Output files based on the barcode
            file1 = output_dir "/" barcode "_R1.fastq";  # File for Read 1
            file2 = output_dir "/" barcode "_R2.fastq";  # File for Read 2
            print read1_buffer[barcode] >> file1;  # Write Read 1 buffer to file
            print read2_buffer[barcode] >> file2;  # Write Read 2 buffer to file
            delete read1_buffer[barcode];  # Clear Read 1 buffer
            delete read2_buffer[barcode];  # Clear Read 2 buffer
        }
    }
}

END {
    for (barcode in barcode_filter) {
        # Output files based on the barcode
        file1 = output_dir "/" barcode "_R1.fastq";  # File for Read 1
        file2 = output_dir "/" barcode "_R2.fastq";  # File for Read 2
        print read1_buffer[barcode] >> file1;  # Write Read 1 buffer to file
        print read2_buffer[barcode] >> file2;  # Write Read 2 buffer to file
    }
}