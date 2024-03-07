library(fasterqParseR)

## SRRs in your fastq file
srr_ids <- c("SRR12345678")

working_dir = "~/Project_X/"
input_dir = "~/Project_X/fasterq_output/"
outdir="~/Project_X/read_lengths/"


# Initialize an empty list to store results
assigned_files_list <- list()

for (srr_id in srr_ids) {
  
  # Run your functions for each SRR ID
  assigned_files <- assignSRAreads(working_dir = working_dir, 
                                   input_dir = input_dir, 
                                   outdir = outdir, 
                                   parallel = TRUE)
  
  # Store the result in the list with the SRR ID as the key
  assigned_files_list[[srr_id]] <- assigned_files
}

# to save the csv so it isnt written over once renaming occurs:
write.csv(assigned_files[,-1, with = FALSE], paste0(outdir, "assigned_SRAreads_final.csv"), row.names = FALSE)

# to rename everything
renameAll(assigned_SRA= assigned_files, input_dir=input_dir, format="cellranger")
