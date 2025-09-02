# If you haven't installed it yet:
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(biomaRt)

# List all available Marts to see the current name (often "ensembl")
listMarts()

# Connect to the main Ensembl Mart
ensembl <- useMart("ensembl")

# List all datasets in the Mart to find the human one (e.g., "hsapiens_gene_ensembl")
listDatasets(ensembl)

# Select the human dataset
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Read the Excel file. The 'skip=1' might be necessary if the first row is the header you showed.
# The exact command depends on how you read it in. Here's an example using readxl:
library(readxl)
genes_df <- read_excel("genes1.xlsx", skip = 1) # Use 'skip=1' to skip the initial header line

# Extract the first column, which contains the transcript IDs
transcript_ids <- genes_df[[1]] # Or use genes_df$`#name` if that's the column name

# (Optional but recommended) Remove the version number for better matching
# biomart can often handle versions, but removing them can sometimes improve match rates.
transcript_ids_no_version <- sub("\\..*", "", transcript_ids) # Removes everything after the dot

# Define what you want to get back (the attributes)
attributes_to_get <- c(
  "ensembl_transcript_id", 
  "ensembl_gene_id",
  "external_gene_name",       # The common gene symbol (e.g., EGFR, BRAF)
  "description",              # A brief description of the gene
  "gene_biotype",             # e.g., "protein_coding", "lncRNA"
  "chromosome_name",          # To double-check it's chr7
  "start_position", 
  "end_position"
)

# Define what you are using as input (the filter)
filters <- "ensembl_transcript_id"

# Now run the query
gene_annotations <- getBM(
  attributes = attributes_to_get,
  filters = filters,
  values = transcript_ids_no_version, # Use the list without version numbers
  mart = ensembl
)

# View the results
head(gene_annotations)
View(gene_annotations) # Opens a spreadsheet-like view