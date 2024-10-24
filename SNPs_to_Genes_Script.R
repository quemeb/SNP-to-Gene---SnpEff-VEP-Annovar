# User-friendly setup for input and output files
input_file <- "path_to_your_input_file.csv"  # <-- Replace this with the path to your input CSV file
output_file <- "path_to_output_file.csv"  # <-- Replace this with the desired output path

# Load required packages
if (!require("httr", quietly = TRUE))
  install.packages("httr")

if (!require("jsonlite", quietly = TRUE))
  install.packages("jsonlite")

library(httr)
library(jsonlite)

# Define the Base URL
Annotations_URL <- "http://annoq.org/api-v2/graphql"

# Convert a list of fields into a GraphQL query format
create_annotations_query_string <- function(annotations) {
  paste(annotations, collapse = "\n")
}

# Perform a GraphQL query
perform_graphql_query <- function(query) {
  response <- POST(
    Annotations_URL, 
    content_type_json(), 
    body = list(query = query),
    encode = "json"
  )
  stop_for_status(response)
  content(response, "text", encoding = "UTF-8")
}

# Desired fields
annotations_to_retrieve <- c(
  "rs_dbSNP151",
  "ANNOVAR_ensembl_Gene_ID",
  "ANNOVAR_ensembl_Closest_gene",
  "SnpEff_ensembl_Gene_ID",
  "VEP_ensembl_Gene_ID",
  "ANNOVAR_refseq_Gene_ID",
  "ANNOVAR_refseq_Closest_gene",
  "SnpEff_refseq_Gene_ID",
  "VEP_refseq_Gene_Name",
  "enhancer_linked_genes"
)

# Get SNP by chromosome, start, and end
snpQuery <- function(chr, start, end, annotations_to_retrieve) {
  annotations_query_string <- create_annotations_query_string(annotations_to_retrieve)
  query <- sprintf('
  query {
    get_SNPs_by_chromosome(chr: "%s", start: %d, end: %d, query_type_option: SNPS) {
      snps {
        %s
      }
    }
  }', chr, start, end, annotations_query_string)
  
  response_content <- perform_graphql_query(query)
  data <- fromJSON(response_content, flatten = TRUE)
  data$data$get_SNPs_by_chromosome$snps
}

# Function to process and extract gene annotations for a single SNP position
process_annotations <- function(chr, start, end) {
  query_result <- snpQuery(chr, start, end, annotations_to_retrieve)
  
  if (!is.null(query_result) && length(query_result) > 0) {
    # Process Ensembl genes
    ensembl_genes <- unlist(query_result[grep("ensembl", names(query_result))])
    ensembl_genes <- ensembl_genes[!is.na(ensembl_genes)]
    ensembl_genes <- unlist(strsplit(as.character(ensembl_genes), "\\|"))
    ensembl_genes <- unlist(strsplit(ensembl_genes, "[-:]"))
    ensembl_genes <- ensembl_genes[!is.na(ensembl_genes)]
    ensembl_genes <- ensembl_genes[grepl("^ENSG", ensembl_genes)]
    ensembl_genes <- unique(ensembl_genes)
    
    # Process RefSeq genes
    refseq_genes <- unlist(query_result[grep("refseq", names(query_result))])
    refseq_genes <- refseq_genes[!is.na(refseq_genes)]
    refseq_genes <- unlist(strsplit(as.character(refseq_genes), "\\|"))
    refseq_genes <- unique(refseq_genes)
    
    # Process enhancers
    enhancers <- query_result$enhancer_linked_genes
    enhancers <- enhancers[!is.na(enhancers)]
    enhancers <- unlist(strsplit(as.character(enhancers), ";"))
    enhancers <- unique(enhancers)
    
    # Process rs_dbSNP151
    rs_dbSNP151 <- query_result$rs_dbSNP151
    rs_dbSNP151 <- unique(rs_dbSNP151)
    
    # Combine and format output
    ensembl_genes <- paste(ensembl_genes[!is.na(ensembl_genes)], collapse = ",")
    refseq_genes <- paste(refseq_genes[!is.na(refseq_genes)], collapse = ",")
    enhancers <- paste(enhancers[!is.na(enhancers)], collapse = ",")
    rs_dbSNP151 <- paste(rs_dbSNP151[!is.na(rs_dbSNP151)], collapse = ",")
    
    return(data.frame(
      chr = chr,
      start = start,
      end = end,
      Ensembl_Genes = ensembl_genes,
      RefSeq_Genes = refseq_genes,
      Enhancers = enhancers,
      rs_dbSNP151 = rs_dbSNP151,
      stringsAsFactors = FALSE
    ))
  } else {
    return(data.frame(
      chr = chr,
      start = start,
      end = end,
      Ensembl_Genes = "",
      RefSeq_Genes = "",
      Enhancers = "",
      rs_dbSNP151 = "",
      stringsAsFactors = FALSE
    ))
  }
}

# Function to match gene annotations back to the original df by SNP position
match_annotations_to_df <- function(df, annotations_df) {
  # Merge the annotations with the original df on SNP position
  merged_df <- merge(df, annotations_df, by = c("chr", "start", "end"), all.x = TRUE)
  return(merged_df)
}

# 1. Read the input CSV file
df <- read.csv(input_file)

# 2. Process the SNPs and generate the annotations
annotations_df <- data.frame(
  chr = character(),
  start = integer(),
  end = integer(),
  Ensembl_Genes = character(),
  RefSeq_Genes = character(),
  Enhancers = character(),
  rs_dbSNP151 = character(),
  stringsAsFactors = FALSE
)

# Loop through each SNP in the input CSV
for (i in 1:nrow(df)) {
  chr <- df$chr[i]
  start <- df$start[i]
  end <- df$end[i]
  annotation <- process_annotations(chr, start, end)
  annotations_df <- rbind(annotations_df, annotation)
}

# 3. Match annotations back to the original dataframe
final_df <- match_annotations_to_df(df, annotations_df)

# 4. Save the output to a CSV file
write.csv(final_df, output_file, row.names = FALSE)

