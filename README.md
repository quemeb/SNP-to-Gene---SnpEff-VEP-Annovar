# SNP to Gene - SnpEff-VEP

### Overview
This repository contains tools for mapping SNPs (Single Nucleotide Polymorphisms) from Genome Build 37 to corresponding genes in the Ensembl and RefSeq databases using SnpEff, VEP, and Annovar. Additionally, it retrieves enhancer-linked genes from Peregrine.

The project provides two approaches for users:
1. **R Shiny Application**: An interactive web-based app to upload SNP data and download results.
2. **R Script**: A fully script-based approach for advanced users who want to handle everything within R.

### Contents
- **app_updated.R**: The Shiny application file that runs a web-based interface for users to upload files and download processed SNP annotations.
- **geneIDs_AnnoQ_updated.R**: The backend script that processes SNPs through the annotation pipeline.
- **SNPs_to_Genes_Script.R**: A standalone R script for users who prefer to handle everything directly in R without using the Shiny app.

### Try the Shiny App Online
You can test the Shiny app live without needing to set up anything locally by visiting the following link:
[https://quemeb.shinyapps.io/SNP-to-genes/](https://quemeb.shinyapps.io/SNP-to-genes/)

#### Features of the Shiny App:
- **Upload CSV Files**: You can upload SNP data in CSV format with required columns (`chr`, `start`, `end`) from Genome Build 37.
- **Download Annotated Results**: Once the file is processed, you can download the SNP-to-gene mapping in CSV or plain text format.

To run the Shiny app locally:
1. Clone this repository.
2. Run the `app_updated.R` file, which will launch the interactive interface.
3. The Shiny app uses `geneIDs_AnnoQ_updated.R` as the backend to process the annotations.

### For Hardcore Programmers: Use the R Script
If you prefer to avoid using the Shiny app and would rather run the entire process within R, you can use the **SNPs_to_Genes_Script.R** script.

#### How to Use the R Script:
1. **Prepare Your Input File**: The input file should be in CSV format and must include the following columns:
    - `chr`: Chromosome number
    - `start`: SNP start position
    - `end`: SNP end position

2. **Edit the Script**:
    - Open the `SNPs_to_Genes_Script.R` file and replace the following lines with the paths to your input and output files:
    ```r
    input_file <- "path_to_your_input_file.csv"  # <-- Replace this with your input file path
    output_file <- "path_to_output_file.csv"  # <-- Replace this with your desired output file path
    ```

3. **Run the Script**:
    - The script will read the input file, process the SNPs using the annotations pipeline (Ensembl, RefSeq, Peregrine), and output the annotated results into a new CSV file.

### Dependencies
Both the Shiny app and the R script require the following R packages:
- `httr`
- `jsonlite`

You can install the required packages by running:
```r
install.packages(c("httr", "jsonlite"))
