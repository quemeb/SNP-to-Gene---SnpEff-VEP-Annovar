import pandas as pd
import requests

# Function to send request to the PANTHER API
def query_panther_api(gene_id, organism="9606"):
    base_url = "https://pantherdb.org/services/oai/pantherdb/geneinfo"
    params = {
        "geneInputList": gene_id,
        "organism": organism
    }
    
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        return response.json()
    else:
        return None

# Function to process the API response and extract pathways and accessions
def process_api_response(api_response):
    pathways = {
        "PANTHER": set(),
        "Reactome": set()
    }
    accessions = set()

    if 'mapped_genes' in api_response['search']:
        genes = api_response['search']['mapped_genes']['gene']
        
        if isinstance(genes, dict):
            genes = [genes]
        
        for gene in genes:
            accession = gene.get('accession')
            accessions.add(accession)
            
            annotation_list = gene.get('annotation_type_list', {}).get('annotation_data_type', [])
            for annotation_type in annotation_list:
                if isinstance(annotation_type, dict):
                    content_type = annotation_type.get('content')
                    annotations = annotation_type.get('annotation_list', {}).get('annotation', [])
                
                    if isinstance(annotations, dict):
                        annotations = [annotations]
                    
                    for annotation in annotations:
                        pathway_entry = f"{annotation.get('name')} ({annotation.get('id')})"
                        if content_type == "ANNOT_TYPE_ID_PANTHER_PATHWAY":
                            pathways["PANTHER"].add(pathway_entry)
                        elif content_type == "ANNOT_TYPE_ID_REACTOME_PATHWAY":
                            pathways["Reactome"].add(pathway_entry)  # Corrected this line

    return list(accessions), '; '.join(pathways["PANTHER"]), '; '.join(pathways["Reactome"])


# Function to read the CSV, process each row, and handle multiple accessions
def process_csv(file_path):
    df = pd.read_csv(file_path)
    processed_rows = []

    for idx, row in df.iterrows():
        chromosome = row['Chromosome']
        start = row['Start']
        end = row['End']
        
        # Split by comma for each column to get multiple gene IDs
        ensembl_genes = str(row['Ensembl Gene']).split(',')
        refseq_genes = str(row['RefSeq Gene']).split(',')
        enhancers = str(row['Enhancers']).split(',')
        
        # Collect all gene IDs
        gene_ids = [gene.strip() for gene_list in [ensembl_genes, refseq_genes, enhancers] for gene in gene_list if pd.notnull(gene)]

        # Dictionary to store pathways and gene IDs per accession
        accessions_dict = {}

        # Query the API for each gene ID individually
        for gene_id in gene_ids:
            api_response = query_panther_api(gene_id)
            if api_response:
                accessions, panther_pathways, reactome_pathways = process_api_response(api_response)
                
                # Add pathways and gene data grouped by unique accession
                for accession in accessions:
                    if accession not in accessions_dict:
                        accessions_dict[accession] = {
                            'Chromosome': chromosome,
                            'Start': start,
                            'End': end,
                            'Ensembl Gene': [],
                            'RefSeq Gene': [],
                            'Enhancers': [],
                            'PANTHER Pathways': panther_pathways,
                            'Reactome Pathways': reactome_pathways,
                            'Gene ID': accession
                        }
                    
                    # Append gene to appropriate list based on its type
                    if gene_id in ensembl_genes:
                        accessions_dict[accession]['Ensembl Gene'].append(gene_id)
                    elif gene_id in refseq_genes:
                        accessions_dict[accession]['RefSeq Gene'].append(gene_id)
                    elif gene_id in enhancers:
                        accessions_dict[accession]['Enhancers'].append(gene_id)
            else:
                print(f"Failed to get API response for gene ID {gene_id}")

        # Add rows for each unique accession found for the SNP
        for accession, acc_info in accessions_dict.items():
            new_row = {
                'Chromosome': acc_info['Chromosome'],
                'Start': acc_info['Start'],
                'End': acc_info['End'],
                'Ensembl Gene': '; '.join(set(acc_info['Ensembl Gene'])),
                'RefSeq Gene': '; '.join(set(acc_info['RefSeq Gene'])),
                'Enhancers': '; '.join(set(acc_info['Enhancers'])),
                'Gene ID': acc_info['Gene ID'],
                'PANTHER Pathways': acc_info['PANTHER Pathways'],
                'Reactome Pathways': acc_info['Reactome Pathways']
            }
            processed_rows.append(new_row)

    # Create a new DataFrame from processed rows and save to CSV
    result_df = pd.DataFrame(processed_rows)
    result_df.to_csv('C:/Users/bryan/Desktop/Analysis package/Updated_Table_with_Accessions.csv', index=False)
    print("Processing complete. Pathways and accessions added to the CSV.")

# Example usage
file_path = 'C:/Users/bryan/Desktop/Analysis package/Table-2024-07-23.csv'
process_csv(file_path)
