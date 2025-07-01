import os
import csv
import requests
import pandas as pd

path_to_uniprot_tsv = "data/TAIR10/ArabidopsisUniprot.tsv"
#should have columns: Entry, AlphaFoldDB, Sequence
output_idr_file_path_csv = "data/TAIR10/IDRs.csv"

#determine low-confidence AlphaFold positions
def DL_range(enumerated_values, min_length):
    ranges = []
    current_range = []

    # Iterate through enumerated_values
    for index, value in enumerated_values:
        if value in ['D', 'L']:
            if not current_range:
                current_range = [index, index]  # Start a new range
            else:
                current_range[1] = index  # Extend the current range
        else:
            if current_range:
                # Check if the range length meets the minimum length
                if current_range[1] - current_range[0] + 1 >= min_length:
                    ranges.append(tuple(current_range))  # Save the completed range
                current_range = []  # Reset the range

    # Append the last range if it exists and meets the minimum length
    if current_range and current_range[1] - current_range[0] + 1 >= min_length:
        ranges.append(tuple(current_range))
    return ranges

### API to get IDRs from an AlphaFold ID and a protein sequence
def alphafold_idr(alphafold_ID, protein_sequence, min_IDR_length = 15):
    url = f'https://alphafold.ebi.ac.uk/files/AF-{alphafold_ID}-F1-confidence_v4.json'
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        
        if 'confidenceCategory' in data and 'confidenceScore' in data:
            confidence_values = data['confidenceCategory']
            #get positions with continuous D and L regions (IDRs)
            enumerated_values = list(enumerate(confidence_values))
            IDR_positions = DL_range(enumerated_values, min_IDR_length)
            #get confidence scores for the positions
            confidence_scores = data['confidenceScore']
            IDR_scores = []
            IDR_seqs = []
            IDR_starts = []
            IDR_ends = []
            for start, end in IDR_positions:
                IDR_scores.append(confidence_scores[start:end + 1])
                IDR_seqs.append(protein_sequence[start:end + 1])
                IDR_starts.append(start)
                IDR_ends.append(end)
            return IDR_starts, IDR_ends, IDR_seqs, IDR_scores
    #error
    return None, None, None, None


###run 
uniprot = pd.read_csv(path_to_uniprot_tsv, sep="\t")
uniprot = uniprot[uniprot['AlphaFoldDB'].notnull()] #remove empty AlphaFoldDB entries
uniprot = uniprot[["Entry", "AlphaFoldDB", "Sequence"]].drop_duplicates()

# Check if the file exists and is not empty
if os.path.exists(output_idr_file_path_csv) and os.path.getsize(output_idr_file_path_csv) > 0:
    # Read the existing entries from IDRs.csv
    written_entries = set()
    with open(output_idr_file_path_csv, "r") as idr_file:
        reader = csv.DictReader(idr_file, delimiter="\t")
        for row in reader:
            written_entries.add(row["Entry"])  # Collect already written entries
    # Filter the uniprot DataFrame to exclude already written entries
    uniprot = uniprot[~uniprot["Entry"].isin(written_entries)]
    print(f"Filtered uniprot DataFrame to exclude {len(written_entries)} already written entries.")
else:
    with open(output_idr_file_path_csv, "w") as w:
        w.write("Entry\tStart\tEnd\tSequence\tConfidenceScore\n")  # Write header for the output file

for index, row in uniprot.iterrows():
    alphafold_ID = row['AlphaFoldDB'].replace(";", "")
    sequence = row['Sequence']
    IDR_starts, IDR_ends, IDR_seq, IDR_scores = alphafold_idr(alphafold_ID, sequence)
    
    if IDR_starts:  # If there are IDRs found
        with open(output_idr_file_path_csv, "a") as f:
            for start, end, seq, score in zip(IDR_starts, IDR_ends, IDR_seq, IDR_scores):
                f.write(f"{row['Entry']}\t{start}\t{end}\t{seq}\t{score}\n")
