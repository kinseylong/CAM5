import os
import csv
import requests
import pandas as pd
import math
import metapredict
from concurrent.futures import ProcessPoolExecutor, as_completed

path_to_uniprot_tsv = "data/TAIR10/ArabidopsisUniprot.tsv"
output_idr_file_path_csv = "data/TAIR10/IDRs_metapredict.csv"

def DL_range(enumerated_values, min_length):
    ranges = []
    current_range = []
    for index, value in enumerated_values:
        if value in ['D', 'L']:
            if not current_range:
                current_range = [index, index]
            else:
                current_range[1] = index
        else:
            if current_range:
                if current_range[1] - current_range[0] + 1 >= min_length:
                    ranges.append(tuple(current_range))
                current_range = []
    if current_range and current_range[1] - current_range[0] + 1 >= min_length:
        ranges.append(tuple(current_range))
    return ranges

def alphafold_idr(alphafold_ID, protein_sequence, min_IDR_length=15):
    url = f'https://alphafold.ebi.ac.uk/files/AF-{alphafold_ID}-F1-confidence_v4.json'
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if 'confidenceCategory' in data and 'confidenceScore' in data:
                confidence_values = data['confidenceCategory']
                enumerated_values = list(enumerate(confidence_values))
                IDR_positions = DL_range(enumerated_values, min_IDR_length)
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
    except Exception as e:
        print(f"Error fetching {alphafold_ID}: {e}")
    return None, None, None, None

def process_entry(row):
    alphafold_ID = row['AlphaFoldDB'].replace(";", "")
    sequence = row['Sequence']
    entry = row['Entry']
    IDR_starts, IDR_ends, IDR_seq, IDR_scores = alphafold_idr(alphafold_ID, sequence)
    results = []
    if IDR_starts:
        for start, end, seq, score in zip(IDR_starts, IDR_ends, IDR_seq, IDR_scores):
            try:
                disorder = metapredict.percent_disorder(seq)
            except Exception as e:
                disorder = math.nan
            results.append((entry, start, end, seq, score, disorder))
    return results

# Load and filter uniprot
uniprot = pd.read_csv(path_to_uniprot_tsv, sep="\t")
uniprot = uniprot[uniprot['AlphaFoldDB'].notnull()]
uniprot = uniprot[["Entry", "AlphaFoldDB", "Sequence"]].drop_duplicates()

# Remove already written entries
if os.path.exists(output_idr_file_path_csv) and os.path.getsize(output_idr_file_path_csv) > 0:
    written_entries = set()
    with open(output_idr_file_path_csv, "r") as idr_file:
        reader = csv.DictReader(idr_file, delimiter="\t")
        for row in reader:
            written_entries.add(row["Entry"])
    uniprot = uniprot[~uniprot["Entry"].isin(written_entries)]
    print(f"Filtered uniprot DataFrame to exclude {len(written_entries)} already written entries.")
else:
    with open(output_idr_file_path_csv, "w") as w:
        w.write("Entry\tStart\tEnd\tSequence\tConfidenceScore\tMetapredictDisorderPercent\n")

# Parallel processing
rows = uniprot.to_dict(orient="records")

num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
with ProcessPoolExecutor(max_workers=num_workers) as executor, open(output_idr_file_path_csv, "a") as f:
    futures = [executor.submit(process_entry, row) for row in rows]
    for future in as_completed(futures):
        results = future.result()
        for entry, start, end, seq, score, metapredict_score in results:
            f.write(f"{entry}\t{start}\t{end}\t{seq}\t{score}\t{metapredict_score}\n")