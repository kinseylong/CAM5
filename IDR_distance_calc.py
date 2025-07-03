import os
import csv
import pandas as pd
import Levenshtein
from concurrent.futures import ProcessPoolExecutor, as_completed

target_sequence = "MAAKRSSNSAEYKEKNGRRKSHCRIL"
                            
path_to_tsv = "data/TAIR10/IDRs.csv"
output_filepath_csv = "data/TAIR10/ArabidopsisReducedAlphabet.csv"

##reduced alphabets
four = {"L": "A", "V": "A", "I": "A", "M": "A", "C": "A",
        "A": "B", "G": "B", "S": "B", "T": "B", "P": "B",
        "F": "C", "Y": "C", "W": "C",
        "E": "D", "D": "D", "N": "D", "Q": "D", "K": "D", "R": "D", "H": "D"}

eight = {"L": "A", "V": "A", "I": "A", "M": "A", "C": "A",
         "A": "B", "G": "B",
         "S": "C", "T": "C",
         "P": "D",
         "F": "E", "Y": "E", "W": "E",
         "E": "F", "D": "F", "N": "F", "Q": "F",
         "K": "G", "R": "G",
         "H": "H"}

ten = {"L": "A", "V": "A", "I": "A", "M": "A",
       "C": "B",
       "A": "C",
       "G": "D",
       "S": "E", "T": "E",
       "P": "F",
       "F": "G", "Y": "G", "W": "G",
       "E": "H", "D": "H", "N": "H", "Q": "H",
       "K": "I", "R": "I",
       "H": "J"}

twelve = {"L": "A", "V": "A", "I": "A", "M": "A",
          "C": "B",
          "A": "C",
          "G": "D",
          "S": "E", "T": "E",
          "P": "F",
          "F": "G", "Y": "G",
          "W": "H",
          "E": "I", "Q": "I",
          "D": "J", "N": "J",
          "K": "K", "R": "K",
          "H": "L"}

fifteen = {"L": "A", "V": "A", "I": "A", "M": "A",
           "C": "B",
           "A": "C",
           "G": "D",
           "S": "E",
           "T": "F",
           "P": "G",
           "F": "H", "Y": "H",
           "W": "I",
           "E": "J",
           "Q": "K",
           "D": "L",
           "N": "M",
           "K": "N", "R": "N",
           "H": "O"}

eighteen = {"L": "A", "M": "A",
            "V": "B", "I": "B",
            "C": "C",
            "A": "D",
            "G": "E",
            "S": "F",
            "T": "G",
            "P": "H",
            "F": "I",
            "Y": "J",
            "W": "K",
            "E": "L",
            "D": "M",
            "N": "N",
            "Q": "O",
            "K": "P", 
            "R": "Q",
            "H": "R"}

#convert sequence into reduced alphabet sequences
def reduce_alphabet(sequence):
    def translate(sequence, reduction_dict):
        return ''.join([reduction_dict.get(res, res) for res in sequence])
    
    reduced_seq_four = translate(sequence, four)
    reduced_seq_eight = translate(sequence, eight)
    reduced_seq_ten = translate(sequence, ten)
    reduced_seq_twelve = translate(sequence, twelve)
    reduced_seq_fifteen = translate(sequence, fifteen)
    reduced_seq_eighteen = translate(sequence, eighteen)

    return reduced_seq_four, reduced_seq_eight, reduced_seq_ten, reduced_seq_twelve, reduced_seq_fifteen, reduced_seq_eighteen


#calculate levenshetin distance using sliding window, returns min distance and corresponding window
def sliding_window_levenshtein(str1, str2, window_size=None, step=2):
    if window_size is None:
        window_size = min(len(str1), len(str2))

    if len(str1) <= len(str2):
        target = str2
        query = str1
    else:
        target = str1
        query = str2

    distances = []
    for i in range(0, len(target) - window_size + 1, step):
        window = target[i:i + window_size]
        dist = Levenshtein.distance(query, window)
        distances.append((dist, window))

    return sorted(distances)[0]

def process_entry(row):
    entry = row['Entry']
    sequence = row['Sequence']
    
    #reduce alphabet
    four_seq, eight_seq, ten_seq, twelve_seq, fifteen_seq, eighteen_seq = reduce_alphabet(sequence)
    #calculate levenshtein distance against each alphabet
    four_dist = sliding_window_levenshtein(four_seq, target_seq_four)
    eight_dist = sliding_window_levenshtein(eight_seq, target_seq_eight)
    twelve_dist = sliding_window_levenshtein(twelve_seq, target_seq_twelve)
    eighteen_dist = sliding_window_levenshtein(eighteen_seq, target_seq_eighteen)

    results = [entry, sequence, four_dist, eight_dist, twelve_dist, eighteen_dist]
    return results

# Load and filter uniprot
IDRs = pd.read_csv(path_to_tsv, sep="\t")
IDRs = IDRs[["Entry", "Sequence"]].drop_duplicates()

#calculate reduced alphabet sequences of target sequence
target_seq_four, target_seq_eight, target_seq_ten, target_seq_twelve, target_seq_fifteen, target_seq_eighteen = reduce_alphabet(target_sequence)


# Remove already written entries
if os.path.exists(output_filepath_csv) and os.path.getsize(output_filepath_csv) > 0:
    written_entries = set()
    with open(output_filepath_csv, "r") as idr_file:
        reader = csv.DictReader(idr_file, delimiter="\t")
        for row in reader:
            written_entries.add(row["Entry"])
    IDRs = IDRs[~IDRs["Entry"].isin(written_entries)]
    print(f"Filtered uniprot DataFrame to exclude {len(written_entries)} already written entries.")
else:
    with open(output_filepath_csv, "w") as w:
        w.write("Entry\tSequence\tFour Dist\tEight Dist\tTwelve Dist\tEighteen Dist\n")

# Parallel processing
rows = IDRs.to_dict(orient="records")

num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
with ProcessPoolExecutor(max_workers=num_workers) as executor, open(output_filepath_csv, "a") as f:
    futures = [executor.submit(process_entry, row) for row in rows]
    for future in as_completed(futures):
        results = future.result()
        for entry, sequence, four_dist, eight_dist, twelve_dist, eighteen_dist in results:
            f.write(f"{entry}\t{sequence}\t{four_dist}\t{eight_dist}\t{twelve_dist}\t{eighteen_dist}\n")