import os
import pandas as pd
from Bio import SeqIO
from primer import get_primer_stats, specificity
from primer_score import (
    get_scores, convert_to_distance, get_primer_dataframe)

def read_primers(fasta_path):
    primer_data = []
    with open(fasta_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            primer_data.append(
                [record.id, 0, len(record.seq), 0, record.seq ])
    return pd.DataFrame(primer_data, columns=["seq_id", "start",
        "length", "flank", "sequence"] )


def score_existing_primers(
    project_name,
    fasta_path,
    outpath, threads):
    tm_params = {'Na': 50, "K": 0, "Tris": 0, 'Mg': 0, "dNTPs": 0}
    gc_opt = {'max': 60, 'min': 30}
    tm_opt = 55
    primer_weights = {"tm": 1, "gc": 1, "homopolymer": 1, "dimer": 1, "specificity": 1, "degenerate": 1}
    set_weights = [1,1,1]
    background_paths = []


    primers = read_primers(fasta_path)
    print(primers)
    primers = specificity(primers, background_paths, outpath, threads)
    primers = get_primer_stats(primers, tm_params)
    print(primers)
    primer_scores = get_scores(
        convert_to_distance(
            get_primer_dataframe(primers),
            tm_opt, gc_opt), primer_weights)
    primer_scores.to_csv(os.path.join(outpath,
        "{}_primer_scores.csv".format(project_name)))


if __name__ == '__main__':
    score_existing_primers("option_3",
    "../../Data/Primer_Interactions/for_adam/chimp_killer/Option_3.fasta",
    "../../Data/Primer_Interactions/for_adam/chimp_killer/",
    2)