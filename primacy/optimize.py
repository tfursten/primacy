import logging
import pandas as pd
import numpy as np
import itertools as it
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq



def run_optimization(
            scored_primers, required_primers, seed,
            iterations, kmer):
    np.random.seed(seed)
    all_primers = pd.DataFrame([])
    for primer_set in scored_primers:
        all_primers = pd.concat(
            [all_primers,
            pd.read_csv(primer_set, sep="\t")])
    optimization_targets = all_primers['Target'].unique()
    if required_primers:
        all_primers = pd.concat(
            [all_primers,
            pd.read_csv(required_primers, sep="\t")])
    # Randomly select primers from each target
    init_sample = all_primers.groupby(
        "Target").sample(
            n=1, random_state=seed)[['PrimerName', 'Target', 'Seq']]
    selected_targets = init_sample['Target'].values
    selected_primers = init_sample['PrimerName'].values
    selected_seqs = init_sample['Seq'].values
    distance_matrix = initialize_distance_matrix(selected_seqs, kmer)
    # Mask for required targets so they are not chosen for optimzation
    optimization_target_idx = np.where(
        np.isin(selected_targets, optimization_targets))[0]

    iteration = 0
    while iteration < iterations:
        # pick a random column
        iteration += 1
        choice = np.random.choice(optimization_target_idx)
        logger = logging.getLogger(__name__)

        logger.info("Iteration: {}".format(iteration))
        logger.info("Average Similarity: {}".format(np.mean(distance_matrix)))
        logger.info("Median Similarity: {}".format(np.median(distance_matrix)))
        logger.info("Max Similarity: {}".format(np.max(distance_matrix)))
        logger.info("Selecting target: {}".format(selected_targets[choice]))

        # Get all sequences options for this target
        seqs, names = get_seqs_for_target(selected_targets[choice], all_primers)
        result = optimize_target(
            seqs, names, selected_seqs,
            distance_matrix, choice, kmer)
        if result:
            distance_matrix = result['distance_matrix']
            selected_seqs = result['seqs']
            selected_primers[choice] = result['primer_name']
        else:
            logger.info("No improvement")
    
    return all_primers[all_primers['PrimerName'].isin(selected_primers)]



def optimize_target(
    target_seqs, primer_names, current_seqs,
    distance_matrix, idx, kmer):
    """
    Run optimization on selected target until 
    a better primer is found or all have been tested.
    """
    current_mean = distance_matrix[idx].mean()
    for seq, name in zip(target_seqs, primer_names):
        current_seqs[idx] = seq
        dd = update_distance_matrix(
            distance_matrix, current_seqs, idx, kmer
        )
        new_mean = dd[idx].mean()
        if new_mean < current_mean:
            return {
                'distance_matrix': dd,
                'seqs': current_seqs,
                'primer_name': name,
                }
    return None




def get_seqs_for_target(target, all_primers):
    df = all_primers[all_primers['Target'] == target]
    return df['Seq'].values, df['PrimerName'].values


def update_distance_matrix(distance_matrix, seqs, position, kmer):
    """
    Update only distance related to the choses position
    """
    dd = distance_matrix.copy()
    for i in range(len(seqs)):
        d = get_similarity_between_primers(seqs[i], seqs[position], kmer)
        dd[i][position] = d
        dd[position][i] = d
    return dd




def find_worst_pair(distance_matrix, mask):
    # Mask non-optimization targets 
    dd = ma.masked_array(distance_matrix, mask=mask)
    return np.unravel_index(
        np.argmax(dd, axis=None),
        dd.shape), dd.max()



def initialize_distance_matrix(seqs, kmer):
    """
    Get a pairwise distance matrix for seqs
    """
    m = np.zeros((len(seqs), len(seqs)))
    for i in range(len(seqs)):
        for j in range(i, len(seqs)):
            d = get_similarity_between_primers(seqs[i], seqs[j], kmer)
            m[i][j] = d
            m[j][i] = d
    return m

def expand_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = ambiguous_dna_values
    return tuple(map("".join, it.product(*map(d.get, seq))))

def get_similarity_between_primers(p1, p2, kmer):
    """
    Generate a similarity score based on the number of kmers
    that align between primers.
    """
    p2 = str(Seq(p2).reverse_complement())
    longer_primer = p1 if len(p1) >= len(p2) else p2
    shorter_primer = p1 if longer_primer == p2 else p2
    # Get all unique kmers in longer primer, including expanded ambiguous bases.
    longer_kmers = set(
        it.chain.from_iterable(
            [expand_ambiguous_dna(longer_primer[i: i + kmer])
        for i in range(len(longer_primer) - kmer + 1)]))
    hits = 0
    for i in range(len(shorter_primer) - kmer + 1):
        expanded_kmer = expand_ambiguous_dna(
            shorter_primer[i: i + kmer])
        for k in expanded_kmer:
            if k in longer_kmers:
                hits += 1
                break
    return hits / (len(longer_primer) - kmer + 1)


