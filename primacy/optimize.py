import logging
import copy
import pandas as pd
import numpy as np
from numpy.random import default_rng
import itertools as it
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from collection import calculate_hybrid_score

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('primacy')

class Primer(object):
    def __init__(self, seq, name=None, amplicon=None, flank=None):
        self.seq = seq
        self.name = name
        self.amplicon = amplicon
        self.flank = flank

class PrimerPair(object):
    def __init__(self, amplicon_table, max_sz, min_sz, seed=None):
        self.max_sz = max_sz
        self.min_sz = min_sz
        self.amplicon_table = amplicon_table
        self.amplicon = amplicon_table['Amplicon'].unique()
        self.seed = seed
        self.forward_primer, self.reverse_primer = self.get_primer_pair(
            f=None, r=None)
        


    def get_primer_pair(self, f=None, r=None):
        """
        f and r are primer objects
        """
        loop_iter = 0
        if f==None and r==None:
            while True:
                loop_iter += 1
                forward = self.amplicon_table[self.amplicon_table['Flank'] == 'F'].sample(1, random_state=self.seed).iloc[0]
                reverse = self.amplicon_table[(self.amplicon_table['Flank'] == 'R')]
                if self.max_sz != None:
                    reverse = reverse[reverse['Position'] <= (forward['Position'] + self.max_sz)]
                if self.min_sz != None:
                    reverse = reverse[reverse['Position'] >= (forward['Position'] + self.min_sz)]
                if len(reverse):
                    reverse = reverse.sample(1, random_state=self.seed).iloc[0]
                    break
                if loop_iter > 1000:
                    raise ValueError(
                        "Amplicon {0}: Having trouble finding amplicons of the correct size, please adjust parameters.".format(
                        self.amplicon))
        elif f != None:
            forward = self.amplicon_table.loc[f.name]
            reverse = self.amplicon_table[(self.amplicon_table['Flank'] == 'R')]
            if self.max_sz != None:
                reverse = reverse[reverse['Position'] <= (forward['Position'] + self.max_sz)]
            if self.min_sz != None:
                reverse = reverse[reverse['Position'] >= (forward['Position'] + self.min_sz)]
            if len(reverse):
                reverse = reverse.sample(1, random_state=self.seed).iloc[0]
            else:
                # Run again replacing both primers
                return self.get_primer_pair(f=None, r=None)
        elif r != None:
            reverse = self.amplicon_table.loc[r.name]
            forward = self.amplicon_table[(self.amplicon_table['Flank'] == 'F')]
            if self.max_sz != None:
                forward = forward[forward['Position'] <= (reverse['Position'] + self.max_sz)]
            if self.min_sz != None:
                forward = forward[forward['Position'] >= (reverse['Position'] + self.min_sz)]
            if len(forward):
                forward = forward.sample(1, random_state=self.seed).iloc[0]
            else:
                return self.get_primer_pair(f=None, r=None)
        return [Primer(seq=forward['Seq'], name=forward.name), Primer(seq=reverse['Seq'], name=reverse.name)]


def get_result_table(result, primer_options):
    primers = []
    for ppair in result:
        primers.append(ppair.forward_primer.name)
        primers.append(ppair.reverse_primer.name)
    primers = primer_options.loc[primers]
    primers = primers.reset_index()
    primers[['CrossHybridScoreMean', 'CrossHybridScoreMax']] = primers['Seq'].apply(
        lambda x: hybrid_score_stats(x, [Seq(p) for p in primers['Seq'].values]))
    amp_sizes = {
        a: df[df['Flank'] == 'R']['Position'].values[0] - df[df['Flank'] == 'F']['Position'].values[0] + 1
        for a, df in primers.groupby('Amplicon')}    
    primers['AmpSize'] = primers['Amplicon'].apply(lambda x: amp_sizes.get(x))
    return primers    


def run_optimization(
    primer_options, required_primers, max_amp_len, min_amp_len, sample_sz, seed):
    primers = []
    for primer_set in primer_options:
        primers.append(pd.read_csv(primer_set, sep="\t", index_col='PrimerName'))
    primers = pd.concat(primers)
    if required_primers:
        req_primers = list([
            Primer(name=n, seq=r.iloc[0]) for n, r in 
            pd.read_csv(required_primers, sep="\t", header=None, index_col=0, comment="#").iterrows()])
    else:
        req_primers = []
    res = optimization(primers, req_primers, sample_sz, max_amp_len, min_amp_len, seed)
    return get_result_table(res, primers)

    

def hybrid_score_stats(primer, primers):
    scores = []
    for p in primers:
        scores.append(calculate_hybrid_score(primer, p))
    return pd.Series([np.mean(scores), np.max(scores)])


def hybrid_score(new_pair, existing_primers):
    scores = []
    for p in new_pair:
        for e in existing_primers:
            scores.append(calculate_hybrid_score(p, e))
    # minimize both the max and average score
    return np.mean(scores) + np.max(scores)                      
    
    

def get_all_primer_seqs(primer_pairs, required_primers=[]):
    primers = []
    for pair in primer_pairs:
        primers.append(Seq(pair.forward_primer.seq))
        primers.append(Seq(pair.reverse_primer.seq))
    for primer in required_primers:
        primers.append(Seq(primer.seq))
    return primers


def optimization(primers, req_primers, reps, max_sz, min_sz, seed):
    current_panel = []
    rng = np.random.default_rng(seed)
    targets = primers['Amplicon'].unique()
    rng.shuffle(targets)
    for target in targets:
        logger.info("Adding primers for amplicon: {}".format(target))
        amp_table = primers[primers['Amplicon'] == target]
        primer_pairs = [PrimerPair(amp_table, max_sz, min_sz, seed) for rep in range(reps)]
        scores = [hybrid_score(get_all_primer_seqs([p]), get_all_primer_seqs(current_panel, req_primers)) for p in primer_pairs]
        logger.debug("Scores: {}".format(", ".join(list(map(str, scores)))))
        logger.info("Best Score: {}".format(np.min(scores)))
        current_panel.append(primer_pairs[np.argmin(scores)])
    return current_panel
 











