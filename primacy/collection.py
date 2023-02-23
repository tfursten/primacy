import pandas as pd
import logging
import itertools as it
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



def get_primer_collection(
    seq, seq_name, orientation, 
    min_size, max_size):
    return get_primer_stats(
        process_primers(
            seq, seq_name, orientation, min_size, max_size))



def get_all_primers_of_a_size(seq, seq_name, size, orientation):
    primers = []
    seq = seq.upper()
    for i in range(len(seq) - size + 1):
        primer = Seq(seq[i: i + size])
        primer = primer.reverse_complement() if orientation=="R" else primer
        primers.append(
            ["{0}_{1}_{2}".format(seq_name, i, size),
            orientation,
            "".join(primer)])
    return pd.DataFrame(primers, columns=['PrimerName', 'Orientation', 'Seq'])



def process_primers(seq, seq_name, orientation, min_size, max_size):
    primers = pd.DataFrame([])
    # add primers for all size ranges
    for sz in range(min_size, max_size + 1):
        primers = pd.concat(
            [primers,
            get_all_primers_of_a_size(
                seq, seq_name, sz, orientation)])
    primers['Target'] = seq_name
    return primers[['Target', 'PrimerName', 'Orientation', 'Seq']]
        

def calculate_tm(seq, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, dnac1=500, dnac2=500, saltcorr=3):
    return mt.Tm_NN(seq, Na=Na, K=K, dnac1=dnac1, dnac2=dnac2, saltcorr=saltcorr,
                        Tris=Tris, Mg=Mg, dNTPs=dNTPs)

def get_number_of_degenerate_bases(seq):
    bases = ['A', 'C', 'G', 'T']
    return len([base for base in list(seq) if base not in bases])

def calculate_gc(seq):
    return ((seq.count("C") + seq.count('G')) / float(len(seq))) * 100


def calculate_homopolymer(seq):
    homopolymers = [
        list(hp) for run, hp
        in it.groupby(seq)]
    n_homopolymer = 0
    for hp in homopolymers:
        if len(hp) >= 3:
            n_homopolymer += len(hp)
    longest_homopolymer = len(max(homopolymers, key=len))
    return longest_homopolymer

def gc_clamp(seq):
    last_5_bases = seq[-5:]
    gc_clamp = last_5_bases.count("C") + last_5_bases.count("G")
    if 1 < gc_clamp < 3:
        return True
    else:
        return False

def get_primer_stats(df):
    df['Ambiguous'] = df['Seq'].apply(get_number_of_degenerate_bases)
    df['Tm'] = df['Seq'].apply(calculate_tm)
    df['GC%'] = df['Seq'].apply(calculate_gc)
    df['Homopolymer'] = df['Seq'].apply(calculate_homopolymer)
    df['GC_Clamp'] = df['Seq'].apply(gc_clamp)
    return df

