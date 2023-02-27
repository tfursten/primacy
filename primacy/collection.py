import pandas as pd
import logging
import itertools as it
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.IUPACData import ambiguous_dna_values

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Star used for overhangs
ambiguous_dna_values['*'] = "*"


class Amplicon(object):
    def __init__(self, seq):
        """
        Given a sequence object, create an amplicon
        """
        self.name, self.start, self.stop = seq.id.split("_")
        self.start = int(self.start)
        self.stop = int(self.stop)
        self.left_flank = seq.seq[:self.start]
        self.right_flank = seq.seq[self.stop + 1:]
        self.amplicon = seq.seq[self.start: self.stop + 1]
        self.left_flank_sz = len(self.left_flank)
        self.right_flank_sz = len(self.right_flank)
        self.amplicon_sz = len(self.amplicon)
    
        
    def get_left_primers(self, min_sz, max_sz):
        primers = []
        for size in range(min_sz, max_sz + 1):
            for pos in range(self.left_flank_sz - size + 1):
                primers.append(
                    SeqRecord(
                    self.left_flank[pos: pos + size + 1],
                    id="{0}_{1}_{2}_F".format(self.name, pos, size)))
        logger.info(
            "{0}-{1}: {2} primers found with lengths between the range [{3}:{4}]".format(
            self.name, "F", len(primers), min_sz, max_sz 
            ))
        return primers

    def get_right_primers(self, min_sz, max_sz):
        # reverse complement
        primers = []
        for size in range(min_sz, max_sz + 1):
            for pos in range(self.right_flank_sz - size + 1):
                primers.append(
                    SeqRecord(
                    self.right_flank[pos: pos + size + 1].reverse_complement(),
                    id="{0}_{1}_{2}_R".format(self.name, pos, size)))
        logger.info(
            "{0}-{1}: {2} primers found with lengths between the range [{3}:{4}]".format(
            self.name, "R", len(primers), min_sz, max_sz 
            ))
        return primers



def get_primer_collection(
    multifasta, min_size, max_size, tm_max,
    tm_min, gc_max, gc_min,  gc_clamp_max, gc_clamp_min,
    ambiguous_max, poly_max,
    max_self_comp_score, hairpin_max_score,
    na_conc, k_conc, mg_conc, tris, dntps, saltcorr):
    results = []
    with open(multifasta, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            amp = Amplicon(record)
            results.append(
                process_primers(
                amp.get_left_primers(min_size, max_size),
                amp.name, "F", tm_max, tm_min, gc_max,
                gc_min,  gc_clamp_max, gc_clamp_min,
                ambiguous_max, poly_max,
                max_self_comp_score, hairpin_max_score,
                na_conc, k_conc, mg_conc, tris, dntps, saltcorr))
            results.append(
                process_primers(
                amp.get_right_primers(min_size, max_size),
                amp.name, "R", tm_max, tm_min, gc_max, 
                gc_min, gc_clamp_max, gc_clamp_min,
                ambiguous_max, poly_max,
                max_self_comp_score, hairpin_max_score,
                na_conc, k_conc, mg_conc, tris, dntps, saltcorr))
    results = pd.concat(results)
    results = results[
        [
        'PrimerName', 'Amplicon', 'Flank',
        'Position', 'Seq', 'Length', 'AmbiguousBases',
        'Tm', 'GC', 'LongestHomopolymer', 'SelfHybridScore',
        'HairpinScore', 'GCClamp']]
    return results


def get_primer_position(primer_name):
    """
    Return the far end position of primers based on name.
    Primer names are formated as "{AMPID}_{START}_{LENGTH}_{FLANK}"
    """
    _, idx, length, flank = primer_name.split("_")
    if flank == "F":
        return int(idx)
    else:
        return int(idx) + int(length)




def process_primers(primer_list, amp_name, flank, tm_max, tm_min, gc_max,
                gc_min, gc_clamp_max, gc_clamp_min,
                ambiguous_max, poly_max,
                max_self_comp_score, hairpin_max_score,
                na_conc, k_conc, mg_conc, tris, dntps, saltcorr):
    """
    Take list of primer sequences. Return a dataframe of filtered
    primers with stats.
    """
    primers = pd.DataFrame(
        [[p.id, p.seq, flank, amp_name] for p in primer_list],
        columns=["PrimerName", "Seq", "Flank", "Amplicon"])
    primers['Length'] = primers['Seq'].apply(lambda x: len(x))
    # Ambiguous bases (do first to reduce number of Tm calculations)
    start_sz = primers.shape[0]
    primers['Position'] = primers['PrimerName'].apply(get_primer_position)
    primers['AmbiguousBases'] = primers['Seq'].apply(get_number_of_degenerate_bases)
    primers = primers[primers['AmbiguousBases'] <= ambiguous_max]
    logger.info(
        "{0}-{1}: {2} primer(s) had more than {3} ambiguous bases.".format(
        amp_name, flank, start_sz - primers.shape[0], ambiguous_max))
    # GC
    start_sz = primers.shape[0]
    primers['GC'] = primers['Seq'].apply(calculate_gc)
    primers = primers[(primers['GC'] >= gc_min) & (primers['GC'] <= gc_max)]
    logger.info(
        "{0}-{1}: {2} primers had GC content outside the range [{3}:{4}]".format(
        amp_name, flank, start_sz - primers.shape[0], gc_min, gc_max
        ))
    # Homopolymers
    start_sz = primers.shape[0]
    primers['LongestHomopolymer'] = primers['Seq'].apply(calculate_homopolymer)
    primers = primers[primers['LongestHomopolymer'] <= poly_max]
    logger.info(
        "{0}-{1}: {2} primers had homopolymers longer than {3}".format(
        amp_name, flank, start_sz - primers.shape[0], poly_max))
    # Selfcomplementary
    start_sz = primers.shape[0]
    primers['SelfHybridScore'] = primers['Seq'].apply(calculate_self_hybrid_score)
    primers = primers[primers['SelfHybridScore'] <= max_self_comp_score]
    logger.info(
        "{0}-{1}: {2} primers had self hybridization scores greater than {3}".format(
        amp_name, flank, start_sz - primers.shape[0], max_self_comp_score))
    # Hairpin
    start_sz = primers.shape[0]
    primers['HairpinScore'] = primers['Seq'].apply(calculate_hairpin)
    primers = primers[primers['HairpinScore'] <= hairpin_max_score]
    logger.info(
        "{0}-{1}: {2} primers had hairpin scores greater than {3}".format(
        amp_name, flank, start_sz - primers.shape[0], hairpin_max_score))
    # GC clamp
    start_sz = primers.shape[0]
    primers['GCClamp'] = primers['Seq'].apply(gc_clamp)
    primers = primers[(primers['GCClamp'] >= gc_clamp_min) & (primers['GCClamp'] <= gc_clamp_max)]
    logger.info(
        "{0}-{1}: {2} primers had GCClamp with CG content outside the range [{3}:{4}]".format(
        amp_name, flank, start_sz - primers.shape[0], gc_clamp_min, gc_clamp_max))
    # TM
    start_sz = primers.shape[0]
    primers['Tm'] = primers['Seq'].apply(
        lambda x: calculate_tm(
        x, Na=na_conc, K=k_conc, Mg=mg_conc, Tris=tris, dNTPs=dntps, saltcorr=saltcorr))
    primers = primers[(primers['Tm'] >= tm_min) & (primers['Tm'] <= tm_max)]
    logger.info(
        "{0}-{1}: {2} primers had melting temperatures (Tm) outside the range [{3}:{4}]".format(
        amp_name, flank, start_sz - primers.shape[0], tm_min, tm_max))
    if primers.shape[0] == 0:
        logger.warn("Amplicon {0} {1} Flank: No primers passed filter! Use less stringent parameters.".format(
            amp_name, flank
        ))
    return primers


def expand_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    return tuple(map("".join, it.product(*map(ambiguous_dna_values.get, seq))))

def calculate_tm(seq, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, dnac1=500, dnac2=500, saltcorr=1):
    expanded_seqs = expand_ambiguous_dna(seq)
    avg_tm = []
    for seq in expanded_seqs:    
        avg_tm.append(mt.Tm_NN(seq, nn_table=mt.DNA_NN1, Na=Na, K=K, dnac1=dnac1, dnac2=dnac2, saltcorr=saltcorr,
                        Tris=Tris, Mg=Mg, dNTPs=dNTPs))
        if len(avg_tm) > 100:
            break
    return sum(avg_tm) / len(avg_tm)

def get_number_of_degenerate_bases(seq):
    bases = ['A', 'C', 'G', 'T']
    return len([base for base in list(seq) if base not in bases])

def calculate_gc(seq):
    gc_count = 0
    for base in seq:
        base = ambiguous_dna_values[base]
        if "C" in base or "G" in base:
            gc_count += 1
    return (gc_count / len(seq)) * 100


def calculate_homopolymer(seq):
    current_base = None
    current_run = 0
    max_run = 1
    # TODO: adjust to work with ambiguous bases
    for base in seq:
        if current_base == None:
            current_base = set(ambiguous_dna_values[base])
        else:
            intersection = set(current_base).intersection(set(ambiguous_dna_values[base]))
            if len(intersection):
                current_run += 1
                current_base = intersection
            else:
                current_base = set(ambiguous_dna_values[base])
                max_run = max(max_run, current_run)
                current_run = 1
    return max_run


def calculate_self_hybrid_score(primer):
    return calculate_hybrid_score(primer, primer)

def calculate_hybrid_score(primer_1, primer_2):
    p1 = primer_1 + ("*" * (len(primer_2) - 1))
    p2 = "*" * (len(primer_1) - 1) + primer_2.reverse_complement()
    max_score = 0
    while len(p1):
        overlap_filter = [0 if "*" in primers else 1 for primers in zip(p1, p2)]
        overlap_p1 = list(it.compress(p1, overlap_filter))
        overlap_p2 = list(it.compress(p2, overlap_filter))
        score = score_calculator(overlap_p1, overlap_p2)
        max_score = max(max_score, score)
        p1 = p1[:-1]
        p2 = p2[1:] 
    return max_score

def score_calculator(seq1, seq2):
    score = 0
    for b1, b2 in zip(seq1, seq2):
        if len(set(ambiguous_dna_values[b1]).intersection(set(ambiguous_dna_values[b2]))):
            score += 1
        else:
            score -= 1
    return score


def hairpin_calculator_helper(seq, min_stem_sz, loop_sz):
    max_stem_sz = len(seq) - min_stem_sz - loop_sz
    seq = seq + "*" * len(seq)
    max_score = 0
    for stem_sz in range(min_stem_sz, max_stem_sz + 1):
        p1 = seq[0: stem_sz]
        p2 = seq[stem_sz + loop_sz: stem_sz * 2 + loop_sz]
        score = score_calculator(p1, p2)
        max_score = max(score, max_score)
    return max_score


def calculate_hairpin(seq):
    loop4 = hairpin_calculator_helper(seq, min_stem_sz=2, loop_sz=4)
    loop5 = hairpin_calculator_helper(seq, min_stem_sz=2, loop_sz=5)
    return max(loop4, loop5)
    

def gc_clamp(seq):
    last_5_bases = seq[-5:]
    gc_count = 0
    for base in last_5_bases:
        base = ambiguous_dna_values[base]
        if "C" in base or "G" in base:
            gc_count += 1
    return gc_count

