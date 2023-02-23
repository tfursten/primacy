import pandas as pd
import numpy as np
import logging



def get_range_score(values, max_val, min_val, weight):
    return np.multiply(
        weight,
        np.logical_and(
            np.greater_equal(values, min_val),
            np.less_equal(values, max_val)))

def get_below_cutoff_score(values, cutoff, weight):
    return np.multiply(
        weight,
        np.less_equal(values, cutoff)
    )

def filter_primers_for_tm(primers, tm_min, tm_max):
    # Do a hard cutoff for melting temperature
    # But if no primers remain, increase the Tm
    # range until some primers are present.
    logger = logging.getLogger(__name__)
    filtered_primers = primers[(primers['Tm'] >= tm_min) & (primers['Tm'] <= tm_max)]
    while not filtered_primers.shape[0]:
        tm_min -= 1
        tm_max += 1
        filtered_primers = primers[(primers['Tm'] >= tm_min) & (primers['Tm'] <= tm_max)]
        logger.warning("No primers within Tm range, increasing range to {0}-{1}.".format(tm_min, tm_max))
    return filtered_primers



def get_primer_ranked_primers(primer_collection,
        homopolymer_cutoff, tm_max, tm_min,
        ambiguous_cutoff, gc_max, gc_min,
        homopolymer_weight,
        ambiguous_weight, gc_weight,
        gc_clamp_weight, percentile):
    
    primers = pd.read_csv(primer_collection, sep="\t")
    # Check to see if any primers remain in range
    primers = filter_primers_for_tm(primers, tm_min, tm_max)
    primers['GC_score'] = get_range_score(primers['GC%'], gc_max, gc_min, gc_weight)
    primers['Homopolymer_score'] = get_below_cutoff_score(
        primers['Homopolymer'], homopolymer_cutoff, homopolymer_weight)
    primers['Ambiguous_score'] = get_below_cutoff_score(
        primers['Ambiguous'], ambiguous_cutoff, ambiguous_weight)
    primers['Clamp_score'] = np.multiply(primers['GC_Clamp'], gc_clamp_weight)
    primers['Combined_score'] = primers[[
        'GC_score', 'Homopolymer_score',
        'Ambiguous_score', 'Clamp_score']].sum(axis=1)
    primers = primers[primers['Combined_score'] >= primers['Combined_score'].quantile(percentile/100)]    
    return primers
