import click
import logging
import os
from collection import get_primer_collection
from primer_score import get_primer_ranked_primers
from optimize import run_optimization
from primer_zone import get_primer_zones

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('primacy')

@click.group(context_settings=dict(show_default=True))
@click.option('--debug/--no-debug', default=False)
def cli(debug):
    if debug:
        logger.setLevel('DEBUG')
        logger.debug("DEBUG LOGGING")
    else:
        logger.setLevel("INFO")



@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'MULTIFASTA',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'OUTFILE',
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, writable=True))
@click.option(
    '--max-size', '-mx', default=30, type=click.IntRange(min=9, max=40),
    help="Max primer size."
)
@click.option(
    '--min-size', '-mn', default=18, type=click.IntRange(min=9, max=40),
    help="Min primer size."
)
@click.option(
    '--tm-max', '-tmx', default=65,
    type=click.IntRange(min=0),
    help="Set upper bound on melting temperature. Tm is calculated as the average for ambiguous primers."
)
@click.option(
    '--tm-min', '-tmn', default=50,
    type=click.IntRange(min=0),
    help="Set lower bound on melting temperature. Tm is calculated as the average for ambiguous primers."
)
@click.option(
    '--ambiguous-max', '-amb', default=5,
    type=click.IntRange(min=0),
    help="Set maximum number of ambiguous bases."
)
@click.option(
    '--gc-max', '-gcmx', default=60,
    type=click.IntRange(min=0, max=100),
    help="Set upper bound for GC%. GC% is the maximum value for ambiguous primers."
)
@click.option(
    '--gc-min', '-gcmn', default=30,
    type=click.IntRange(min=0, max=100),
    help="Set upper bound for GC%. GC% is the maximum value for ambiguous primers."
)
@click.option(
    '--gc-clamp-max', '-clmx', default=2,
    type=click.IntRange(min=0, max=5),
    help="Set upper bound on number of G and C bases in last 5 bases at 3' end of primer. More than 2 C/G bases at the end of the primer could cause issues with specificity. Set max to 5 and min to 0 to ignore GC clamp."
)
@click.option(
    '--gc-clamp-min', '-clmn', default=1,
    type=click.IntRange(min=0, max=5),
    help="Set lower bound on number of G and C bases in last 5 bases at 3' end of primer. At least 1 C/G base at the end of the primer could improve specificity."
)
@click.option(
    '--poly-max', '-plymx', default=5,
    type=click.IntRange(min=0),
    help="Set maximum homopolymer size. Homopolymers are the maximum for ambiguous primers."
)
@click.option(
    '--max_self-comp-score', '-scmx', default=8,
    type=click.IntRange(min=0),
    help="Set maximum self-complement score. Score is the maximum for ambiguous primers."
)
@click.option(
    "--hairpin-max-score", '-hpmx', default=5,
    type=click.INT,
    help="Maximum hybrid score for hairpin secondary structures. Score is the maximum for ambiguous primers."
)
@click.option(
    "--na-conc", '-na', default=50, type=click.FLOAT,
    help="Milimolar concentration of Na for calculating Tm (using Biopython nearest neighbor calculation)."
)
@click.option(
    "--k-conc", '-k', default=0, type=click.FLOAT,
    help="Milimolar concentration of K for calculating Tm."
)
@click.option(
    "--mg-conc", '-mg', default=0, type=click.FLOAT,
    help="Milimolar concentration of Mg for calculating Tm."
)
@click.option(
    "--tris", '-tc', default=0, type=click.FLOAT,
    help="Milimolar concentration of Tris for calculating Tm."
)
@click.option(
    "--dntps", '-dtc', default=0, type=click.FLOAT,
    help="Milimolar concentration of dNTPs for calculating Tm."
)
@click.option(
    "--saltcorr", '-sc', default=4,
    help="Salt correction value (see Biopython melting temp calculator)."
)
def primer_collection(
    multifasta, outfile, max_size, min_size, tm_max,
    tm_min, gc_max, gc_min, gc_clamp_max, gc_clamp_min, 
    ambiguous_max, poly_max, max_self_comp_score, hairpin_max_score,
    na_conc, k_conc, mg_conc, tris, dntps, saltcorr):
    """
    Provide a multifasta file of the target amplicon with flanking regions with headers
    in the format: >{TARGET_NAME}_{AMP_START}_{AMP_STOP}. All primers in the flanking 
    regions will be collected and filtered using provided parameters. The output is a table
    with all primers that pass the filters with their stats.     
    """
    if max_size < min_size:
        raise ValueError("Max primer size must be larger than min primer size.")
    if gc_max < gc_min:
        raise ValueError("Max GC\% must be higher than min GC%.")
    if tm_max < tm_min:
        raise ValueError("Max Tm must be higher than min Tm.")
    df = get_primer_collection(
        multifasta, min_size, max_size, tm_max,
    tm_min, gc_max, gc_min, gc_clamp_max, gc_clamp_min, 
    ambiguous_max, poly_max, max_self_comp_score, hairpin_max_score,
    na_conc, k_conc, mg_conc, tris, dntps, saltcorr)
    df.to_csv(outfile, sep='\t', index=False)
    

    

@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'PRIMERS', nargs=-1
)
@click.option(
    '--outfile', '-o',
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True)
)
@click.option(
    '--required-primers', '-r',
    type=click.Path(file_okay=True, dir_okay=False, readable=True),
    help="Table (tsv) containing primer name and sequences that should be included in solution optimization."
)
@click.option(
    '--max-amp-len', '-mx',
    type=click.IntRange(min=0), default=None,
    help="Maximum amplicon length."
)
@click.option(
    '--min-amp-len', '-mn',
    type=click.IntRange(min=0), default=None,
    help="Minimum amplicon length."
)
@click.option(
    '--pop-size', '-s', default=100,
    type=click.INT,
    help="Size of population for genetic optimization algorithm."
)
@click.option(
    '--iterations', '-it', default=100,
    type=click.IntRange(min=0),
    help="Number of generations to run the genetic optimization algorithm."
)
@click.option(
    '--mutation', '-mu', default=0.05,
    type=click.FloatRange(min=0, max=1),
    help="Set probability of a mutation occuring in offspring."
)
@click.option(
    '--crossover', '-xo', default=0.05,
    type=click.FloatRange(min=0, max=1),
    help="Set propbability of a cross over occuring in offspring."
)
def optimize(
    primers, outfile, required_primers,
    max_amp_len, min_amp_len, pop_size,
    iterations, mutation, crossover):

    res = run_optimization(
            primers, required_primers,
            max_amp_len, min_amp_len, pop_size,
            iterations, mutation, crossover)
    res.to_csv(outfile, sep='\t', index=False)



@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'MULTIFASTA',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'OUTFILE',
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, writable=True))
@click.option(
    '--amp-start', '-b',
    type=click.IntRange(min=0), required=True,
    help="Start position (inclusive) of the region to amplify"
)
@click.option(
    '--amp-stop', '-e',
    type=click.IntRange(min=0), required=True,
    help="End position (inclusive) of the region to amplify"
)
@click.option(
    "--flank-size", '-f',
    type=click.IntRange(min=0), default=100,
    help="Maximum size of the upstream and downstream region.")
@click.option(
    "--amp-name", '-n',
    type=click.STRING, default="",
    help="Name added to header in fasta file. Header is formated as '>{AMP_NAME}_{AMP_START}_{AMP_STOP}'. "
)
@click.option(
    "--ignore-genome", '-x',
    multiple=True, default=[],
    help="Ignore genome when identifying ambiguous bases. Repeat option to ignore multiple genomes."
)
@click.option(
    "--ignore-percent", '-p',
    type=click.FloatRange(min=0, max=100), default=0,
    help="Ignore variants with frequencies less than IGNORE_PERCENT"
)
def primer_zones(
        multifasta, outfile, amp_start, amp_stop, flank_size, amp_name, ignore_genome, ignore_percent):
    """
    Extract ampicon sequence and upstream and downstream primer zones from fasta given 
    the start and stop positions of the amplicon and flank size. A multifasta alignment 
    file can be provided and ambiguous bases will be included where there is variation
    in the alignment. Troublesome genomes (those that introduce too many ambiguous bases)
    can be ignored using the ignore_genome argument. IMPORTANT NOTE: this is meant to handle
    only signle nucleotide polymorphisms, insertions and deletions in the primer zones
    will need to dealt with separately. For this reason, the alignments should be to a reference
    sequence that has no gaps in the sequence.
    """
    get_primer_zones(
        multifasta, outfile, amp_start, amp_stop,
        flank_size, amp_name, ignore_genome, ignore_percent)




if __name__ == '__main__':
    cli()