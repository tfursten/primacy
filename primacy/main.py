import click
import logging
import os
from collection import get_primer_collection
from primer_score import get_primer_ranked_primers
from optimize import run_optimization
from primer_zone import get_primer_zones

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('vast')

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
    'SEQ_NAME',
    type=click.STRING
)
@click.argument(
    'SEQ',
    type=click.STRING
)
@click.argument(
    'OUTPUT',
    type=click.Path(
        file_okay=True, dir_okay=False, writable=True
    )
)
@click.option(
    '--orientation', '-o', required=True,
    type=click.Choice(['F', 'R'], case_sensitive=False),
    help="Select orientation of primers. If R is selected, primers will be the reverse complement of the provided sequence"
)
@click.option(
    '--max_size', '-mx', default=25, type=click.IntRange(min=9, max=40), show_default=True,
    help="Max primer size."
)
@click.option(
    '--min_size', '-mn', default=12, type=click.IntRange(min=9, max=40), show_default=True,
    help="Min primer size."
)
def primer_collection(seq, seq_name, orientation, output, max_size, min_size):
    """
    Provided a sequence for a primer region,
    returns a table of all possible primers (within user specified sizes) with stats
    for Tm, homopolymer runs, GC%, 3' GC clamp, and number of ambiguous bases.
    
    """
    if max_size < min_size:
        raise ValueError("Max primer size must be larger than min primer size.")
    df = get_primer_collection(
        seq, seq_name, orientation, min_size, max_size)
    df.to_csv(output, sep='\t', index=False)
    
@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'primer_collection',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=True))
@click.argument(
    'outfile',
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True)
)
@click.option(
    '--homopolymer_cutoff', '-h', default=5, show_default=True,
    type=click.IntRange(min=0),
    help="Set max size of a homopolymer run."
)
@click.option(
    '--tm_max', '-tmx', default=60, show_default=True,
    type=click.IntRange(min=0),
    help="Set upper bound on melting temperature."
)
@click.option(
    '--tm_min', '-tmn', default=50, show_default=True,
    type=click.IntRange(min=0),
    help="Set lower bound on melting temperature."
)
@click.option(
    '--ambiguous_cutoff', '-amb', default=5, show_default=True,
    type=click.IntRange(min=0),
    help="Set maximum number of ambiguous bases."
)
@click.option(
    '--gc_max', '-gcmx', default=60, show_default=True,
    type=click.IntRange(min=0, max=100),
    help="Set upper bound for GC%."
)
@click.option(
    '--gc_min', '-gcmn', default=40,
    type=click.IntRange(min=0, max=100), show_default=True,
    help="Set upper bound for GC%."
)
@click.option(
    '--homopolymer_weight', '-hw', default=1, show_default=True,
    type=click.IntRange(min=1),
    help="Set a weight for the importance of homopolymers in scoring."
)
@click.option(
    '--ambiguous_weight', '-aw', default=1, show_default=True,
    type=click.IntRange(min=1),
    help="Set a weight for the importance of ambiguous bases in scoring."
)
@click.option(
    '--gc_weight', '-gcw', default=1, show_default=True,
    type=click.IntRange(min=1),
    help="Set a weight for the importance of GC content in scoring."
)
@click.option(
    '--gc_clamp_weight', '-gccw', default=1, show_default=True,
    type=click.IntRange(min=1),
    help="Set a weight for the importance of a GC clamp in scoring."
)
@click.option(
    '--percentile', '-p', default=90, show_default=True,
    type=click.FloatRange(min=0, max=100),
    help="Return the top N percentile of primers."
)
def primer_ranking(
    primer_collection, outfile,
    homopolymer_cutoff, tm_max, tm_min,
    ambiguous_cutoff, gc_max, gc_min,
    homopolymer_weight,
    ambiguous_weight, gc_weight,
    gc_clamp_weight, percentile):
    """
    Rank and weigh primers, return only the top scoring primers. Primers that fall outside of
    the provided cutoffs are given a score of 0 for that feature but may still be considered
    if they remain in the top percentile overall.
    """
    primers = get_primer_ranked_primers(primer_collection,
        homopolymer_cutoff, tm_max, tm_min,
        ambiguous_cutoff, gc_max, gc_min,
        homopolymer_weight,
        ambiguous_weight, gc_weight,
        gc_clamp_weight, percentile)
    primers.to_csv(outfile, sep="\t", index=False)
    

@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'SCORED_PRIMERS', nargs=-1
)
@click.option(
    '--outfile', '-o',
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True)
)
@click.option(
    '--required_primers', '-r',
    type=click.Path(file_okay=True, dir_okay=False, readable=True),
    help="List of primers that must be included in the solution. Tab separated header should include: Target PrimerName Orientation Seq"
)
@click.option(
    '--seed', '-s',
    type=click.INT,
    help="Provide a seed for random number generator."
)
@click.option(
    '--iterations', '-it', default=100, show_default=True,
    type=click.FloatRange(min=0),
    help="Set number of iterations to run"
)
@click.option(
    '--kmer', '-k', default=3, show_default=True,
    type=click.IntRange(min=3),
    help="Set kmer size for finding cross hybridization score."
)
def optimize(
    scored_primers, outfile, required_primers, seed,
    iterations, kmer):

    res = run_optimization(
            scored_primers, required_primers, seed,
            iterations, kmer)
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