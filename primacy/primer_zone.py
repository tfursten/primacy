import logging
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.SeqRecord import SeqRecord

LOGGER = logging.getLogger('primacy')


def remove_ignored_genomes(align, ignore_genomes):
    align_minus_ignored = MultipleSeqAlignment([])
    for genome in align:
        if genome.id in ignore_genomes:
            LOGGER.info("Removing genome from alignment {}".format(genome.id))
        else:
            align_minus_ignored.add_sequence(genome.id, str(genome.seq).upper())
    return align_minus_ignored
    
        
def get_zone_positions(multifasta, amp_start, amp_stop, flank_size):
    """
    Get zone start and stop point by adding flank size to the amp_start and amp_stop.
    If zone is larger than the sequence, adjust the zone to the max length.
    """
    with open(multifasta, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seq_len = len(record.seq)
            break
    zone_start = max(0, amp_start - flank_size)
    zone_stop = min(seq_len, amp_stop + flank_size)
    return zone_start, zone_stop, amp_start - zone_start, amp_stop - zone_start


def get_consensus_sequence(align):
    seq = []
    ambiguous_map = {v: k for k, v in ambiguous_dna_values.items() if v != 'X'}
    for position in range(align.get_alignment_length()):
        seq.append(
            ambiguous_map.get("".join(sorted(set(align[:, position].replace("X", "")))), 'N'))
    return Seq("".join(seq))



def get_primer_zones(
    multifasta, outfile, amp_start, amp_stop, flank_size, amp_name, ignore_genomes):
    zone_start, zone_stop, amp_start, amp_stop = get_zone_positions(
        multifasta, amp_start, amp_stop, flank_size)
    align = remove_ignored_genomes(
        AlignIO.read(multifasta, "fasta")[:, zone_start : zone_stop + 1], ignore_genomes)
    LOGGER.info(
        "Total target length is {0}".format(
        align.get_alignment_length()))
    LOGGER.info(
        "Total number of genomes in alignment: {0}".format(len(align))
    )
    seq_name = "_".join([str(val) for val in [amp_name, amp_start, amp_stop] if val])
    with open(outfile, 'w') as outfile:
        record = SeqRecord(
            get_consensus_sequence(align),
            id=seq_name,
            name=seq_name,
            description=seq_name,
        )
        SeqIO.write(record, outfile, 'fasta')
    



    
   

