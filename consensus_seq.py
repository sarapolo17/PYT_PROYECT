import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

def get_consensus_seq(aln_file):

    """Creates a consensus sequence
    from an input MSA file in CLUSTAL format.
    If the percentage of the most common residue type is
    greater than 70%, then we will add that
    residue type; otherwise an X will be added."""
    
    alignment = AlignIO.read(aln_file, 'clustal')
    summary_align = AlignInfo.SummaryInfo(alignment)
    con_seq = summary_align.dumb_consensus()
    
    return (con_seq)
    


