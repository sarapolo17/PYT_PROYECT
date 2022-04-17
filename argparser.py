###########################################
############ COMMAND LINE ARGS ############
###########################################

import argparse
parser = argparse.ArgumentParser (description = """This program calculates a flexibility score for the amino acids of protein sequence or a family of proteins.
                                                It returns a parseable text output file and a graphical representation of the scores. """)
parser.add_argument('-s', '--sequence',
                    dest = "input_sequence",
                    action = "store",
                    default = None,
                    help = "Input FASTA file with a single protein sequence.")

parser.add_argument('-f', '--family',
                    dest = "input_family",
                    action = "store",
                    default = None,
                    help = "Input file with MSA of a protein family in CLUSTAL format.")

parser.add_argument('-o', '--output',
                    dest = "outfile",
                    action = "store",
                    default = "flexibility_results",
                    help = "Flexibility analysis output files.")

options = parser.parse_args()

input_sequence = options.input_sequence
input_family = options.input_family
output_prefix = options.outfile
