#####################################
############ CLUSTAL MSA ############
#####################################

import sys
import os
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

def run_clustalo (input_FASTA, out_MSA_file, MSA_format = "fasta"):

    """ Wrapper to call CLUSTAL OMEGA using Python. Takes an input 
    file with FASTA sequences and creates an output file with the 
    MSA in the specified format (default: FASTA) """

    # write command

    clustalo_cline = ClustalOmegaCommandline(infile = input_FASTA, 
                                             outfile = out_MSA_file, 
                                             outfmt = MSA_format)

    # execute command

    clustalo_cline_exec = str(clustalo_cline)
    os.system(clustalo_cline_exec)  


def read_msa (MSA_file, MSA_format = "fasta"):

        """ Reads an MSA, separates it into records and returns
        it as an ID -> sequence dictionary """

        mult_aln = AlignIO.read(MSA_file, MSA_format)

        # inform of the results of the alignment

        sys.stderr.write("MSA created correctly with CLUSTAL OMEGA\n")
        sys.stderr.write(str(mult_aln) + "\n")

        # separate alignment records, create and return dict

        msa_records = {} 

        for record in mult_aln:
            msa_records[record.id] = record.seq

        return(msa_records)

def assign_msa_record_indices(msa_records):

    """ 
    Funtion that takes MSA records and returns a dictionary where each key
    has a list of tupples associated:

    * sequence identifier -> list of tupples

        * [x]: position in the MSA.

            * [0]: residue.

            * [1]: position in the MSA.

            * [2]: position in the PDB file.
    """

    msa_indices_dict = {} # dictionary for output

    for record_id, sequence in msa_records.items(): 

        msa_indices_dict[record_id] = [] # initialize dict for that record id

        msa_count = -1 # start counters so that first pos is 0
        pdb_count = -1 

        for aminoacid in sequence:
            msa_count += 1 # always increase msa count

            if aminoacid != '-': # not a gap
                 pdb_count += 1
                 msa_indices_dict[record_id].append((aminoacid, msa_count, pdb_count))

    return(msa_indices_dict)

