#########################################
############ PSI-BLAST + IDs ############
#########################################

import os
import re
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast import Record
from Bio.Blast import NCBIXML
from Bio import SearchIO


def run_psiblast_homologues_PSSM (FASTA_seq, 
                                  input_iters,
                                  database,
                                  in_pssm_filename,
                                  out_pssm_filename,
                                  out_hits_filename,
                                  output_format = 6
                                  ):

    """ Wrapper for running PSI-BLAST from the command line. It takes a sequence 
    in FASTA format and returns an output file with the results of the search in the
    database and the PSSM corresponding to the MSA with the hits """

    # write command

    if in_pssm_filename: #  input PSSM provided
        psiblast_cline = NcbipsiblastCommandline(num_iterations = input_iters,
                                                 db = database,
                                                 in_pssm = in_pssm_filename,
                                                 out_pssm = out_pssm_filename,
                                                 out = out_hits_filename,
                                                 outfmt = output_format)
                                            
    else: # no input PSSM provided
        psiblast_cline = NcbipsiblastCommandline(query = FASTA_seq, 
                                                 num_iterations = input_iters,
                                                 db = database,
                                                 out_pssm = out_pssm_filename,
                                                 out = out_hits_filename,
                                                 outfmt = output_format)

    # execute command
    
    psiblast_cline_exec = str(psiblast_cline)
    os.system(psiblast_cline_exec)        

