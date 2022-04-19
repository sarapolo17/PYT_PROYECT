#########################################
############## MAIN PROGRAM #############
#########################################
### Import modules

import sys
import os
import argparse
from FlexiScor import consensus_seq
from FlexiScor import psi_BLAST_PSSM
from FlexiScor import top_hits_pdb
from FlexiScor import normalized_b_values
from FlexiScor import msa_clustal
from FlexiScor import calculate_flexibility
from FlexiScor import output_text_graph
from FlexiScor import database

if __name__ == "__main__":
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

    ### Create databases if they do not exist
    db_file_pdb=["./PDB_FASTA/PDB_FASTA.db.phr", "./PDB_FASTA/PDB_FASTA.db.pin", "./PDB_FASTA/PDB_FASTA.db.psq", "./PDB_FASTA/PDB_FASTA.db"]
    for file_db in db_file_pdb:
        if os.path.exists(file_db):
            continue
        else:
            database.database_pdb()

    db_file_uni=["./UniProt/Uniprot.db.fasta.phr", "./UniProt/Uniprot.db.fasta.pin", "./UniProt/Uniprot.db.fasta.psq", "./UniProt/Uniprot.db.fasta"]
    for file_uni in db_file_uni:
        if os.path.exists(file_uni):
            continue
        else:
            database.database_uniprot()
            
    if ((input_sequence and not input_family) or (not input_sequence and input_family)): # xor, only one input

        # correct input

        sys.stderr.write("Correct input. Processing...\n")

        # remove directory if it already exists

        if os.path.isdir("flexibility"):
            os.system("rm -r flexibility")

        # create new directory for results

        os.mkdir("flexibility")

        # move input to new directory

        if input_sequence:
            os.system("cp " + input_sequence + " ./flexibility")

        else:
             os.system("cp " + input_family + " ./flexibility")

        os.chdir("./flexibility")

    elif (input_sequence and input_family):
        raise SystemExit("Incorrect input. Please, make sure you introduce only one type of input (protein sequence or MSA of protein family in CLUSTAL format).")

    else:
        raise SystemExit("No input. Please, introduce an input option.")

    ### The input is a MSA of a protein family

    if input_family: # convert protein family to consensus sequence
        input_file = open('consensus_seq.fa', 'w')
        input_file.write(">query\n")

        try:
            input_file.write(str(consensus_seq.get_consensus_seq(input_family)))
            input_file.close()
            input_sequence = input_file.name

        except:
            input_file.close()
            raise SystemExit("Cannot obtain consensus sequence for protein family. Are you sure you introduced a MSA in CLUSTAL format?")

    ### PSI-BLAST
    ## Get PSSM

    psi_BLAST_PSSM.run_psiblast_homologues_PSSM(FASTA_seq = input_sequence,
                                                input_iters = 5,
                                                database = "../UniProt/UniProt.db.fasta",
                                                in_pssm_filename = None,
                                                out_pssm_filename = "psiblast_uniprot_5.pssm",
                                                out_hits_filename = "psiblast_uniprot_5.out",
                                                output_format = 6) # tabular

    ## Search in PDB with PSSM

    psi_BLAST_PSSM.run_psiblast_homologues_PSSM(FASTA_seq = input_sequence,
                                                input_iters = 1,
                                                database = "../PDB_FASTA/PDB_FASTA.db",
                                                in_pssm_filename = "./psiblast_uniprot_5.pssm",
                                                out_pssm_filename = "psiblast_pdb_1.pssm",
                                                out_hits_filename = "psiblast_pdb_1.out",
                                                output_format = 6) # tabular

    ### Get hits from PDB

    PDB_scores_to_download = top_hits_pdb.get_IDs_from_blastp_PDB("psiblast_pdb_1.out",
                                                                   number_hits = 10)

    list_id_PDB = top_hits_pdb.download_pdb (PDB_scores_to_download.keys(), os.getcwd() + "/PDB_downloads/")

    top_hits_pdb.select_chain_from_pdb(list_id_PDB, os.getcwd() + "/PDB_downloads/", os.getcwd() + "/PDB_downloads/split/")

    ## Extract protein sequences from the PDB files
    # File for CLUSTAL MSA

    hits_for_MSA_file = open('top_hits_split_chain.fa', 'w')

    all_seq_Bnorm = {}

    for PDB_split_file in (os.listdir(os.getcwd() + "/PDB_downloads/split")):
        seq_B = normalized_b_values.get_seq_B_from_PDB(os.getcwd() + "/PDB_downloads/split/" + PDB_split_file)

        # Normalize and append to big dictionary

        seq_B_norm = normalized_b_values.normalized_b_values(seq_B)
        all_seq_Bnorm[PDB_split_file.split(".")[0]] = seq_B_norm

    # Eliminate PDB homologs that could not be normalized form dictionary

    for PDB_homolog, seq_B_norm in all_seq_Bnorm.copy().items():
        for aa_Bnorm_tup in seq_B_norm["B_val_list"]:

            if aa_Bnorm_tup[1] == "?":
                sys.stderr.write("Eliminating one PDB homolog because of error during normalization\n")
                del all_seq_Bnorm[PDB_homolog]
                break
    # Create input for MSA: templates

    for PDB_homolog in all_seq_Bnorm:
        hits_for_MSA_file.write (">" + str(PDB_homolog) + "\n")
        hits_for_MSA_file.write (all_seq_Bnorm[PDB_homolog]["FASTA_seq"] + "\n")

    sys.stderr.write("Using " + str(len(all_seq_Bnorm)) + " PDB homologs with normalized B-values\n")

    # Add query to MSA input

    hits_for_MSA_file.write(">query\n")

    if input_family:

        hits_for_MSA_file.write(str(consensus_seq.get_consensus_seq(input_family)))

    elif input_sequence:
        fd = open (input_sequence, "r")
        for line in fd:
            if not line.startswith(">"):
                hits_for_MSA_file.write(line)

        fd.close()

    hits_for_MSA_file.close()

    ### CLUSTAL MSA

    msa_clustal.run_clustalo(hits_for_MSA_file.name, "PDB_MSA_clustalo.out.fasta")

    ### Put indices to MSA records
    ## Read MSA

    msa_records = msa_clustal.read_msa("PDB_MSA_clustalo.out.fasta")

    ## Assign indices

    msa_records_with_indices = msa_clustal.assign_msa_record_indices(msa_records)

    ## Calculate flexibility per position

    sys.stderr.write("Calculating flexibility...\n")

    incomplete_flexibility = calculate_flexibility.calculate_flex_pos(all_seq_Bnorm, msa_records_with_indices, PDB_scores_to_download)
    complete_flexibility = calculate_flexibility.complete_missing_scores(incomplete_flexibility)

    ## Create parseable text output

    try:
        output_text_graph.create_output_text(output_prefix, complete_flexibility)

    except:
        pass

    else:
        sys.stderr.write("The text output file has been created successfully!\n")

    ## Draw flexibility graph

    try:
        output_text_graph.draw_flex_line_col(output_prefix, complete_flexibility)

    except:
        pass

    else:
        sys.stderr.write("The plot has been created successfully!\n")
