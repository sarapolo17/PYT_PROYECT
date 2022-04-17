########################################################
############ GET SEQ and B-FACTORS FROM PDB ############
########################################################

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one

def get_seq_B_from_PDB (PDB_file):

    """ 
    Reads a PDB file and extracts the B factors for
    the alpha carbons of each amino acid. Returns a dictionary 
    with two keys:

        * "FASTA_seq": the sequence from the crystallized protein
                       in FASTA format

        * "B_val_list": list of tupples with the code and the B-value
                        for each alpha carbon of the amino acids in the
                        PDB file
    """

    # parse file

    parser = PDBParser()
    pdb_struct = parser.get_structure("struct_ID" ,PDB_file)

    # access atoms
    # get information

    seq_Bval_dict = {"FASTA_seq": "", "B_val_list": []}

    for model in pdb_struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "CA": # only for alpha carbons

                        # FASTA sequence

                        seq_Bval_dict["FASTA_seq"] += three_to_one(residue.get_resname())

                        # (residue, B-value) tupple

                        seq_Bval_dict["B_val_list"].append((three_to_one(residue.get_resname()), atom.get_bfactor()))
    
    return(seq_Bval_dict)
