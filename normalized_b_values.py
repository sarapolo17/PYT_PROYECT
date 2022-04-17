########################################################
############ GET SEQ and B-FACTORS FROM PDB ############
########################################################

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one
import statistics

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
    list_residues=[]
    list_b_factor_raw=[]

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

                        # # (residue, B-value) tupple
                        # list_residues.append(three_to_one(residue.get_resname()))
                        # list_b_factor_raw.append(atom.get_bfactor())
                        # # list with all the b-values before normalizing

                        seq_Bval_dict["B_val_list"].append((three_to_one(residue.get_resname()), atom.get_bfactor()))

    return(seq_Bval_dict)

def normalized_b_values(dict_seq_bval):
    """ 
    Takes a dictionary with 2 keys as input:
    * "FASTA_seq": the sequence from the crystallized protein
                   in FASTA format

    * "B_val_list": list of tupples with the code and the B-value
                    for each alpha carbon of the amino acids in the
                    PDB file

    It returns a dictionary with 2 keys:
    * "FASTA_seq": the sequence from the crystallized protein
                   in FASTA format

    * "B_val_list": list of tupples with the code and the B-value _normalized_
                    for each alpha carbon of the amino acids in the
                    PDB file
    """
    
    list_seq_bval = dict_seq_bval["B_val_list"]
    
    # value of dict with key "B_val_list"
    
    list_b_val_before = [] # list b-values before normalizing
    list_b_val_norm = [] # list b-values after normalizing
    list_res = []
    res_b_val_norm = []
    
    for residue_b_val in list_seq_bval:
        list_res.append(residue_b_val[0])
        list_b_val_before.append(residue_b_val[1])
        
    for b_val_not_norm in list_b_val_before:
        b_val_norm = round((b_val_not_norm-statistics.mean(list_b_val_before))/statistics.stdev(list_b_val_before), 5)
        # b val norm with 5 decimal positions
        list_b_val_norm.append(b_val_norm)
        res_b_val_norm = list(zip(list_res, list_b_val_norm))
        
        # join two list in a third one
    dict_seq_bval["B_val_list"] = ""
    # delete previous value of the dictionary
    dict_seq_bval["B_val_list"] = res_b_val_norm
    # asign value to dict key
    return (dict_seq_bval)


