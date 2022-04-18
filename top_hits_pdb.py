####################################
############ PDB CHAINS ############
####################################
'''Functions for reading and processing the BLAST output against the PDB:
getting identifiers and scores, downloading PDB files and splitting chains'''

import os
import re
import sys
from urllib import request
from Bio import PDB


def get_IDs_from_blastp_PDB (blastp_outFile, number_hits):

    """ Parses  BLASTP output file and retrieves the PDB IDs and the scores
    for the hits. Returns them as a dictionary: ID -> score """

    # regular expresion for PDB IDs

    PDB_ID_score = {}
    
    n = 0 # counter for number of hits

    with open(blastp_outFile) as bpfile:
        for line in bpfile:
            records = line.split()
            PDB_ID = records[1]
            score = records[11]
            
            n+=1
            
            if (PDB_ID and n <= number_hits): # check that we have ID
                PDB_ID_score[PDB_ID.upper()] = score
                
        if PDB_ID_score: # check dictionary not empty
            return(PDB_ID_score)  
            
        else:
            raise SystemExit('No homologues found in UniProt. Exiting the program.')    

        

def download_pdb (list_id_chain, pdb_location):

    """ Separates the PDB IDs and chains and saves them as a tupple to a list (output).
    Then downloads the PDB files of the first hits of the psi-blast and saves them 
    to the folder 'pdb' """

    # Create tupple list

    list_id_chain_separated = [] # initialize list for (id, chain) tupples

    if os.path.isdir(pdb_location): # checks if directory exists
        pass
    else:
        os.mkdir(pdb_location) # create directory if it does not exist

    for hit_top in list_id_chain:
        PDB_all = hit_top.split('_') # separate id and chains

        id = PDB_all[0].upper() # pdb id
        chain = PDB_all[1].upper() # pdb homolog chain
        tuple_id_chain = (id, chain) # (id, chain) tupple

        list_id_chain_separated.append(tuple_id_chain) # append to list with tupples
    
    # Download PDB files

        PDB_file_location = pdb_location + id + '.pdb' # directory for downloaded files

        try: # search in PDB
            request.urlretrieve('https://files.rcsb.org/download/' + id + '.pdb', PDB_file_location) # download file
            sys.stderr.write("%s PDB file downloaded successfully\n" %id)

        except:
            raise SystemExit("ID PDB %s was not found." %id)
    
    return (list_id_chain_separated)


def select_chain_from_pdb(id_chain_pdb, pdb_location, pdb_split_location):

    """ Split the PDB files, creating new ones that contain only the data 
    from the selected (homologous) chains """

    if os.path.isdir(pdb_split_location): # checks if directory exists
        pass
    else:
        os.mkdir(pdb_split_location) # creates it if it does not exist

    for file_name in id_chain_pdb:

        current_file = pdb_location + file_name[0] + '.pdb' # name of downloaded file
        selected_chain = file_name[1] # chain
        output_file = pdb_split_location + file_name[0] + "_" + selected_chain + '.pdb' # name output file
        
        found_chain = False

        with open (output_file, 'w') as pdb_split: # create for writing
            with open (current_file) as pdb_not_split: # open for reading
                for line in pdb_not_split:
                    if line.startswith('ATOM') and line[21] == selected_chain:
                        # select residues of the chain
                        pdb_split.write(line) # write to the output file
                        found_chain = True # change Boolean
                        
                    elif line.startswith("TER") and found_chain: # avoid repetition in output
                        break


        sys.stderr.write("%s PDB file split successfully, chain %s selected \n" %(file_name[0], file_name[1]))
