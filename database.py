from urllib import request
import os
import sys
import gzip
import shutil

def database_pdb():
    sys.stderr.write("Creating PDB database...\n")
    os.system("mkdir PDB_FASTA")
    os.chdir("./PDB_FASTA")
    request.urlretrieve("https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz", 'PDB_FASTA.db.gz')
    with gzip.open('PDB_FASTA.db.gz', 'rb') as f_in:
        with open('PDB_FASTA.db', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.system("makeblastdb -in PDB_FASTA.db -parse_seqids -blastdb_version 5 -dbtype prot")
    os.chdir("./..")




def database_uniprot():
    sys.stderr.write("Creating UniProt database...\n")
    os.system("mkdir UniProt")
    os.chdir("./UniProt")
    request.urlretrieve("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", 'UniProt.db.fasta.gz')
    with gzip.open('UniProt.db.fasta.gz', 'rb') as f_in:
        with open('UniProt.db.fasta', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.system("makeblastdb -in UniProt.db.fasta -parse_seqids -blastdb_version 5 -dbtype prot")
    os.chdir("./..")
