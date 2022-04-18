#########################################
############## SET UP FILE ##############
#########################################

""" Set up file with information needed to run the program """
from setuptools import setup, find_packages, dist

with open('requires.txt') as f:
    requirements = f.readlines()

setup(
    name='FlexiScore',
    version='1.0',
    packages=find_packages(),
    url='https://github.com/sarapolo17/PYT_PROYECT',
    author='Vicente Ángel Ledesma Martín',
    author_email='vicenteangel.ledesma01@estudiant.upf.edu and sara.polo02@estudiant.upf.edu',
    description='FlexiScore, calculating flexibilities',
    install_requires=requirements,
    python_requires='>=3.9',
    scripts=['main.py', 'argparser.py', 'database.py', 'calculate_flexibility.py', 'consensus_seq.py', 'get_Bfact_seq_from_PDB.py',
    'msa_clustal.py', 'normalized_b_values.py', 'output_text_graph.py', 'psi_BLAST_PSSM.py', 'top_hits_pdb.py']
)
