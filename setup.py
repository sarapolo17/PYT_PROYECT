################################
############ SET UP ############
################################

""" Set up file with with information and needs to run the program""" 
from setuptools import setup, find_packages, dist

requirements = ['Bio==1.3.8', 'matplotlib==3.4.3','numpy==1.20.3', 'setuptools==58.0.4']

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
