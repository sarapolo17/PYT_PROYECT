################################
############ SET UP ############
################################

""" Set up file with with information and needs to run the program"""
from setuptools import setup, find_packages, dist

requirements = ['Bio==1.3.8', 'matplotlib==3.4.3','numpy==1.20.3', 'setuptools==58.0.4', 'testresources']

setup(
    name='FlexiScor',
    version='1.0',
    packages=['FlexiScor'],
    url='https://github.com/sarapolo17/PYT_PROYECT',
    author='Vicente Ángel Ledesma Martín',
    author_email='vicenteangel.ledesma01@estudiant.upf.edu and sara.polo02@estudiant.upf.edu',
    description='FlexiScor, calculating flexibilities',
    install_requires=requirements,
    python_requires='>=3.9',
    py_modules=['main', 'database', 'calculate_flexibility', 'consensus_seq', 'get_Bfact_seq_from_PDB',
    'msa_clustal', 'normalized_b_values', 'output_text_graph', 'psi_BLAST_PSSM', 'top_hits_pdb']
)
