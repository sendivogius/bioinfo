import os

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
GENOMES_DIR = os.path.join(DATA_DIR, 'genomes')

print('imported', GENOMES_DIR)