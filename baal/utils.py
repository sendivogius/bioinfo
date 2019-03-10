def is_dna(dna):
    return set(dna) <= set('ACTG')


def get_invalid_chars(dna):
    return set(dna) - set('ACTG')


def get_genome(species):
    from .. import DATA_DIR
    print(DATA_DIR)
    raw_lines = open(rf'{DATA_DIR}\genomes\{species}.txt').readlines()
    dna = ''.join(line.strip() for line in raw_lines)
    wrong_chars = get_invalid_chars(dna)
    if wrong_chars:
        print(f'DNA for {species} contain invalid characters: {wrong_chars}')
    return dna
