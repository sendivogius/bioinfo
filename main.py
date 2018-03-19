import utils


def find_ori(file, txt=''):
    if file:
        dna = ''.join(f.strip() for f in open(file).readlines())
        # utils.plot_skew(dna)
        min_skew = utils.get_min_skew_posiion(dna)[0]
        print(min_skew)
        L = 500
        dna_box = dna[min_skew-L:min_skew+L]
    else:
        dna_box = txt
    k = 9
    d = 1
    most_freq_kmer = utils.frequent_kmers(dna_box, k, d , True)
    print('Found ORI', most_freq_kmer)


if __name__ == "__main__":
    find_ori('data\\Salmonella_enterica.txt')
    find_ori('data\\E_coli.txt')

