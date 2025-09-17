from codon_dbg import DBG
from codon_utils import read_data
import sys

# Codon does not support sys.setrecursionlimit
# stack size is fixed by the system

if __name__ == "__main__":
    argv = sys.argv

    # relative_path: str = "../data/"
    folder: str = argv[1]
    input_path: str = folder
    # input_path: str = relative_path + folder

    short1, short2, long1 = read_data(input_path)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])
    
    outpath: str = input_path + "/contig.fasta"
    with open(outpath, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            # f.write('>contig_%d\n' % i)
            f.write(f'>contig_{i}\n')
            f.write(c + '\n')