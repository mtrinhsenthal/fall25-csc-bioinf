from dbg import DBG
from utils import read_data
import sys
import os

# Codon does not support sys.setrecursionlimit
# (stack size is fixed by the system, so we drop this)
# do i need to write dfs iteratively rather than recursively?

def main(argv: list[str]) -> None:
    if len(argv) < 2:
        print("Usage: ./program <input_folder>")
        return

    short1, short2, long1 = read_data(os.path.join("./", argv[1]))

    k: int = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    outpath: str = os.path.join("./", argv[1], "contig.fasta")
    with open(outpath, "w") as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write(">contig_%d\n" % i)
            f.write(c + "\n")


if __name__ == "__main__":
    main(sys.argv)
