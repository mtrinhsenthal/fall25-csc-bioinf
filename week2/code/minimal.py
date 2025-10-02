import __init__ as motifs
from typing import List, Dict, Optional

def read(handle: File):
    motif_number = 0
    record = Record()
    _read_version(record, handle)
    _read_alphabet(record, handle)
    _read_background(record, handle)

    while True:
        for line in handle:
            if line.startswith("MOTIF"):
                break
        else:
            return record
        name = line.split()[1]
        motif_number += 1
        length, num_occurrences, evalue = _read_motif_statistics(handle)
        counts = _read_lpm(record, handle, length, num_occurrences)

        motif = motifs.Motif(alphabet=record.alphabet, counts=counts)
        motif.background = record.background
        motif.length = motif.counts.length
        # motif.num_occurrences = num_occurrences
        # motif.evalue = evalue
        motif.name = name
        record.append(motif)
        assert len(record) == motif_number
    return record

class Record():
    version: str
    datafile: str
    command: str
    alphabet: Optional[str]
    background: Dict[str, float]
    sequences: List[str]
    data: List[motifs.Motif]

    def __init__(self):
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = None
        self.background = {}
        self.sequences = []
        self.data = List[motifs.Motif]()

    def append(self, item: motifs.Motif):
        self.data.append(item)

    def __len__(self):
        return len(self.data)
    
    def __iter__(self):
        return self.data.__iter__()
    
    def __getitem__(self, key):
        if isinstance(key, str):
            for motif in self.data:
                if motif.name == key:
                    return motif
        else:
            return self.data.__getitem__(key)

def _read_version(record: Record, handle):
    for line in handle:
        if line.startswith("MEME version"):
            break
    else:
        raise ValueError("Improper input file")
    line = line.strip().split()
    record.version = line[2]

def _read_alphabet(record: Record, handle: File):
    for line in handle:
        if line.startswith("ALPHABET"):
            break
    else:
        raise ValueError("Unexpected end of stream")
    if not line.startswith("ALPHABET= "):
        raise ValueError(f"Line does not start with 'ALPHABET':\n{line}")
    line = line.strip().replace("ALPHABET= ", "")
    if line == "ACGT":
        al = "ACGT"
    elif line == "ACGU":
        al = "ACGU"
    else:
        raise ValueError("Only parsing of DNA and RNA motifs is implemented")
    record.alphabet = al

def _read_background(record: Record, handle: File):
    for line in handle:
        if line.startswith("Background letter frequencies"):
            background_freqs = []
            for line in handle:
                line = line.rstrip()
                if line:
                    background_freqs.extend([float(freq) for i,freq in enumerate(line.split(" ")) if i % 2 == 1])
                else:
                    break
            if not background_freqs:
                raise ValueError(
                    "Unexpected end of stream: Expected to find line starting background frequencies."
                )
            break
    else:
        raise ValueError(
            "Improper input file. File should contain a line starting background frequencies."
        )
    record.background = dict(zip(record.alphabet, background_freqs))

def _read_motif_statistics(handle):
    for line in handle:
        if line.startswith("letter-probability matrix:"):
            break
    num_occurrences = int(line.split("nsites=")[1].split()[0]) if line.find("nsites=") != -1 else 20
    length = int(line.split("w=")[1].split()[0]) if line.find("w=") != -1 else None
    evalue = float(line.split("E=")[1].split()[0]) if line.find("E=") != -1 else 0.0
    return length, num_occurrences, evalue

def _read_lpm(record, handle, length, num_occurrences):
    counts = [[], [], [], []]
    for line in handle:
        freqs = line.split()
        if len(freqs) != 4:
            break
        counts[0].append(round(float(freqs[0]) * num_occurrences))
        counts[1].append(round(float(freqs[1]) * num_occurrences))
        counts[2].append(round(float(freqs[2]) * num_occurrences))
        counts[3].append(round(float(freqs[3]) * num_occurrences))
        if length and len(counts[0]) == length:
            break
    c = dict(zip(record.alphabet, counts))
    return c