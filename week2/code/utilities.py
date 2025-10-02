from matrix import FrequencyPositionMatrix

######### JASPAR FORMAT #########
def jaspar_write(motifs, format):
    """Return the representation of motifs in "pfm" or "jaspar" format."""
    letters = "ACGT"
    lines: List[str] = []
    if format == "pfm":
        motif = motifs[0]
        counts: FrequencyPositionMatrix = motif.counts
        for letter in letters:
            terms = []
            for i in range(motif.length):
                terms.append(f"{counts[letter][i]:6.2f}")
            line = f"{' '.join(terms)}\n"
            lines.append(line)
    elif format == "jaspar":
        for m in motifs:
            counts: FrequencyPositionMatrix = m.counts
            # try:
            #     matrix_id = m.matrix_id
            # except AttributeError:
            #     matrix_id = None
            line = f">None {m.name}\n"
            lines.append(line)
            for letter in letters:
                terms = []
                for i in range(m.length):
                    terms.append(f"{counts[letter][i]:6.2f}")
                line = f"{letter} [{' '.join(terms)}]\n"
                lines.append(line)
    else:
        raise ValueError(f"Unknown JASPAR format {format}")

    # Finished; glue the lines together
    text = "".join(lines)
    return text

######### TRANSAC FORMAT #########

MULTIPLE_VALUE_KEYS = {"BF", "OV", "HP", "BS", "HC", "DT", "DR", "CC"}
REFERENCE_KEYS = {"RX", "RA", "RT", "RL"}
def transac_write(motifs):
    blocks = []
    sections = (
        ("AC", "AS"),  # Accession
        ("ID",),  # ID
        ("DT", "CO"),  # Date, copyright
        ("NA",),  # Name
        ("DE",),  # Short factor description
        ("TY",),  # Type
        ("OS", "OC"),  # Organism
        ("HP", "HC"),  # Superfamilies, subfamilies
        ("BF",),  # Binding factors
        ("P0",),  # Frequency matrix
        ("BA",),  # Statistical basis
        ("BS",),  # Factor binding sites
        ("CC",),  # Comments
        ("DR",),  # External databases
        ("OV", "PV"),  # Versions
    )
    for m in motifs:
        lines = []
        for section in sections:
            blank = False
            for key in section:
                if key == "P0":
                    length = m.length
                    if length == 0: continue
                    sequence = m.degenerate_consensus
                    letters = sorted(m.alphabet)
                    line = "      ".join(["P0"] + letters)
                    lines.append(line)
                    for i in range(length):
                        line = " ".join([f"{i+1:0>2}"] + [f"{m.counts[_][i]:6.20g}" for _ in letters]) + f"      {sequence[i]}"
                        lines.append(line)
                    blank = True
            if blank:
                line = "XX"
                lines.append(line)
        line = "//"
        lines.append(line)
        block = "\n".join(lines) + "\n"
        blocks.append(block)
    text = "".join(blocks)
    return text

######### CLUSTERBUSTER FORMAT #########
def clusterbuster_write(motifs, precision = 0):
    """Return the representation of motifs in Cluster Buster position frequency matrix format.

    By default (`precision=0`) Cluster Buster position frequency matrices will be written
    with integer values.
    If a higher precision value is set, Cluster Buster position frequency matrices will be
    written as floats with `x` decimal places.
    """
    lines = []
    for m in motifs:
        lines.append(f">{m.name}\n")
        for ACGT_counts in zip(
            m.counts["A"], m.counts["C"], m.counts["G"], m.counts["T"]
        ):
            line = "\t".join([f"{round(val, precision)}" for val in ACGT_counts]) + "\n"
            lines.append(line)

    # Finished; glue the lines together.
    text = "".join(lines)

    return text