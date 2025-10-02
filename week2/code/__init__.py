from python import Bio.Align as Align
import numpy as np
import matrix
from typing import Optional, Dict, List
from python import Bio.Seq as Seq
from python import urllib.parse as parse
from python import urllib.request as request
import utilities

def create(instances, alphabet="ACGT"):
    alignment = Align.Alignment(instances)
    return Motif(alphabet=alphabet, alignment=alignment)

def parse(handle, fmt, strict=True):
    """Parse an output file from a motif finding program.

    Currently supported formats (case is ignored):
     - AlignAce:         AlignAce output file format
     - ClusterBuster:    Cluster Buster position frequency matrix format
     - XMS:              XMS matrix format
     - MEME:             MEME output file motif
     - MINIMAL:          MINIMAL MEME output file motif
     - MAST:             MAST output file motif
     - TRANSFAC:         TRANSFAC database file format
     - pfm-four-columns: Generic position-frequency matrix format with four columns. (CIS-BP, HOMER, HOCOMOCO, Neph, Tiffin)
     - pfm-four-rows:    Generic position-frequency matrix format with four row. (ScerTF, YeTFaSCo, hDPI, iDMMPMM, FlyFactorSurvey, Cys2His2 Zinc Finger Proteins PWM Predictor)
     - pfm:              JASPAR-style position-frequency matrix
     - jaspar:           JASPAR-style multiple PFM format
     - sites:            JASPAR-style sites file

    As files in the pfm and sites formats contain only a single motif,
    it is easier to use Bio.motifs.read() instead of Bio.motifs.parse()
    for those.

    For example:

    >>> from Bio import motifs
    >>> with open("motifs/alignace.out") as handle:
    ...     for m in motifs.parse(handle, "AlignAce"):
    ...         print(m.consensus)
    ...
    TCTACGATTGAG
    CTGCACCTAGCTACGAGTGAG
    GTGCCCTAAGCATACTAGGCG
    GCCACTAGCAGAGCAGGGGGC
    CGACTCAGAGGTT
    CCACGCTAAGAGAAGTGCCGGAG
    GCACGTCCCTGAGCA
    GTCCATCGCAAAGCGTGGGGC
    GAGATCAGAGGGCCG
    TGGACGCGGGG
    GACCAGAGCCTCGCATGGGGG
    AGCGCGCGTG
    GCCGGTTGCTGTTCATTAGG
    ACCGACGGCAGCTAAAAGGG
    GACGCCGGGGAT
    CGACTCGCGCTTACAAGG

    If strict is True (default), the parser will raise a ValueError if the
    file contents does not strictly comply with the specified file format.
    """
    fmt = fmt.lower()
    if fmt == "minimal":
        import minimal
        return minimal.read(handle)
    else:
        raise ValueError(f"Unknown format {fmt}")

class Motif:
    name: str
    counts: Optional[matrix.FrequencyPositionMatrix]
    length: int
    alignment: Optional[pyobj]
    alphabet: str
    _pseudocounts: Dict[str, float]
    _background: Optional[Dict[str, float]]
    _mask: Optional[List[int]]

    def __init__(self, alphabet="ACGT", alignment: Optional[pyobj] =None, counts: Optional[Dict[str, List[float]]]=None):
        self.name = ""

        if counts is not None and alignment is not None:
            raise ValueError("Specify either counts or an alignment, don't specify both")
        elif counts is not None:
            self.alignment = None
            self.counts = matrix.FrequencyPositionMatrix(alphabet=alphabet, values=counts)
            self.length = self.counts.length
        elif alignment is not None:
            length = alignment.length
            frequencies: Dict[str, List[float]] = {}
            for letter in alphabet:
                if letter not in list(alignment.frequencies.keys()):
                    frequencies[letter] = [0.0 for _ in range(length)]
                else:
                    frequencies[letter] = [float(alignment.frequencies[letter][i]) for i in range(length)]
            self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
            self.alignment = alignment
            self.length = length
        else:
            self.counts = None
            self.alignment = None
            self.length = 0
        self.alphabet = alphabet
        self._pseudocounts = dict.fromkeys(self.alphabet, 0.0)
        self._background = None
        self.__set_background(value=None)
        self._mask = None
        self.__set_mask(None)
    
    def __init__(self, alphabet="ACGT", alignment: Optional[pyobj] =None, counts: Optional[matrix.FrequencyPositionMatrix]=None):
        self.name = ""

        if counts is not None and alignment is not None:
            raise ValueError("Specify either counts or an alignment, don't specify both")
        elif counts is not None:
            self.alignment = None
            self.counts = counts
            self.length = self.counts.length
        elif alignment is not None:
            length = alignment.length
            frequencies: Dict[str, List[float]] = {}
            for letter in alphabet:
                if letter not in list(alignment.frequencies.keys()):
                    frequencies[letter] = [0.0 for _ in range(length)]
                else:
                    frequencies[letter] = [float(alignment.frequencies[letter][i]) for i in range(length)]
            self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
            self.alignment = alignment
            self.length = length
        else:
            self.counts = None
            self.alignment = None
            self.length = 0
        self.alphabet = alphabet
        self._pseudocounts = dict.fromkeys(self.alphabet, 0.0)
        self._background = None
        self.__set_background(value=None)
        self._mask = None
        self.__set_mask(None)
    
    def __len__(self):
        return 0 if self.length is None else self.length
    
    def __format__(self, format_spec: Optional[str] = "default", **kwargs):
        if format_spec in ["pfm", "jaspar"]:
            motifs = [self]
            return utilities.jaspar_write(motifs, format_spec)
        elif format_spec == "transfac":
            motifs = [self]
            return utilities.transac_write(motifs)
        elif format_spec == "clusterbuster":
            motifs = [self]
            return utilities.clusterbuster_write(motifs, **kwargs)
        elif format_spec == "default":
            # Follow python convention and default to using __str__
            return str(self)
        else:
            raise ValueError(f"Unknown format type {format_spec}")
    
    def format(self, format_spec):
        return self.__format__(format_spec)
    
    @property
    def pseudocounts(self):
        return self._pseudocounts
    
    @pseudocounts.setter
    def pseudocounts(self, value):
        self._pseudocounts = {}
        if isinstance(value, Dict[str, float]):
            self._pseudocounts = {letter: value[letter] for letter in self.alphabet}
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = dict.fromkeys(self.alphabet, value)
    
    @property
    def background(self):
        return self._background
    
    def __set_background(self, value = None):
        if isinstance(value, Dict[str, float]): self._background = {letter: value[letter] for letter in self.alphabet}
        elif value is None: self._background = dict.fromkeys(self.alphabet, 1.0)
        elif isinstance(value, float) or isinstance(value, int):
            if not self._has_dna_alphabet() and not self._has_rna_alphabet():
                raise ValueError(
                    "Setting the background to a single value only works for DNA and RNA"
                    "motifs (in which case the value is interpreted as the GC content)"
                )
            value = float(value)
            T_or_U = "T" if self._has_dna_alphabet() else "U"
            d = {}
            d['A'] = (1.0 - value) / 2.0
            d["C"] = value / 2.0
            d["G"] = value / 2.0
            d[T_or_U] = (1.0 - value) / 2.0
            self._background = d
        
        # normalize
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total
    
    @background.setter
    def background(self, value = None):
        self.__set_background(value)

    def _has_dna_alphabet(self):
        return sorted(self.alphabet) == ["A", "C", "G", "T"]
    
    def _has_rna_alphabet(self):
        return sorted(self.alphabet) == ["A", "C", "G", "U"]
    
    @property
    def relative_entropy(self):
        background = self.background
        pseudocounts = self.pseudocounts
        alphabet = self.alphabet
        counts = self.counts
        length = self.length
        values = np.zeros(length)
        if self.alignment is None:
            total = np.array(
                [
                    sum(counts[c][i] + pseudocounts[c] for c in alphabet)
                    for i in range(length)
                ]
            )
            for letter in alphabet:
                frequencies = []
                for i in range(length):
                    frequencies.append(counts[letter][i])
                frequencies = np.array(frequencies) + pseudocounts[letter]
                mask = frequencies > 0
                frequencies = frequencies[mask] / total[mask]
                values[mask] += frequencies * np.log2(frequencies / background[letter])
        else:
            total = np.zeros(length)
            for letter in alphabet:
                frequencies = []
                for i in range(length):
                    frequencies.append(counts[letter][i])
                total += np.array(frequencies) + pseudocounts[letter]
            for letter in alphabet:
                frequencies = []
                for i in range(length):
                    frequencies.append(counts[letter][i])
                frequencies = np.array(frequencies) + pseudocounts[letter]
                mask = frequencies > 0
                frequencies = frequencies[mask] / total[mask]
                values[mask] += frequencies * np.log2(frequencies / background[letter])
        return values
    

    @property
    def pwm(self):
        """Calculate and return the position weight matrix for this motif."""
        return self.counts.normalize(self._pseudocounts)

    @property
    def pssm(self):
        """Calculate and return the position specific scoring matrix for this motif."""
        return self.pwm.log_odds(self._background)

    def __str__(self, masked=False):
        """Return string representation of a motif."""
        text = ""
        if self.alignment is not None:
            i = 0
            for m in self.alignment:
                if i != 0: text += '\n'
                text += m
                i += 1

        if masked:
            for i in range(self.length):
                if self._mask[i]:
                    text += "*"
                else:
                    text += " "
            text += "\n"
        return text
    
    @property
    def mask(self):
        return self._mask
    
    def __set_mask(self, mask = None):
        if self.length is None or self.length == 0:
            #  TODO:
            self._mask = []
        if mask is None:
            self._mask = [1 for _ in range(self.length)]
        elif len(mask) != self.length:
            raise ValueError(
                f"The length ({len(mask)}) of the mask is inconsistent with the length ({self.length}) of the motif"
            )
        elif isinstance(mask, str):
            _mask = []
            for char in mask:
                if char == "*":
                    _mask.append(1)
                elif char == " ":
                    _mask.append(0)
                else:
                    raise ValueError(
                        f"Mask should contain only '*' or ' ' and not a '{char}'"
                    )
            self._mask = _mask
        else:
            self._mask = [int(bool(c)) for c in mask]
    
    @mask.setter
    def mask(self, mask = None):
        self.__set_mask(mask)

    @property
    def consensus(self):
        """Return the consensus sequence."""
        return self.counts.consensus

    @property
    def anticonsensus(self):
        """Return the least probable pattern to be generated from this motif."""
        return self.counts.anticonsensus

    @property
    def degenerate_consensus(self):
        """Return the degenerate consensus sequence.

        Following the rules adapted from
        D. R. Cavener: "Comparison of the consensus sequence flanking
        translational start sites in Drosophila and vertebrates."
        Nucleic Acids Research 15(4): 1353-1361. (1987).

        The same rules are used by TRANSFAC.
        """
        return self.counts.degenerate_consensus

    def reverse_complement(self):
        """Return the reverse complement of the motif as a new motif."""
        alphabet = self.alphabet
        if not self._has_dna_alphabet() and not self._has_rna_alphabet():
            raise ValueError(
                "Calculating reverse complement only works for DNA and RNA motifs"
            )
        T_or_U = "T" if self._has_dna_alphabet() else "U"
        if self.alignment is not None:
            alignment = self.alignment.reverse_complement()
            sequences = alignment.sequences
            if T_or_U == "U":
                seqs = []
                for s in sequences: seqs.append(s.replace("T", "U"))
                res = Motif(alphabet=alphabet, alignment=Align.Alignment(sequences=seqs))
            else: res = Motif(alphabet=alphabet, alignment=alignment)
        else:  # has counts
            counts = {
                "A": self.counts[T_or_U][::-1],
                "C": self.counts["G"][::-1],
                "G": self.counts["C"][::-1],
                T_or_U: self.counts["A"][::-1],
            }
            res = Motif(alphabet=alphabet, counts=counts)
        res._mask = self._mask[::-1]
        res.background = {
            "A": self.background[T_or_U],
            "C": self.background["G"],
            "G": self.background["C"],
            T_or_U: self.background["A"],
        }
        res.pseudocounts = {
            "A": self.pseudocounts[T_or_U],
            "C": self.pseudocounts["G"],
            "G": self.pseudocounts["C"],
            T_or_U: self.pseudocounts["A"],
        }
        return res
    
    # Catch all calls to getitem where key is not a slice;
    # avoid compile-time errors
    def __getitem__(self, key):
        raise TypeError("motif indices must be slices")
    
    def __getitem__(self, key: slice):
        alphabet = self.alphabet
        if self.alignment is None:
            alignment = None
            if self.counts is None:
                counts = None
            else:
                temp = {letter: self.counts[letter][key] for letter in alphabet}
                if not isinstance(temp, Dict[str, List[float]]): raise TypeError()
                else: counts = temp
        else:
            alignment = self.alignment[:, key]
            counts = None
        motif = Motif(alphabet=alphabet, alignment=alignment, counts=counts)
        motif.mask = self.mask[key]
        if alignment is None and counts is None:
            try:
                length = self.length
            except AttributeError:
                pass
            else:
                motif.length = len(range(*key.indices(length)))
        motif.pseudocounts = self.pseudocounts.copy()
        motif.background = self.background.copy()
        return motif
    
    def weblogo(self, fname, fmt="PNG", **kwds):
        """Download and save a weblogo using the Berkeley weblogo service.

        Requires an internet connection.

        The parameters from ``**kwds`` are passed directly to the weblogo server.

        Currently, this method uses WebLogo version 3.3.
        These are the arguments and their default values passed to
        WebLogo 3.3; see their website at http://weblogo.threeplusone.com
        for more information::

            'stack_width' : 'medium',
            'stacks_per_line' : '40',
            'alphabet' : 'alphabet_dna',
            'ignore_lower_case' : True,
            'unit_name' : "bits",
            'first_index' : '1',
            'logo_start' : '1',
            'logo_end': str(self.length),
            'composition' : "comp_auto",
            'percentCG' : '',
            'scale_width' : True,
            'show_errorbars' : True,
            'logo_title' : '',
            'logo_label' : '',
            'show_xaxis': True,
            'xaxis_label': '',
            'show_yaxis': True,
            'yaxis_label': '',
            'yaxis_scale': 'auto',
            'yaxis_tic_interval' : '1.0',
            'show_ends' : True,
            'show_fineprint' : True,
            'color_scheme': 'color_auto',
            'symbols0': '',
            'symbols1': '',
            'symbols2': '',
            'symbols3': '',
            'symbols4': '',
            'color0': '',
            'color1': '',
            'color2': '',
            'color3': '',
            'color4': '',

        """
        if set(self.alphabet) == set("ACDEFGHIKLMNPQRSTVWY"):
            alpha = "alphabet_protein"
        elif set(self.alphabet) == set("ACGU"):
            alpha = "alphabet_rna"
        elif set(self.alphabet) == set("ACGT"):
            alpha = "alphabet_dna"
        else:
            alpha = "auto"

        frequencies = self.format("transfac")
        url = "https://weblogo.threeplusone.com/create.cgi"
        values = {
            "sequences": frequencies,
            "format": fmt.lower(),
            "stack_width": "medium",
            "stacks_per_line": "40",
            "alphabet": alpha,
            "ignore_lower_case": str(True),
            "unit_name": "bits",
            "first_index": "1",
            "logo_start": "1",
            "logo_end": str(self.length),
            "composition": "comp_auto",
            "percentCG": "",
            "scale_width": str(True),
            "show_errorbars": str(True),
            "logo_title": "",
            "logo_label": "",
            "show_xaxis": str(True),
            "xaxis_label": "",
            "show_yaxis": str(True),
            "yaxis_label": "",
            "yaxis_scale": "auto",
            "yaxis_tic_interval": "1.0",
            "show_ends": str(True),
            "show_fineprint": str(True),
            "color_scheme": "color_auto",
            "symbols0": "",
            "symbols1": "",
            "symbols2": "",
            "symbols3": "",
            "symbols4": "",
            "color0": "",
            "color1": "",
            "color2": "",
            "color3": "",
            "color4": "",
        }
        # values.update({k: "" if v is False else str(v) for k, v in kwds._asdict().items()})
        data = parse.urlencode(values).encode("utf-8")
        req = request.Request(url, data)
        response = request.urlopen(req)
        with open(fname, "wb") as f:
            im = response.read()
            f.write(im)

def write(motifs, fmt, **kwargs):
    """Return a string representation of motifs in the given format.

    Currently supported formats (case is ignored):
     - clusterbuster: Cluster Buster position frequency matrix format
     - pfm : JASPAR simple single Position Frequency Matrix
     - jaspar : JASPAR multiple PFM format
     - transfac : TRANSFAC like files

    """
    fmt = fmt.lower()
    if fmt in ("pfm", "jaspar"):
        return utilities.jaspar_write(motifs, fmt)
    elif fmt == "transfac":
        return utilities.transac_write(motifs)
    elif fmt == "clusterbuster":
        return utilities.clusterbuster_write(motifs, **kwargs)
    else:
        raise ValueError(f"Unknown format type {fmt}")