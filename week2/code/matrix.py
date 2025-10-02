from typing import Dict, List, Optional, Union
from python import Bio.Seq as Seq
import math
import numpy as np
import _pwm
from python import numbers

class GenericPositionMatrix:
    alphabet: str
    length: int
    data: Dict[str, List[float]]

    def __init__(self, alphabet: str, values: Dict[str, List[int]]):
        self.data: Dict[str, List[float]] = Dict[str, List[float]]()
        self.alphabet = alphabet

        length = None
        for letter in alphabet:
            vals: List[float] = [float(v) for v in values[letter]]
            if length is None:
                self.length = len(vals)
            elif length != len(vals):
                raise Exception("data has inconsistent lengths")
            self.data[letter] = vals
    
    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        self.data: Dict[str, List[float]] = Dict[str, List[float]]()
        self.alphabet = alphabet

        length = None
        for letter in alphabet:
            vals: List[float] = [float(v) for v in values[letter]]  # normalize here
            if length is None:
                self.length = len(vals)
            elif length != len(vals):
                raise Exception("data has inconsistent lengths")
            self.data[letter] = vals

    def __str__(self):
        """Return a string containing nucleotides and counts of the alphabet in the Matrix."""
        words = [f"{i:6d}" for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            words = [
                f"{val:6.2f}" if not math.isinf(val) else f"{val:>6}" 
                for val in self.data[letter]
            ]

            line = f"{letter}: " + " ".join(words)
            lines.append(line)
        text = "\n".join(lines) + "\n"
        return text

    # Bypass __getitem__ for numeric computations
    # Always returns a float
    def get_value(self, letter: str, i: int) -> float:
        return self.data[letter][i]

    # Always returns a list of floats
    def get_column(self, i: int) -> List[float]:
        return [self.data[letter][i] for letter in self.alphabet]

    # Always returns a dict mapping letters -> float
    def get_column_dict(self, i: int, letters: List[str]) -> Dict[str, float]:
        return {letter: self.data[letter][i] for letter in letters}

    # Always returns a dict mapping letters -> list of floats
    def get_rows(self, indices: List[int], letters: List[str]) -> Dict[str, List[float]]:
        return {letter: [self.data[letter][j] for j in indices] for letter in letters}

    def __getitem__(self, key):
        if isinstance(key, tuple) and len(key) == 2:
            key1, key2 = key

            # First dimension
            if isinstance(key1, slice):
                indices1 = range(*key1.indices(len(self.alphabet)))
                letters1 = [self.alphabet[i] for i in indices1]
            elif isinstance(key1, int):
                letters1 = [self.alphabet[key1]]
            elif isinstance(key1, str) and len(key1) == 1:
                letters1 = [key1]
            else:
                raise KeyError(f"Cannot understand first key {key1}")

            # Second dimension
            if isinstance(key2, slice):
                indices2 = range(*key2.indices(self.length))
                return self.get_rows(list(indices2), letters1)
            elif isinstance(key2, int):
                return self.get_column_dict(key2, letters1)
            else:
                raise KeyError(f"Cannot understand second key {key2}")

        # Single-dimension keys
        elif isinstance(key, int):
            return self.data[self.alphabet[key]]
        elif isinstance(key, str) and len(key) == 1:
            return self.data[key]
        elif isinstance(key, slice):
            indices = range(*key.indices(len(self.alphabet)))
            letters = [self.alphabet[i] for i in indices]
            return {letter: self.data[letter] for letter in letters}
        else:
            raise KeyError(f"Unsupported key type: {key}")



    @property
    def consensus(self):
        sequence = ""
        for i in range(self.length):
            max_count = float('-inf')
            for letter in self.alphabet:
                count = self[letter][i]
                if count > max_count:
                    max_count = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq.Seq(sequence)
    
    @property
    def anticonsensus(self):
        sequence = ""
        for i in range(self.length):
            min_count = float('inf')
            for letter in self.alphabet:
                count = self[letter][i]
                if count < min_count:
                    min_count = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq.Seq(sequence)

    @property
    def degenerate_consensus(self):
        """Return the degenerate consensus sequence."""
        degenerate_nucleotide: dict[str, str] = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "U": "U",
            "AC": "M",
            "AG": "R",
            "AT": "W",
            "AU": "W",
            "CG": "S",
            "CT": "Y",
            "CU": "Y",
            "GT": "K",
            "GU": "K",
            "ACG": "V",
            "ACT": "H",
            "ACU": "H",
            "AGT": "D",
            "AGU": "D",
            "CGT": "B",
            "CGU": "B",
            "ACGT": "N",
            "ACGU": "N",
        }
        sequence = ""
        for i in range(self.length):

            def get(nucleotide):
                return self[nucleotide][i]  # noqa: B023

            nucleotides = sorted(self.data, key=get, reverse=True)

            counts = [self[c][i] for c in nucleotides]
            # Follow the Cavener rules:
            if counts[0] > sum(counts[1:]) and counts[0] > 2 * counts[1]:
                key = nucleotides[0]
            elif 4 * sum(counts[:2]) > 3 * sum(counts):
                key = "".join(sorted(nucleotides[:2]))
            elif counts[3] == 0:
                key = "".join(sorted(nucleotides[:3]))
            else:
                key = "ACGT"
            
            # nucleotide = degenerate_nucleotide.get(key, key)
            if key in degenerate_nucleotide:
                nucleotide = degenerate_nucleotide[key]
            else:
                nucleotide = key

            sequence += nucleotide
        return Seq.Seq(sequence)

    def calculate_consensus(self, substitution_matrix=None, plurality=None, identity=0, setcase=None):
        alphabet = self.alphabet
        if set(alphabet).union(set("ACGTUN-")) == set("ACGTUN-"):
            undefined = "N"
        else:
            undefined = "X"
        if substitution_matrix is None:
            if plurality is not None:
                raise ValueError(
                    "plurality must be None if substitution_matrix is None"
                )
            sequence = ""
            for i in range(self.length):
                maximum: float = 0.0
                total: float = 0.0
                for letter in alphabet:
                    count = self[letter][i]
                    total += count
                    if count > maximum:
                        maximum = count
                        consensus_letter = letter
                if maximum < identity * total:
                    consensus_letter = undefined
                else:
                    if setcase is None:
                        setcase_threshold = total / 2
                    else:
                        setcase_threshold = setcase * total
                    if maximum <= setcase_threshold:
                        consensus_letter = consensus_letter.lower()
                sequence += consensus_letter
        else:
            raise NotImplementedError(
                "calculate_consensus currently only supports substitution_matrix=None"
            )
        return sequence

    @property
    def gc_content(self):
        """Compute the fraction GC content."""
        alphabet = self.alphabet
        gc_total = 0.0
        total = 0.0
        for i in range(self.length):
            for letter in alphabet:
                if letter in "CG":
                    gc_total += self[letter][i]
                total += self[letter][i]
        return gc_total / total
    
    def reverse_complement(self):
        values = {}
        if self.alphabet == "ACGU":
            values["A"] = self["U"][::-1]
            values["U"] = self["A"][::-1]
        else:
            values["A"] = self["T"][::-1]
            values["T"] = self["A"][::-1]
        values["G"] = self["C"][::-1]
        values["C"] = self["G"][::-1]
        alphabet = self.alphabet
        return self.__class__(alphabet, values)

    def __getalphabet__(self):
        return self.alphabet
    
    def __getlength__(self):
        return self.length
    
    def __getdata__(self):
        return self.data

    def calculate(self, sequence: str):
        if sorted(self.alphabet) != ['A', 'C', 'G', 'T']:
            raise ValueError(f"PSSM has wrong alphabet: {self.alphabet} - Use only with DNA motifs")

        n = len(sequence)
        m = self.length

        scores = np.empty(n - m + 1, np.float32)
        logodds = np.array(
            [[self[letter][i] for letter in "ACGT"] for i in range(m)], float
        )
        _pwm.calculate(sequence, logodds, scores)
        return scores
    
class FrequencyPositionMatrix(GenericPositionMatrix):
    alphabet: str
    length: int

    def __init__(self, alphabet: str, values: Dict[str, List[int]]):
        super().__init__(alphabet=alphabet, values=values)
        self.length = super().__getlength__()
        self.alphabet = super().__getalphabet__()
    
    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        super().__init__(alphabet=alphabet, values=values)
        self.length = super().__getlength__()
        self.alphabet = super().__getalphabet__()

    def normalize(self, pseudocounts = None):
        counts: Dict[str, List[float]] = {}
        if pseudocounts is None:
            for letter in self.alphabet:
                counts[letter] = [0.0] * self.length
        elif isinstance(pseudocounts, dict[str, float]) or isinstance(pseudocounts, dict[str, int]):
            for letter in self.alphabet:
                counts[letter] = [float(pseudocounts[letter])] * self.length
        else:
            for letter in self.alphabet:
                counts[letter] = [float(pseudocounts)] * self.length
        for i in range(self.length):
            for letter in self.alphabet:
                counts[letter][i] += self[letter][i]
        # Actual normalization is done in the PositionWeightMatrix initializer
        return PositionWeightMatrix(self.alphabet, counts=counts)
    
class PositionWeightMatrix(GenericPositionMatrix):
    length: int
    alphabet: str

    def __init__(self, alphabet: str, counts: Dict[str, List[float]]):
        super().__init__(alphabet=alphabet, values=counts)
        self.length = super().__getlength__()
        self.alphabet = alphabet

        for i in range(self.length):
            total = sum(self[letter][i] for letter in alphabet)
            for letter in alphabet:
                self[letter][i] /= total
    
    def __init__(self, alphabet: str, counts: Dict[str, List[int]]):
        super().__init__(alphabet=alphabet, values=counts)
        self.length = super().__getlength__()
        self.alphabet = alphabet

        for i in range(self.length):
            total = sum(self[letter][i] for letter in alphabet)
            for letter in alphabet:
                self[letter][i] /= total
    
    def log_odds(self, background: Optional[Dict[str, float]]=None):
        values: Dict[str, List[float]] = {}
        alphabet = self.alphabet

        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())

        for letter in alphabet:
            background[letter] /= total
            values[letter] = []
        for i in range(self.length):
            for letter in alphabet:
                b = background[letter]

                if b > 0:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.log(p / b, 2)
                    else:
                        logodds = -math.inf
                else:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.inf
                    else:
                        logodds = math.nan
                values[letter].append(logodds)
        pssm = PositionSpecificScoringMatrix(alphabet=alphabet, values=values)
        return pssm

class PositionSpecificScoringMatrix(GenericPositionMatrix):
    alphabet: str
    length: int

    def __init__(self, alphabet: str, values: Dict[str, List[int]]):
        super().__init__(alphabet=alphabet, values=values)
        self.alphabet = alphabet
        self.length = super().__getlength__()

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        super().__init__(alphabet=alphabet, values=values)
        self.alphabet = alphabet
        self.length = super().__getlength__()    
    
    def calculate(self, sequence: str):
        if sorted(self.alphabet) != ['A', 'C', 'G', 'T']:
            raise ValueError(f"PSSM has wrong alphabet: {self.alphabet} - Use only with DNA motifs")

        n = len(sequence)
        m = self.length

        scores = np.empty(n - m + 1, np.float32)
        logodds = np.array(
            [[self[letter][i] for letter in "ACGT"] for i in range(m)], float
        )
        _pwm.calculate(sequence, logodds, scores)
        return scores

    def search(self, sequence: str, threshold: float=0.0, both: bool=True, chunksize: int=10**6):
        """Find hits with PWM score above given threshold.

        A generator function, returning found hits in the given sequence
        with the pwm score higher than the threshold.
        """
        sequence = sequence.upper()
        seq_len = len(sequence)
        motif_l = self.length
        chunk_starts = np.arange(0, seq_len, chunksize)
        if both:
            rc = self.reverse_complement()
        for chunk_start in chunk_starts:
            subseq = sequence[chunk_start : chunk_start + chunksize + motif_l - 1]
            pos_scores = self.calculate(subseq)
            pos_ind = pos_scores >= threshold
            pos_positions = np.where(pos_ind)[0] + chunk_start
            pos_scores = pos_scores[pos_ind]
            if both:
                neg_scores = rc.calculate(subseq)
                neg_ind = neg_scores >= threshold
                neg_positions = np.where(neg_ind)[0] + chunk_start
                neg_scores = neg_scores[neg_ind]
            else:
                neg_positions = np.empty((0), dtype=int)
                neg_scores = np.empty((0), dtype=np.float32)
            chunk_positions = np.append(pos_positions, neg_positions - seq_len)
            chunk_scores = np.append(pos_scores, neg_scores)
            order = np.argsort(np.append(pos_positions, neg_positions))
            chunk_positions = chunk_positions[order]
            chunk_scores = chunk_scores[order]
            yield from zip(chunk_positions, chunk_scores)

    @property
    def max(self):
        """Maximal possible score for this motif.

        returns the score computed for the consensus sequence.
        """
        score = 0.0
        letters = self.alphabet
        for position in range(self.length):
            score += max(self[letter][position] for letter in letters)
        return score

    @property
    def min(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
        score = 0.0
        letters = self.alphabet
        for position in range(self.length):
            score += min(self[letter][position] for letter in letters)
        return score

    @property
    def gc_content(self):
        """Compute the GC-ratio."""
        return super().gc_content

    def mean(self, background: Optional[Dict[str, float]]=None):
        """Return expected value of the score of a motif."""
        logodds: float

        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self.alphabet:
            background[letter] /= total
        sx = 0.0
        for i in range(self.length):
            for letter in self.alphabet:
                logodds = self.get_value(letter, i) # created a fn specifically for this, maybe not ideal
                if math.isnan(logodds):
                    continue
                if math.isinf(logodds) and logodds < 0:
                    continue
                b = background[letter]
                p = b * math.pow(2, logodds)
                sx += p * logodds
        return sx

    def std(self, background: Optional[Dict[str, float]]=None):
        """Return standard deviation of the score of a motif."""
        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self.alphabet:
            background[letter] /= total
        variance = 0.0
        for i in range(self.length):
            sx = 0.0
            sxx = 0.0
            for letter in self.alphabet:
                logodds = self.get_value(letter, i)
                if math.isnan(logodds):
                    continue
                if math.isinf(logodds) and logodds < 0:
                    continue
                b = background[letter]
                p = b * math.pow(2, logodds)
                sx += p * logodds
                sxx += p * logodds * logodds
            sxx -= sx * sx
            variance += sxx
        variance = max(variance, 0)  # to avoid roundoff problems
        return math.sqrt(variance)

    def dist_pearson(self, other):
        """Return the similarity score based on pearson correlation for the given motif against self.

        We use the Pearson's correlation of the respective probabilities.
        """
        if self.alphabet != other.alphabet:
            raise ValueError("Cannot compare motifs with different alphabets")

        max_p = -2.0
        for offset in range(-self.length + 1, other.length):
            if offset < 0:
                p = self.dist_pearson_at(other, -offset)
            else:  # offset>=0
                p = other.dist_pearson_at(self, offset)
            if max_p < p:
                max_p = p
                max_o = -offset
        return 1 - max_p, max_o

    def dist_pearson_at(self, other, offset):
        """Return the similarity score based on pearson correlation at the given offset."""
        letters = self.alphabet
        sx = 0.0  # \sum x
        sy = 0.0  # \sum y
        sxx = 0.0  # \sum x^2
        sxy = 0.0  # \sum x \cdot y
        syy = 0.0  # \sum y^2
        norm = max(self.length, offset + other.length) * len(letters)
        for pos in range(min(self.length - offset, other.length)):
            # xi = [self[letter, pos + offset] for letter in letters]
            # yi = [other[letter, pos] for letter in letters]
            xi = [self.get_value(letter, pos + offset) for letter in letters]
            yi = [other.get_value(letter, pos) for letter in letters]

            sx += sum(xi)
            sy += sum(yi)
            sxx += sum(x * x for x in xi)
            sxy += sum(x * y for x, y in zip(xi, yi))
            syy += sum(y * y for y in yi)
        sx /= norm
        sy /= norm
        sxx /= norm
        sxy /= norm
        syy /= norm
        numerator = sxy - sx * sy
        denominator = math.sqrt((sxx - sx * sx) * (syy - sy * sy))
        return numerator / denominator

    def distribution(self, background: Optional[Dict[str,float]]=None, precision=10**3):
        """Calculate the distribution of the scores at the given precision."""
        from .thresholds import ScoreDistribution

        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self.alphabet:
            background[letter] /= total
        return ScoreDistribution(precision=precision, pssm=self, background=background)
