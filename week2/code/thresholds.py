"""Approximate calculation of appropriate thresholds for motif finding."""

# only used in PositionSpecificScoringMatrix class (matrix.py) 
class ScoreDistribution:
    """Class representing approximate score distribution for a given motif.

    Utilizes a dynamic programming approach to calculate the distribution of
    scores with a predefined precision. Provides a number of methods for calculating
    thresholds for motif occurrences.
    """
    min_score: float
    interval: float
    n_points: int
    ic: float
    step: float
    mo_density: List[float]
    bg_density: List[float]

    def __init__(self, motif=None, precision: int=10**3, pssm=None, background: Optional[Dict[str, float]] =None):
        """Initialize the class."""
        
        # i think this can be removed, since we only ever call this from Matrix.distribution, which always passes pssm and never motif
        if pssm is None:
            assert motif is not None, "motif must be provided if pssm is None"
            self.min_score = min(0.0, motif.min_score())
            self.interval = max(0.0, motif.max_score()) - self.min_score
            self.n_points = precision * motif.length
            self.ic = motif.ic()
        else:
            assert background is not None, "background must be provided if pssm is not None"
            self.min_score = min(0.0, pssm.min) # minimal possible score for this motif
            self.interval = max(0.0, pssm.max) - self.min_score # the range of possible scores for this motif
            self.n_points = precision * pssm.length
            self.ic = pssm.mean(background) # information content (how much a motif deviates from random)

        # initialize 2 probability distributions
        # mo_density: for motif containing sequences
        # bg_density: for background sequences
        self.step = self.interval / (self.n_points - 1)
        self.mo_density = [0.0] * self.n_points
        self.mo_density[-self._index_diff(self.min_score)] = 1.0 # sets the bin corresponding to the min possible score to have a probability of 1.0
        self.bg_density = [0.0] * self.n_points
        self.bg_density[-self._index_diff(self.min_score)] = 1.0
        if pssm is None: # pssm should never be None in our case -- so I think we do not need to implement modify
            for lo, mo in zip(motif.log_odds(), motif.pwm()):
                self.modify(lo, mo, motif.background)
        else:
            # loop through each position in motif, create distributions for this position, and update existing distributions
            for position in range(pssm.length):
                mo_new = [0.0] * self.n_points
                bg_new = [0.0] * self.n_points
                # lo = pssm[:, position]
                lo = pssm.get_column_dict(position, pssm.alphabet)

                for letter, score in lo.items():
                    bg = background[letter]
                    # mo = pow(2, pssm[letter, position]) * bg
                    mo = pow(2, pssm.get_value(letter, position)) * bg

                    numeric_score = pssm.get_value(letter, position)  # # use get_value for numeric access to the PSSM, always float
                    d = self._index_diff(numeric_score)
                    # d = 1
                    # d = self._index_diff(score)

                    for i in range(self.n_points):
                        mo_new[self._add(i, d)] += self.mo_density[i] * mo
                        bg_new[self._add(i, d)] += self.bg_density[i] * bg
                self.mo_density = mo_new
                self.bg_density = bg_new

    def _index_diff(self, x: float, y: float=0.0) -> int:
        return int((x - y + 0.5 * self.step) // self.step)

    def _add(self, i: int, j: int) -> int:
        return max(0, min(self.n_points - 1, i + j))

    # I think we can get rid of this, since ScoreDistribution is only ever called with pssm
    def modify(self, scores: Dict[str, float], mo_probs: Dict[str, float], bg_probs: Dict[str, float]) -> None:
        """Modify motifs and background density."""
        mo_new = [0.0] * self.n_points
        bg_new = [0.0] * self.n_points
        for k, v in scores.items():
            d = self._index_diff(v)
            for i in range(self.n_points):
                mo_new[self._add(i, d)] += self.mo_density[i] * mo_probs[k]
                bg_new[self._add(i, d)] += self.bg_density[i] * bg_probs[k]
        self.mo_density = mo_new
        self.bg_density = bg_new

    def threshold_fpr(self, fpr: float) -> float:
        """Approximate the log-odds threshold which makes the type I error (false positive rate)."""
        i = self.n_points
        prob = 0.0
        while prob < fpr:
            i -= 1
            prob += self.bg_density[i]
        return self.min_score + i * self.step

    def threshold_fnr(self, fnr: float) -> float:
        """Approximate the log-odds threshold which makes the type II error (false negative rate)."""
        i = -1
        prob = 0.0
        while prob < fnr:
            i += 1
            prob += self.mo_density[i]
        return self.min_score + i * self.step

    def threshold_balanced(self, rate_proportion: float=1.0, return_rate: bool=False):
        """Approximate log-odds threshold making FNR equal to FPR times rate_proportion."""
        i = self.n_points
        fpr = 0.0
        fnr = 1.0
        while fpr * rate_proportion < fnr:
            i -= 1
            fpr += self.bg_density[i]
            fnr -= self.mo_density[i]
        if return_rate:
            return self.min_score + i * self.step, fpr
        else:
            return self.min_score + i * self.step

    def threshold_patser(self):
        """Threshold selection mimicking the behaviour of patser (Hertz, Stormo 1999) software.

        It selects such a threshold that the log(fpr)=-ic(M)
        note: the actual patser software uses natural logarithms instead of log_2, so the numbers
        are not directly comparable.
        """
        return self.threshold_fpr(fpr=2**-self.ic)