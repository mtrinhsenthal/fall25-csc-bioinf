# Copyright 2008 by Bartek Wilczynski.  All rights reserved.
# Revisions copyright 2019 by Victor Lin.
# Adapted from test_Mymodule.py by Jeff Chang.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for motifs module."""

import math
import tempfile
import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingExternalDependencyError

    raise MissingExternalDependencyError(
        "Install numpy if you want to use Bio.motifs."
    ) from None

from Bio import motifs
from Bio.Seq import Seq


class TestBasic(unittest.TestCase):
    """Basic motif tests."""

    def test_format(self):
        m = motifs.create([Seq("ATATA")])
        m.name = "Foo"
        s1 = format(m, "pfm")
        expected_pfm = """  1.00   0.00   1.00   0.00  1.00
  0.00   0.00   0.00   0.00  0.00
  0.00   0.00   0.00   0.00  0.00
  0.00   1.00   0.00   1.00  0.00
"""
        s2 = format(m, "jaspar")
        expected_jaspar = """>None Foo
A [  1.00   0.00   1.00   0.00   1.00]
C [  0.00   0.00   0.00   0.00   0.00]
G [  0.00   0.00   0.00   0.00   0.00]
T [  0.00   1.00   0.00   1.00   0.00]
"""
        self.assertEqual(s2, expected_jaspar)
        s3 = format(m, "transfac")
        expected_transfac = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(s3, expected_transfac)
        self.assertRaises(ValueError, format, m, "foo_bar")

    def test_relative_entropy(self):
        m = motifs.create([Seq("ATATA"), Seq("ATCTA"), Seq("TTGTA")])
        self.assertEqual(len(m.alignment), 3)
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        self.assertEqual(m.pseudocounts, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
            )
        )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.8186697601117167,
                        1.7369655941662063,
                        0.5419780939258206,
                        1.7369655941662063,
                        1.7369655941662063,
                    ]
                ),
            )
        )
        m.background = None
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        pseudocounts = math.sqrt(len(m.alignment))
        m.pseudocounts = {
            letter: m.background[letter] * pseudocounts for letter in "ACGT"
        }
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.3532586861097656,
                        0.7170228827697498,
                        0.11859369972847714,
                        0.7170228827697498,
                        0.7170228827697499,
                    ]
                ),
            )
        )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.19727984803857979,
                        0.561044044698564,
                        0.20984910512125132,
                        0.561044044698564,
                        0.5610440446985638,
                    ]
                ),
            )
        )

    def test_reverse_complement(self):
        """Test if motifs can be reverse-complemented."""
        background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        pseudocounts = 0.5
        m = motifs.create([Seq("ATATA")])
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        expected_forward = """\
P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(received_forward, expected_forward)
        expected_forward_pwm = """\
        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.17   0.50   0.17   0.50   0.17
"""
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        expected_reverse = """\
P0      A      C      G      T
01      0      0      0      1      T
02      1      0      0      0      A
03      0      0      0      1      T
04      1      0      0      0      A
05      0      0      0      1      T
XX
//
"""
        self.assertEqual(received_reverse, expected_reverse)
        expected_reverse_pwm = """\
        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.50   0.17   0.50   0.17   0.50
"""
        self.assertEqual(str(m.pwm), expected_reverse_pwm)
        # Same but for RNA motif.
        background_rna = {"A": 0.3, "C": 0.2, "G": 0.2, "U": 0.3}
        pseudocounts = 0.5
        m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        expected_forward_rna_counts = """\
        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   0.00   1.00   0.00   1.00   0.00
"""
        self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
        expected_forward_rna_pwm = """\
        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.17   0.50   0.17   0.50   0.17
"""
        self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
        expected_reverse_rna_counts = """\
        0      1      2      3      4
A:   0.00   1.00   0.00   1.00   0.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   1.00   0.00   1.00   0.00   1.00
"""
        self.assertEqual(
            str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
        )
        expected_reverse_rna_pwm = """\
        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.50   0.17   0.50   0.17   0.50
"""
        self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)
        # Same thing, but now start with a motif calculated from a count matrix
        m = motifs.create([Seq("ATATA")])
        counts = m.counts
        m = motifs.Motif(counts=counts)
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        self.assertEqual(received_forward, expected_forward)
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        self.assertEqual(received_reverse, expected_reverse)
        self.assertEqual(str(m.pwm), expected_reverse_pwm)
        # Same, but for RNA count matrix
        m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
        counts = m_rna.counts
        m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
        self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
        self.assertEqual(
            str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
        )
        self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)