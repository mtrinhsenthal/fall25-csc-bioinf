import math
import unittest
import numpy as np
import __init__ as motifs
from python import Bio.Seq as Seq
from typing import Tuple, NoneType
import time

tests: list[tuple[str, function[Tuple, NoneType]]] = []

def test(name: str):
    def decorator(fn: function[Tuple, NoneType]) -> function[Tuple, NoneType]:
        tests.append((name, fn))
        return fn
    return decorator

@test("test_format")
def test_format():
        tc = unittest.TestCase()
        m = motifs.create([Seq.Seq("ATATA")])
        m.name = "Foo"
        s1 = m.format(format_spec="pfm")

        expected_pfm = """  1.00   0.00   1.00   0.00   1.00
  0.00   0.00   0.00   0.00   0.00
  0.00   0.00   0.00   0.00   0.00
  0.00   1.00   0.00   1.00   0.00
"""

        s2 = m.format(format_spec="jaspar")
        expected_jaspar = """>None Foo
A [  1.00   0.00   1.00   0.00   1.00]
C [  0.00   0.00   0.00   0.00   0.00]
G [  0.00   0.00   0.00   0.00   0.00]
T [  0.00   1.00   0.00   1.00   0.00]
"""
        tc.assertEqual(s1, expected_pfm)
        tc.assertEqual(s2, expected_jaspar)
        tc.assertRaises(ValueError, lambda : m.format(format_spec="foo_bar"))

@test("test_format_transfac")
def test_format_transfac():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq("ATATA")])
    m.name = "Foo"
    s = m.format(format_spec="transfac")
    expected_transfac = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
    tc.assertEqual(s, expected_transfac)

@test("test_format_clusterbuster")
def test_format_clusterbuster():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq("ATATA")])
    m.name = "Foo"
    s = m.format(format_spec="clusterbuster")

    expected = """>Foo
1	0	0	0
0	0	0	1
1	0	0	0
0	0	0	1
1	0	0	0
"""
    tc.assertEqual(s, expected)

@test("test_relative_entropy_alignment")
def test_relative_entropy_alignment():
    m = motifs.create([Seq.Seq("ATATA"), Seq.Seq("ATCTA"), Seq.Seq("TTGTA")])
    tc = unittest.TestCase()
    tc.assertEqual(len(m.alignment), 3)
    tc.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    tc.assertEqual(m.pseudocounts, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})

    tc.assertTrue(
        np.allclose(
            m.relative_entropy,
            np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
        )
    )

    m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
    tc.assertTrue(
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
    tc.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    pseudocounts = math.sqrt(len(m.alignment))
    m.pseudocounts = {
        letter: m.background[letter] * pseudocounts for letter in "ACGT"
    }
    tc.assertTrue(
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
    tc.assertTrue(
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

@test("test_relative_entropy_counts")
def test_relative_entropy_counts():
    m = motifs.Motif(counts={'A': [1.0, 2.0, 3.0], 'C': [2.0, 3.0, 4.0], 'G': [1.0, 2.0, 3.0], 'T': [2.0, 3.0, 4.0]})
    tc = unittest.TestCase()
    tc.assertTrue(
        np.allclose(
            m.relative_entropy,
            np.array(
                [
                    0.08170417,
                    0.02904941,
                    0.01477186,
                ]
            ),
        )
    )

@test("test_pwm")
def test_pwm():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    expected = """        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
T:   0.00   1.00   0.00   1.00   0.00
"""
    tc.assertEqual(expected, m.pwm.__str__())

@test("test_pssm")
def test_pssm():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    expected="""        0      1      2      3      4
A:   2.00   -inf   2.00   -inf   2.00
C:   -inf   -inf   -inf   -inf   -inf
G:   -inf   -inf   -inf   -inf   -inf
T:   -inf   2.00   -inf   2.00   -inf
"""
    tc.assertEqual(expected, m.pssm.__str__())

@test("test_str")
def test_str():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    tc.assertEqual("ATATA", m.__str__())

    m.mask = "* * *"
    tc.assertEqual('ATATA* * *\n', m.__str__(masked=True))

@test("test_mask")
def test_mask():
    tc = unittest.TestCase()
    m = motifs.create([Seq.Seq(("ATATA"))])
    tc.assertEqual([1] * m.length, m.mask)

    m.mask = "* * *"
    tc.assertEqual([1, 0, 1, 0, 1], m.mask)

    m.mask = [2, 0, 3, 0, 1]
    tc.assertEqual([1, 0, 1, 0, 1], m.mask)

    def exception(): m.mask = "abcab"
    tc.assertRaises(ValueError, exception)

    def exception(): m.mask = [1,2]
    tc.assertRaises(ValueError, exception)

@test("test_reverse_complement")
def test_reverse_complement():
    tc = unittest.TestCase()
    background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
    pseudocounts = 0.5
    m = motifs.create([Seq.Seq(("ATATA"))])
    m.background = background
    m.pseudocounts = pseudocounts

    received_forward = m.format(format_spec="transfac")
    expected_forward = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
    tc.assertEqual(received_forward, expected_forward)

    expected_forward_pwm = """        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.17   0.50   0.17   0.50   0.17
"""
    tc.assertEqual(str(m.pwm), expected_forward_pwm)
    
    rc = m.reverse_complement()

    received_reverse = rc.format(format_spec="transfac")
    expected_reverse = """P0      A      C      G      T
01      0      0      0      1      T
02      1      0      0      0      A
03      0      0      0      1      T
04      1      0      0      0      A
05      0      0      0      1      T
XX
//
"""
    tc.assertEqual(received_reverse, expected_reverse)

    expected_reverse_pwm = """        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.50   0.17   0.50   0.17   0.50
"""
    tc.assertEqual(expected_reverse_pwm, str(rc.pwm))


    background_rna = {"A": 0.3, "C": 0.2, "G": 0.2, "U": 0.3}
    pseudocounts = 0.5
    m_rna = motifs.create([Seq.Seq("AUAUA")], alphabet="ACGU")
    m_rna.background = background_rna
    m_rna.pseudocounts = pseudocounts
    expected_forward_rna_counts = """        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   0.00   1.00   0.00   1.00   0.00
"""
    tc.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
    
    expected_forward_rna_pwm = """        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.17   0.50   0.17   0.50   0.17
"""
    tc.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
    expected_reverse_rna_counts = """        0      1      2      3      4
A:   0.00   1.00   0.00   1.00   0.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   1.00   0.00   1.00   0.00   1.00
"""
    tc.assertEqual(
        str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
    )
    expected_reverse_rna_pwm = """        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.50   0.17   0.50   0.17   0.50
"""
    tc.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

    m = motifs.create([Seq.Seq("ATATA")])
    counts = m.counts
    m = motifs.Motif(counts=counts)
    m.background = background
    m.pseudocounts = pseudocounts
    received_forward = m.format(format_spec="transfac")
    tc.assertEqual(received_forward, expected_forward)
    tc.assertEqual(str(m.pwm), expected_forward_pwm)
    m = m.reverse_complement()
    received_reverse = m.format(format_spec="transfac")
    tc.assertEqual(received_reverse, expected_reverse)
    tc.assertEqual(str(m.pwm), expected_reverse_pwm)
    # Same, but for RNA count matrix
    m_rna = motifs.create([Seq.Seq("AUAUA")], alphabet="ACGU")
    counts = m_rna.counts
    m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
    m_rna.background = background_rna
    m_rna.pseudocounts = pseudocounts
    tc.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
    tc.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
    tc.assertEqual(
        str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
    )
    tc.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)

@test("test_consensus")
def test_consensus():
    tc = unittest.TestCase()
    m = motifs.create([
        Seq.Seq("ATATA"),
        Seq.Seq("ATCTA"),
        Seq.Seq("TTGTA"),
        Seq.Seq("ATATA"),
    ])

    expected_consensus = "ATATA"
    tc.assertEqual(str(m.consensus), expected_consensus)

@test("test_anticonsensus")
def test_anticonsensus():
    tc = unittest.TestCase()
    m = motifs.create([
        Seq.Seq("ATCGA"),
        Seq.Seq("ATCGA"),
        Seq.Seq("GGGTG"),
        Seq.Seq("GGGTG"),
        Seq.Seq("CCACC"),
        Seq.Seq("CCACC"),
        Seq.Seq("TATAT")
    ])

    expected_anticonsensus = "TATAT"
    tc.assertEqual(str(m.anticonsensus), expected_anticonsensus)

@test("test_degenerate_consensus")
def test_degenerate_consensus():
    tc = unittest.TestCase()
    m = motifs.create([
        Seq.Seq("ATATA"),
        Seq.Seq("ATCTA"),
        Seq.Seq("TTGTA"),
        Seq.Seq("ATGTA")
    ])
    # Position-wise breakdown:
    # 1: A (3), T (1) → A
    # 2: T (4) → T
    # 3: A (1), C (1), G (2) → V (A/C/G)
    # 4: T (4) → T
    # 5: A (4) → A

    expected_degenerate_consensus = "ATVTA"
    tc.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)

@test("test_degenerate_consensus_with_ties")
def test_degenerate_consensus_with_ties():
    tc = unittest.TestCase()
    m = motifs.create([
        Seq.Seq("A"),
        Seq.Seq("C"),
        Seq.Seq("G"),
        Seq.Seq("T"),
    ])
    expected_degenerate_consensus = "N"  # all bases equally represented
    tc.assertEqual(str(m.degenerate_consensus), expected_degenerate_consensus)

@test("test_degenerate_consensus_rna")
def test_degenerate_consensus_rna():
    tc = unittest.TestCase()
    m_rna = motifs.create([
        Seq.Seq("AUAUA"),
        Seq.Seq("AUCUA"),
        Seq.Seq("UUGUA"),
        Seq.Seq("AUGUA")
    ], alphabet="ACGU")

    expected_degenerate_consensus = "AUVUA"
    tc.assertEqual(str(m_rna.degenerate_consensus), expected_degenerate_consensus)

@test("test_getitem")
def test_getitem():
    tc = unittest.TestCase()
    m = motifs.create(["AACGCCA", "ACCGCCC", "AACTCCG"])
    expected="""AACGCC\nACCGCC\nAACTCC"""
    tc.assertEqual(str(m[:-1]), expected)

    expected="\n\n"
    tc.assertEqual(str(m[0:0]), expected)

    expected="A\nA\nA"
    tc.assertEqual(str(m[0:1]), expected)
    
    # EXCEPTIONS
    tc.assertRaises(TypeError, lambda: m[0])
    tc.assertRaises(TypeError, lambda: m['A'])

@test("test_minimal_parser_1")
def test_minimal_parser_1():
    tc = unittest.TestCase()
    """Parse motifs/minimal_test.meme file."""
    with open("week2/tests/minimal_test.meme") as stream:
        record = motifs.parse(stream, "minimal")
    tc.assertEqual(record.version, "4")
    tc.assertEqual(record.alphabet, "ACGT")
    tc.assertEqual(len(record.sequences), 0)
    tc.assertEqual(record.command, "")
    tc.assertEqual(len(record), 3)
    motif = record[0]
    tc.assertEqual(motif.name, "KRP")
    # tc.assertEqual(record["KRP"], motif)
    # tc.assertEqual(motif.num_occurrences, 17)
    tc.assertEqual(motif.length, 19)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 4.1e-09, places=10)
    tc.assertEqual(motif.alphabet, "ACGT")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "TGTGATCGAGGTCACACTT")
    tc.assertEqual(motif.degenerate_consensus, "TGTGANNNWGNTCACAYWW")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    1.1684297174927525,
                    0.9432809925744818,
                    1.4307101633876265,
                    1.1549413780465179,
                    0.9308256303218774,
                    0.009164393966550805,
                    0.20124190687894253,
                    0.17618542656995528,
                    0.36777933103380855,
                    0.6635834532368525,
                    0.07729943368061855,
                    0.9838293592717438,
                    1.72489868427398,
                    0.8397561713453014,
                    1.72489868427398,
                    0.8455332015343343,
                    0.3106481207768122,
                    0.7382733641762232,
                    0.537435993300495,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "TGATCGA")
    motif = record[1]
    tc.assertEqual(motif.name, "IFXA")
    # tc.assertEqual(record["IFXA"], motif)
    # tc.assertEqual(motif.num_occurrences, 14)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 3.2e-35, places=36)
    tc.assertEqual(motif.alphabet, "ACGT")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "TACTGTATATATATCCAG")
    tc.assertEqual(motif.degenerate_consensus, "TACTGTATATAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.9632889858595118,
                    1.02677956765017,
                    2.451526420551951,
                    1.7098384161433415,
                    2.2598671267551107,
                    1.7098384161433415,
                    1.02677956765017,
                    1.391583804103081,
                    1.02677956765017,
                    1.1201961888781142,
                    0.27822438781180836,
                    0.36915366971717867,
                    1.7240522753630425,
                    0.3802185945622609,
                    0.790937683007783,
                    2.451526420551951,
                    1.7240522753630425,
                    1.3924085743645374,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CTGTATA")

@test("test_minimal_parser_2")
def test_minimal_parser_2():
    tc = unittest.TestCase()
    with open("week2/tests/minimal_test.meme") as stream:
        record = motifs.parse(stream, "minimal")
    motif = record[2]
    tc.assertEqual(motif.name, "IFXA_no_nsites_no_evalue")
    # tc.assertEqual(record["IFXA_no_nsites_no_evalue"], motif)
    # tc.assertEqual(motif.num_occurrences, 20)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 0.0, places=36)
    tc.assertEqual(motif.alphabet, "ACGT")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "TACTGTATATATATCCAG")
    tc.assertEqual(motif.degenerate_consensus, "TACTGTATATAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.99075309,
                    1.16078104,
                    2.45152642,
                    1.70983842,
                    2.25986713,
                    1.70983842,
                    1.16078104,
                    1.46052586,
                    1.16078104,
                    1.10213019,
                    0.29911041,
                    0.36915367,
                    1.72405228,
                    0.37696488,
                    0.85258086,
                    2.45152642,
                    1.72405228,
                    1.42793329,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CTGTATA")

@test("test_minimal_parser_rna")
def test_minimal_parser_rna():
    tc = unittest.TestCase()
    """Test if Bio.motifs can parse MEME output files using RNA."""
    with open("week2/tests/minimal_test_rna.meme") as stream:
        record = motifs.parse(stream, "minimal")
    tc.assertEqual(record.version, "4")
    tc.assertEqual(record.alphabet, "ACGU")
    tc.assertEqual(len(record.sequences), 0)
    tc.assertEqual(record.command, "")
    tc.assertEqual(len(record), 3)
    motif = record[0]
    tc.assertEqual(motif.name, "KRP_fake_RNA")
    # tc.assertEqual(record["KRP_fake_RNA"], motif)
    # tc.assertEqual(motif.num_occurrences, 17)
    tc.assertEqual(motif.length, 19)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 4.1e-09, places=10)
    tc.assertEqual(motif.alphabet, "ACGU")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "UGUGAUCGAGGUCACACUU")
    tc.assertEqual(motif.degenerate_consensus, "UGUGANNNWGNUCACAYWW")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    1.1684297174927525,
                    0.9432809925744818,
                    1.4307101633876265,
                    1.1549413780465179,
                    0.9308256303218774,
                    0.009164393966550805,
                    0.20124190687894253,
                    0.17618542656995528,
                    0.36777933103380855,
                    0.6635834532368525,
                    0.07729943368061855,
                    0.9838293592717438,
                    1.72489868427398,
                    0.8397561713453014,
                    1.72489868427398,
                    0.8455332015343343,
                    0.3106481207768122,
                    0.7382733641762232,
                    0.537435993300495,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "UGAUCGA")
    motif = record[1]
    tc.assertEqual(motif.name, "IFXA_fake_RNA")
    # tc.assertEqual(record["IFXA_fake_RNA"], motif)
    # tc.assertEqual(motif.num_occurrences, 14)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 3.2e-35, places=36)
    tc.assertEqual(motif.alphabet, "ACGU")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "UACUGUAUAUAUAUCCAG")
    tc.assertEqual(motif.degenerate_consensus, "UACUGUAUAUAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.9632889858595118,
                    1.02677956765017,
                    2.451526420551951,
                    1.7098384161433415,
                    2.2598671267551107,
                    1.7098384161433415,
                    1.02677956765017,
                    1.391583804103081,
                    1.02677956765017,
                    1.1201961888781142,
                    0.27822438781180836,
                    0.36915366971717867,
                    1.7240522753630425,
                    0.3802185945622609,
                    0.790937683007783,
                    2.451526420551951,
                    1.7240522753630425,
                    1.3924085743645374,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CUGUAUA")

    motif = record[2]
    tc.assertEqual(motif.name, "IFXA_no_nsites_no_evalue_fake_RNA")
    # tc.assertEqual(record["IFXA_no_nsites_no_evalue_fake_RNA"], motif)
    # tc.assertEqual(motif.num_occurrences, 20)
    tc.assertEqual(motif.length, 18)
    tc.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
    tc.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
    tc.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
    tc.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
    # tc.assertAlmostEqual(motif.evalue, 0.0, places=36)
    tc.assertEqual(motif.alphabet, "ACGU")
    tc.assertIsNone(motif.alignment)
    tc.assertEqual(motif.consensus, "UACUGUAUAUAUAUCCAG")
    tc.assertEqual(motif.degenerate_consensus, "UACUGUAUAUAHAWMCAG")
    tc.assertTrue(
        np.allclose(
            motif.relative_entropy,
            np.array(
                [
                    0.99075309,
                    1.16078104,
                    2.45152642,
                    1.70983842,
                    2.25986713,
                    1.70983842,
                    1.16078104,
                    1.46052586,
                    1.16078104,
                    1.10213019,
                    0.29911041,
                    0.36915367,
                    1.72405228,
                    0.37696488,
                    0.85258086,
                    2.45152642,
                    1.72405228,
                    1.42793329,
                ]
            ),
        )
    )
    tc.assertEqual(motif[2:9].consensus, "CUGUAUA")

@test("test_pwm_getitem")
def test_pwm_getitem():
    tc = unittest.TestCase()
    counts_ = {'A': [2.0, 9.0, 0.0, 1.0, 32.0, 3.0, 46.0, 1.0, 43.0, 15.0, 2.0, 2.0], 'C': [1.0, 33.0, 45.0, 45.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0], 'G': [39.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 43.0], 'T': [4.0, 2.0, 0.0, 0.0, 13.0, 42.0, 0.0, 45.0, 3.0, 30.0, 0.0, 0.0]}
    m = motifs.Motif(counts=counts_)
    counts = m.counts
    ints = range(13)
    i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12 = ints
    # slice, slice
    d = counts[i1::i2, i2:i12:i3]
    tc.assertEqual(len(d), 2)
    tc.assertEqual(len(d["C"]), 4)
    tc.assertEqual(len(d["T"]), 4)
    tc.assertAlmostEqual(d["C"][i0], 45.0)
    tc.assertAlmostEqual(d["C"][i1], 1.0)
    tc.assertAlmostEqual(d["C"][i2], 0.0)
    tc.assertAlmostEqual(d["C"][i3], 1.0)
    tc.assertAlmostEqual(d["T"][i0], 0.0)
    tc.assertAlmostEqual(d["T"][i1], 42.0)
    tc.assertAlmostEqual(d["T"][i2], 3.0)
    tc.assertAlmostEqual(d["T"][i3], 0.0)
    #    slice, int
    d = counts[i1::i2, i4]
    tc.assertEqual(len(d), 2)
    tc.assertAlmostEqual(d["C"], 1.0)
    tc.assertAlmostEqual(d["T"], 13.0)
    
@test("test_pwm_mixed")
def test_pwm_mixed():
    tc = unittest.TestCase()
    counts_ = {'A': [2.0, 9.0, 0.0, 1.0, 32.0, 3.0, 46.0, 1.0, 43.0, 15.0, 2.0, 2.0], 'C': [1.0, 33.0, 45.0, 45.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0], 'G': [39.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 43.0], 'T': [4.0, 2.0, 0.0, 0.0, 13.0, 42.0, 0.0, 45.0, 3.0, 30.0, 0.0, 0.0]}
    m = motifs.Motif(counts=counts_)
    counts = m.counts
    pwm = counts.normalize(pseudocounts=0.25)
    pssm = pwm.log_odds()
    result = pssm.calculate(str(Seq.Seq("AcGTgTGCGtaGTGCGT")))
    tc.assertEqual(6, len(result))
    tc.assertAlmostEqual(float(result[0]), -29.18363571, places=5)
    tc.assertAlmostEqual(float(result[1]), -38.3365097, places=5)
    tc.assertAlmostEqual(float(result[2]), -29.17756271, places=5)
    tc.assertAlmostEqual(float(result[3]), -38.04542542, places=5)
    tc.assertAlmostEqual(float(result[4]), -20.3014183, places=5)
    tc.assertAlmostEqual(float(result[5]), -25.18009186, places=5)

@test("test_pwm_simple")
def test_pwm_simple():
    tc = unittest.TestCase()
    counts = {'A': [2.0, 9.0, 0.0, 1.0, 32.0, 3.0, 46.0, 1.0, 43.0, 15.0, 2.0, 2.0], 'C': [1.0, 33.0, 45.0, 45.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0], 'G': [39.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 43.0], 'T': [4.0, 2.0, 0.0, 0.0, 13.0, 42.0, 0.0, 45.0, 3.0, 30.0, 0.0, 0.0]}
    m = motifs.Motif(counts=counts)
    s = str(Seq.Seq("ACGTGTGCGTAGTGCGT"))
    pwm = m.counts.normalize(pseudocounts=0.25)
    pssm = pwm.log_odds()
    result = pssm.calculate(s)
    tc.assertEqual(6, len(result))
    # The fast C-code in Bio/motifs/_pwm.c stores all results as 32-bit
    # floats; the slower Python code in Bio/motifs/__init__.py uses 64-bit
    # doubles. The C-code and Python code results will therefore not be
    # exactly equal. Test the first 5 decimal places only to avoid either
    # the C-code or the Python code to inadvertently fail this test.
    tc.assertAlmostEqual(float(result[0]), -29.18363571, places=5)
    tc.assertAlmostEqual(float(result[1]), -38.3365097, places=5)
    tc.assertAlmostEqual(float(result[2]), -29.17756271, places=5)
    tc.assertAlmostEqual(float(result[3]), -38.04542542, places=5)
    tc.assertAlmostEqual(float(result[4]), -20.3014183, places=5)
    tc.assertAlmostEqual(float(result[5]), -25.18009186, places=5)


def testRunner():
    passed = total = 0
    start = time.perf_counter()
    for name, test in tests:
        try:
            test()
            print(f"{name} ... ok")
            passed += 1
        except AssertionError:
            print(f"{name} ... not ok")
        except:
            print(f"{name} ... ERROR")
        total += 1
    end = time.perf_counter()
    elapsed_time = end - start
    print(f"\n\n{passed} out of {total} tests passed; tests ran in {elapsed_time:.6f} seconds\n")

testRunner()
