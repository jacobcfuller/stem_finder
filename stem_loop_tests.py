import unittest
import stem_loop as stem


class SearchTestCase(unittest.TestCase):
    """Tests for stem_loop.py
    """

    def test_complement(self):
        """Fails if nucleotides are not translated
        """
        self.assertEqual(stem.complement('ATUCGNatucgnKRSBDMYWVHkm'),
                         'TAAGCNtaagcnMYSVHKRWBDmk')

    def test_comp_error(self):
        """Fails if ValueError is not raised with unknown bases
        """
        with self.assertRaises(ValueError) as cm:
            stem.complement("X")
        self.assertEqual('Unknown bases passed: %s' % ', '.join("X"),
                         str(cm.exception))

    def test_compl_ratio(self):
        """Fails if correct ratio is not returned
        """
        seq = "AAAGCTTAGCGATAAGCTAA"
        self.assertEqual(stem.compl_ratio(seq), 0.75)

    def test_compl_ratio_error(self):
        """Fails if ValueError not raised on odd # input
        """
        seq = "AAAGCTTAGCGATAAGCTA"
        with self.assertRaises(ValueError) as cm:
            stem.compl_ratio(seq)
        self.assertEqual('Seq length must be even', str(cm.exception))

    def test_motif(self):
        """Fails if motif list is not returned
        """
        seq = "AATGGCGATAAAAATGGCGATAAG"
        motif = [(4, "TGGCGATA"), (16, "TGGCGATA")]
        self.assertEqual(stem.find_motif(seq), motif)

    def test_search_big(self):
        """Fails if doesn't find stem loop in input seq >= 100
        """
        seq = "GCCTGGAAAGGC"
        filler = "A"*50
        big_seq = filler + seq + filler
        motif = [(54, 'CTGGAAAG')]
        self.assertEqual(stem.search(motif, big_seq), [seq])

    def test_search_end(self):
        """Fails if doesn't find stem loop at the end of seq > 100
        """
        seq = "GCCTGGAAAGGC"
        filler = "A" * 50
        big_seq = filler*2 + seq
        motif = [(104, 'CTGGAAAG')]
        self.assertEqual(stem.search(motif, big_seq), [seq])

    def test_search_small(self):
        """Fails if doesn't find stem loop in input seq < 100
        """
        seq = "GCCTGGAAAGGC"
        motif = [(4, "CTGGAAAG")]
        self.assertEqual(stem.search(motif, seq), [seq])


if __name__ == '__main__':
    unittest.main()