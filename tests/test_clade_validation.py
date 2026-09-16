import unittest

import numpy as np

from modules.gdm_clade_validation import (
    adjusted_rand_index,
    alignment_distance_matrix,
    clade_sequence_concordance,
    kmedoids,
)


ALIGNMENT = ">A\nAAAAAAAAAAAAAAAAAAAA\n>B\nAAAAAAAAAAAAAAAAAAAT\n>C\nAAAAAAAAAAAAAAAAAATA\n>D\nCCCCCCCCCCCCCCCCCCCC\n>E\nCCCCCCCCCCCCCCCCCCCA\n>F\nCCCCCCCCCCCCCCCCCCAC\n"


class CladeValidationTests(unittest.TestCase):
    def test_clear_two_group_family_is_ml_concordant(self):
        tree = '((A:0.1,B:0.1,C:0.1)0.99:0.5,(D:0.1,E:0.1,F:0.1)0.98:0.5);'
        result = clade_sequence_concordance(tree, ALIGNMENT)
        self.assertEqual(result['recommended_groups'], 2)
        self.assertGreaterEqual(result['ari'], 0.95)
        self.assertGreater(result['silhouette'], 0.5)
        self.assertTrue(str(result['status']).startswith('STRONG'))
        self.assertEqual(set(result['clade_table']['evidence_status']), {'ML-CONCORDANT'})

    def test_discordant_topology_is_not_called_verified(self):
        tree = '((A:0.1,D:0.1,E:0.1)0.99:0.5,(B:0.1,C:0.1,F:0.1)0.98:0.5);'
        result = clade_sequence_concordance(tree, ALIGNMENT, target_groups=2)
        self.assertLess(result['ari'], 0.5)
        self.assertIn('REVIEW', str(result['status']))

    def test_alignment_distance_and_kmedoids_are_deterministic(self):
        names, distance = alignment_distance_matrix(ALIGNMENT)
        labels1, medoids1 = kmedoids(distance, 2)
        labels2, medoids2 = kmedoids(distance, 2)
        self.assertEqual(names, ['A', 'B', 'C', 'D', 'E', 'F'])
        np.testing.assert_array_equal(labels1, labels2)
        self.assertEqual(medoids1, medoids2)
        self.assertGreaterEqual(adjusted_rand_index(labels1, labels2), 0.999)


if __name__ == '__main__':
    unittest.main()
