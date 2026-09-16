import unittest

import matplotlib.pyplot as plt

from modules.gdm_circular import (
    _parse_tree,
    auto_topology_clades,
    circular_phylogeny_figure,
    support_text,
)


TREE = '(((GeneA:0.1,GeneB:0.1)0.95:0.2,(GeneC:0.1,GeneD:0.1)0.90:0.2)0.88:0.3,((GeneE:0.1,GeneF:0.1)0.97:0.2,(GeneG:0.1,GeneH:0.1)0.92:0.2)0.89:0.3);'


class CircularPhylogenyTests(unittest.TestCase):
    def test_auto_groups_partition_all_tips_once(self):
        tree = _parse_tree(TREE)
        groups, leaf_to_group = auto_topology_clades(tree, target_groups=4)
        self.assertEqual(len(groups), 4)
        self.assertEqual(set(leaf_to_group), {f'Gene{x}' for x in 'ABCDEFGH'})
        self.assertTrue(all(str(v).startswith('Auto Clade ') for v in leaf_to_group.values()))

    def test_generic_ids_are_preserved(self):
        tree = '(WRKY_001:0.1,(NAC-X:0.1,Protein_42:0.2)0.8:0.3);'
        fig, table = circular_phylogeny_figure(tree, target_groups=2, method_label='MAFFT + FastTree')
        text = ';'.join(table['member_ids'].astype(str))
        self.assertIn('WRKY_001', text)
        self.assertIn('NAC-X', text)
        self.assertIn('Protein_42', text)
        plt.close(fig)

    def test_support_formats_are_preserved(self):
        tree = _parse_tree('(A:0.1,B:0.1)0.95;')
        self.assertEqual(support_text(tree.root), '0.95')
        tree2 = _parse_tree('(A:0.1,B:0.1)95/100;')
        self.assertEqual(support_text(tree2.root), '95/100')

    def test_midpoint_and_subset_render(self):
        fig, table = circular_phylogeny_figure(
            TREE,
            order=['GeneA', 'GeneB', 'GeneE', 'GeneF'],
            target_groups=2,
            midpoint_root=True,
            show_support=True,
            method_label='Uploaded IQ-TREE tree/report',
        )
        members = ';'.join(table['member_ids'].astype(str))
        for wanted in ['GeneA', 'GeneB', 'GeneE', 'GeneF']:
            self.assertIn(wanted, members)
        for excluded in ['GeneC', 'GeneD', 'GeneG', 'GeneH']:
            self.assertNotIn(excluded, members)
        self.assertGreaterEqual(len(table), 2)
        plt.close(fig)


if __name__ == '__main__':
    unittest.main()
