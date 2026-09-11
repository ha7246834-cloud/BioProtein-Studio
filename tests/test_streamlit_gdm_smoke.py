from pathlib import Path
import unittest

from streamlit.testing.v1 import AppTest


REPO_ROOT = Path(__file__).resolve().parents[1]
PAGE = REPO_ROOT / 'pages' / '2_Gene_Structure_Domain_Motif.py'


class GDMStreamlitSmokeTests(unittest.TestCase):
    def test_page_renders_without_runtime_exception(self):
        at = AppTest.from_file(str(PAGE), default_timeout=30)
        at.run()
        self.assertEqual(
            len(at.exception),
            0,
            msg='GDM Streamlit page raised a runtime exception during initial render.',
        )

    def test_core_controls_are_present(self):
        at = AppTest.from_file(str(PAGE), default_timeout=30)
        at.run()
        self.assertEqual(len(at.exception), 0)

        button_labels = [button.label for button in at.button]
        self.assertIn('🚀 Run complete analysis', button_labels)

        uploader_labels = [uploader.label for uploader in at.file_uploader]
        self.assertIn('Upload primary FASTA', uploader_labels)
        self.assertIn('Optional matching CDS FASTA', uploader_labels)
        self.assertIn('Optional matching genomic FASTA', uploader_labels)

        checkbox_labels = [box.label for box in at.checkbox]
        self.assertIn('Automatically build phylogeny when Newick is not uploaded', checkbox_labels)
        self.assertIn('Run NCBI CDD', checkbox_labels)
        self.assertIn('Run MEME', checkbox_labels)


if __name__ == '__main__':
    unittest.main()
