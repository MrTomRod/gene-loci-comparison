import os
from unittest import TestCase
from gene_loci_comparison import Locus
import matplotlib
import matplotlib.pyplot as plt
from bokeh.plotting import output_file, show, save

matplotlib.rcParams['font.family'] = "PT Sans Narrow"

save_plots = False

assert os.path.isfile('tests/test_loci.py'), f'Please set working directory to git root!'

pgap_file = 'tests/data/PGAP/FAM3257.gbk'
prokka_file = 'tests/data/Prokka/FAM3257.gbk'
new_prokka_file = 'tests/data/Prokka/Lbombicola_ESL0228.gbk'
bad_first_gene_file = 'tests/data/FirstGene/REFERENCE.gbk'
yeast_file = 'tests/data/yeast/R64-3-1.gbk'


class TestLocus(TestCase):
    def test_edge_genes(self):
        Locus(gbk_file=prokka_file, locus_tag='FAM3257_00001')  # at the start of the file
        Locus(gbk_file=prokka_file, locus_tag='FAM3257_03189')  # at the end of the file
        Locus(gbk_file=prokka_file, locus_tag='FAM3257_00099')  # at the end of scf1

    def test_nonexistent_gene(self):
        with self.assertRaises(KeyError):
            Locus(gbk_file=prokka_file, locus_tag='FAM3257_00000')

    def test_pgap(self):
        locus_tag = 'FAM3257_000993'

        locus = Locus(gbk_file=pgap_file, locus_tag=locus_tag)

        locus_tags = locus.locus_tags()
        print('locus_tags:', locus_tags)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_pgap.svg', format='svg')
        else:
            plt.show()

    def test_prokka(self):
        locus_tag = 'FAM3257_00934'

        locus = Locus(gbk_file=prokka_file, locus_tag=locus_tag)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_prokka.svg', format='svg')
        else:
            plt.show()

    def test_new_prokka(self):
        locus_tag = 'Lbombicola_ESL0228_00004'

        locus = Locus(gbk_file=new_prokka_file, locus_tag=locus_tag)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_new_prokka.svg', format='svg')
        else:
            plt.show()

    def test_gc(self):
        locus_tag = 'FAM3257_000993'

        locus = Locus(gbk_file=pgap_file, locus_tag=locus_tag)

        locus_tags = locus.locus_tags()
        print('locus_tags:', locus_tags)

        ax, _ = locus.plot_gc(figure_width=12, window_bp=20)

        if save_plots:
            plt.savefig('tests/output/locus/test_gc.svg', format='svg')
        else:
            plt.show()

    def test_to_string(self):
        locus_tag = 'FAM3257_000993'

        locus = Locus(gbk_file=pgap_file, locus_tag=locus_tag)

        svg_string = locus.plot_to_string(figure_width=12)

        with open('tests/output/locus/test_to_string.svg', 'w') as f:
            f.write(svg_string)

    def test_scf_end(self):
        locus_tag = 'FAM3257_00098'

        locus = Locus(gbk_file=prokka_file, locus_tag=locus_tag)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_scf_end.svg', format='svg')
        else:
            plt.show()

    def test_scf_start(self):
        locus_tag = 'FAM3257_00001'

        locus = Locus(gbk_file=prokka_file, locus_tag=locus_tag)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_scf_start.svg', format='svg')
        else:
            plt.show()

    def test_scf_start_bad_first_gene(self):
        locus = Locus(gbk_file=bad_first_gene_file, locus_tag='REFERENCE.1_000003')

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_scf_start.svg', format='svg')
        else:
            plt.show()

    def test_custom_colors(self):
        locus_tag = 'FAM3257_001019'

        locus = Locus(gbk_file=pgap_file, locus_tag=locus_tag)

        locus_tags = locus.locus_tags()
        print('locus_tags:', locus_tags)

        locus_to_color = dict(
            FAM3257_001014='#1271c3',
            FAM3257_001015='#3171c3',
            FAM3257_001016='#5d71c3',
            FAM3257_001017='#9371c3',
            FAM3257_001018='#b171c3',
            FAM3257_001019='#cb71c3',
            FAM3257_001020='#ea71c3',
            FAM3257_001021='#fd71c3',
            # FAM3257_001021='#fd71c3'  # last gene: white (default color)
        )

        locus.colorize(locus_to_color)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_custom_colors.svg', format='svg')
        else:
            plt.show()

    def test_custom_locus_tags(self):
        locus_tag = 'FAM3257_001019'

        locus = Locus(gbk_file=pgap_file, locus_tag=locus_tag)

        locus_tags = locus.locus_tags()
        print('locus_tags:', locus_tags)

        locus_to_new_name_dict = dict(
            FAM3257_001014='FAM3257_001014 | OG1',
            # FAM3257_001015='FAM3257_001015 | OG0',
            FAM3257_001016='FAM3257_001016 | OG4',
            FAM3257_001017='FAM3257_001017 | OG6',
            FAM3257_001018='FAM3257_001018 | OG8',
            FAM3257_001019='FAM3257_001019 \n OG5',  # NOTE: newlines (\n) do not work well.
            FAM3257_001020='FAM3257_001020 | OG7',
            FAM3257_001021='FAM3257_001021 | OG2',
            FAM3257_001022='FAM3257_001022 | OG2',
        )

        locus = locus.rename_labels(locus_to_new_name_dict)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_custom_locus_tags.svg', format='svg')
        else:
            plt.show()

    def test_bokeh(self):
        locus_tag = 'FAM3257_001020'

        locus = Locus(gbk_file=pgap_file, locus_tag=locus_tag, span=10000)

        bokeh = locus.plot_with_bokeh(figure_width=12, figure_height='auto', viewspan=3000)

        if save_plots:
            output_file(filename='tests/output/locus/test_single_locus_pgap.html', )
            save(bokeh)
        else:
            show(bokeh)

    def test_yeast(self):
        locus_tag = 'R64-3-1_00340'

        locus = Locus(gbk_file=yeast_file, locus_tag=locus_tag, span=5000)

        locus_tags = locus.locus_tags()
        print('locus_tags:', locus_tags)

        ax, _ = locus.plot(figure_width=12)
        ax.figure.tight_layout()

        if save_plots:
            plt.savefig('tests/output/locus/test_yeast.svg', format='svg')
        else:
            plt.show()
