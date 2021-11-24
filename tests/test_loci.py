import os
from unittest import TestCase
from gene_loci_comparison import Loci
import matplotlib
from bokeh.layouts import column
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
    def test_multiple(self):
        loci_of_interest = [
            dict(gbk=prokka_file, gene='FAM3257_00934', title='title1'),
            dict(gbk=pgap_file, gene='FAM3257_000019', title='title2'),
            dict(gbk=pgap_file, gene='FAM3257_001020', title='title3'),
        ]

        locus_to_color_dict = {locus['gene']: '#1984ff' for locus in loci_of_interest}

        loci = Loci.generate(
            loci_of_interest,
            locus_to_color_dict=locus_to_color_dict
        )

        print('locus tags:', loci.locus_tags())

        plot = loci.plot(auto_reverse=False)

        if save_plots:
            plot.savefig('tests/output/loci/test_multiple.svg', format='svg')
        else:
            plot.show()

    def test_multiple_auto_reverse(self):
        # bottom plot should be reversed.
        loci_of_interest = [
            dict(gbk=prokka_file, gene='FAM3257_00934', title='prokka'),
            dict(gbk=pgap_file, gene='FAM3257_000019', title='pgap'),
            dict(gbk=pgap_file, gene='FAM3257_001020', title='pgap'),
        ]

        locus_to_color_dict = {locus['gene']: '#1984ff' for locus in loci_of_interest}

        loci = Loci.generate(
            loci_of_interest,
            locus_to_color_dict=locus_to_color_dict,
            span=4000
        )

        print('locus tags:', loci.locus_tags())

        plot = loci.plot(
            auto_reverse=True,
            annotate_inline=False,  # show labels above gene boxes
        )

        if save_plots:
            plot.savefig('tests/output/loci/test_multiple_auto_reverse.svg', format='svg')
        else:
            plot.show()

    def test_multiple_auto_reverse_gc(self):
        # bottom plot should be reversed.
        loci_of_interest = [
            dict(gbk=prokka_file, gene='FAM3257_00934', title='prokka'),
            dict(gbk=pgap_file, gene='FAM3257_000019', title='pgap'),
            dict(gbk=pgap_file, gene='FAM3257_001020', title='pgap'),
        ]

        locus_to_color_dict = {locus['gene']: '#1984ff' for locus in loci_of_interest}

        loci = Loci.generate(
            loci_of_interest,
            locus_to_color_dict=locus_to_color_dict,
            span=4000
        )

        print('locus tags:', loci.locus_tags())

        plot = loci.plot_gc(
            auto_reverse=True,
            window_bp=200
        )

        if save_plots:
            plot.savefig('tests/output/loci/test_multiple_auto_reverse_gc.svg', format='svg')
        else:
            plot.show()

    def test_multiple_bokeh(self):
        loci_of_interest = [
            dict(gbk=prokka_file, gene='FAM3257_00934', title='title1'),
            dict(gbk=pgap_file, gene='FAM3257_000019', title='title2'),
            dict(gbk=pgap_file, gene='FAM3257_001020', title='title3'),
        ]

        locus_to_color_dict = {locus['gene']: '#1984ff' for locus in loci_of_interest}

        loci = Loci.generate(
            loci_of_interest,
            locus_to_color_dict=locus_to_color_dict,
            span=30000
        )

        plots = loci.plot_bokeh(viewspan=3000, auto_reverse=True)

        bokeh = column(plots)
        bokeh.sizing_mode = 'scale_width'

        if save_plots:
            output_file(filename='tests/output/loci/test_multiple_bokeh.html', )
            save(bokeh)
        else:
            show(bokeh)
