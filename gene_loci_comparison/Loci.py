import os
import matplotlib.pyplot as plt
from bokeh.models import CustomJS

from .utils import JAVASCRIPT_SYNC_SCROLL
from .Locus import Locus

DEFAULT_DESCRIPTIOIN_ORDER = [
    "locus_tag",
    "label",
    "name",
    "gene",
    "product",
    "source",
    "note"
]


class Loci:
    def __init__(self, loci: [Locus]):
        self.loci = loci

    def __str__(self) -> str:
        return f'Loci: {len(self.loci)}'

    def locus_tags(self) -> [[str]]:
        return [locus.locus_tags() for locus in self.loci]

    @staticmethod
    def generate(
            loci_of_interest: [dict],
            span=3000,
            description_order=None,
            add_start_end_feature=True,
            locus_to_color_dict=None,
            default_color='#ffffff',
            strict=False,
    ):
        if description_order is None:
            description_order = DEFAULT_DESCRIPTIOIN_ORDER

        for l in loci_of_interest:
            for key in ('gbk', 'gene', 'title'):
                assert key in l, F'Locus of interest ({l}) lacks key: {key}.'
            assert os.path.isfile(l['gbk']), F'File not found: {l}'
            assert type(l['gene']) is str, F'Gene must be string: {l}'
            assert type(l['title']) is str, F'Title must be string: {l}'

        loci: [Locus] = []
        for l in loci_of_interest:
            locus = Locus(
                gbk_file=l['gbk'],
                locus_tag=l['gene'],
                title=l['title'],
                span=span,
                description_order=description_order,
                add_start_end_feature=add_start_end_feature
            )
            locus.colorize(
                locus_to_color_dict=locus_to_color_dict,
                default_color=default_color,
                strict=strict
            )
            loci.append(locus)

        return Loci(loci)

    def plot(
            self,
            fig_single_height=2,
            fig_width=10,
            auto_reverse=True,
            *args, **kwargs
    ):
        """:returns matplotlib.pyplot module"""

        n_loci = len(self.loci)

        fig = plt.figure(figsize=(fig_width, n_loci * fig_single_height))

        lr_padding = 0.07
        tb_padding = 0
        plt.subplots_adjust(left=lr_padding, bottom=tb_padding, right=1 - lr_padding, top=1 - tb_padding,
                            wspace=None, hspace=None)

        for i, locus in enumerate(self.loci):
            locus: Locus

            ax = fig.add_subplot(n_loci, 1, i + 1)

            ax, _ = locus.plot(auto_reverse=auto_reverse, ax=ax, *args, **kwargs)

        return plt

    def plot_gc(
            self,
            fig_single_height=3,
            fig_width=10,
            auto_reverse=True,
            *args, **kwargs
    ):
        """:returns matplotlib.pyplot module"""

        n_loci = len(self.loci)

        fig, axes = plt.subplots(
            ncols=1, nrows=n_loci * 2,
            constrained_layout=True,
            gridspec_kw=dict(height_ratios=[4, 1] * n_loci),
            figsize=(fig_width, n_loci * fig_single_height)
        )

        lr_padding = 0.07
        tb_padding = 0
        plt.subplots_adjust(left=lr_padding, bottom=tb_padding, right=1 - lr_padding, top=1 - tb_padding,
                            wspace=None, hspace=None)

        for i, locus in enumerate(self.loci):
            locus: Locus

            ax1 = axes[i * 2]
            ax2 = axes[i * 2 + 1]
            ax1, ax2 = locus.plot_gc(auto_reverse=auto_reverse, ax1=ax1, ax2=ax2, *args, **kwargs)

        return plt

    def plot_bokeh(self, figure_width=12, single_figure_height='auto', viewspan=None,
                   auto_reverse=True):

        plots = []
        for current_record in self.loci:
            current_record: Locus
            p_curr = current_record.plot_bokeh(figure_width=figure_width, figure_height=single_figure_height,
                                               viewspan=viewspan, auto_reverse=auto_reverse)
            # tags: [gene_location, is_backward]
            p_curr.tags = [current_record.gene_location, current_record.is_backward]

            plots.append(p_curr)

            # Connect plot to first plot
            if len(plots) == 1:
                # assign variable first_plot, do nothing else
                p_first = p_curr
            else:
                reverse = (auto_reverse and p_first.tags[1] != p_curr.tags[1])
                for attr in ['start', 'end']:
                    p_first.x_range.js_on_change(attr, CustomJS(args=dict(
                        reverse=reverse,
                        x_range=p_curr.x_range,
                        my_center=p_first.tags[0],
                        other_center=p_curr.tags[0]
                    ), code=JAVASCRIPT_SYNC_SCROLL))

                    p_curr.x_range.js_on_change(attr, CustomJS(args=dict(
                        reverse=reverse,
                        x_range=p_first.x_range,
                        my_center=p_curr.tags[0],
                        other_center=p_first.tags[0]
                    ), code=JAVASCRIPT_SYNC_SCROLL))

        return plots
