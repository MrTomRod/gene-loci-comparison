from io import BytesIO
import numpy as np
import matplotlib.pyplot as plt
from Bio.SeqFeature import SeqFeature
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.models import Range1d, TapTool, CustomJS
from bokeh.plotting._tools import process_tools_arg

from .utils import get_locus_tag, get_scaffold_and_geneposition, JAVASCRIPT_TAP_CALLBACK
from .CustomBiopythonTranslator import CustomBiopythonTranslator


class Locus:
    default_description_order = [
        "locus_tag",
        "label",
        "name",
        "gene",
        "product",
        "source",
        "note"
    ]

    def __init__(self, gbk_file, locus_tag, title=None, span=3000, add_start_end_feature=True,
                 description_order: [str] = default_description_order):
        self.title = title
        self.gbk_file = gbk_file
        self.locus_tag = locus_tag
        self.span = span

        self.scaffold, self.gene_location = get_scaffold_and_geneposition(gbk_file, locus_tag)

        self.scaffold_id = self.scaffold.id

        unique_features: {(int, int, str)} = set()

        def add_unique(f: SeqFeature) -> bool:
            if 'locus_tag' not in f.qualifiers:
                return False
            feature = (f.location.nofuzzy_start, f.location.nofuzzy_end, f.qualifiers['locus_tag'][0])
            if feature in unique_features:
                return False
            else:
                unique_features.add(feature)
                return True

        self.graphic_record: GraphicRecord = CustomBiopythonTranslator(
            label_fields=description_order,
            features_filters=[add_unique, lambda f: f.type != 'source'],
            features_properties=lambda f: dict(qualifiers=f.qualifiers)
        ).translate_record(self.scaffold)

        self.scaffold_start = 0
        self.scaffold_end = self.graphic_record.sequence_length

        self.crop_window = self._crop_coordinates()

        self.graphic_record = self.graphic_record.crop(self.crop_window)

        if add_start_end_feature:
            self._add_start_and_end_feature()

    def __str__(self) -> str:
        return f'Locus: {self.title} ({self.locus_tag})'

    @property
    def is_backward(self) -> bool:
        for feature in self.graphic_record.features:
            if get_locus_tag(feature) == self.locus_tag:
                assert feature.strand in [1, -1]
                return feature.strand == -1
        raise KeyError("Error: Could not find locus_tag in graphic features")

    def locus_tags(self) -> [str]:
        return [get_locus_tag(f) for f in self.graphic_record.features]

    def rename_labels(self, locus_to_new_name_dict, strict=False, remove_unspecified=False) -> GraphicRecord:
        for f in self.graphic_record.features:
            locus_tag = get_locus_tag(f)
            if locus_tag in locus_to_new_name_dict:
                f.label = locus_to_new_name_dict[locus_tag]
            else:
                if strict:
                    raise KeyError(F'Label {f.label} not found in locus_to_color_dict!')
                if remove_unspecified:
                    f.label = None
        return self.graphic_record

    def colorize(self, locus_to_color_dict, strict=False,
                 default_color='#ffffff'):
        for f in self.graphic_record.features:
            locus_tag = get_locus_tag(f)
            if locus_tag in locus_to_color_dict:
                f.color = locus_to_color_dict[locus_tag]
            else:
                if strict:
                    raise KeyError(F'Locus tag {locus_tag} not found in locus_to_color_dict!')
                f.color = default_color

    def plot(self, auto_reverse=True, add_title=True, title_kwargs=dict(x=0.5, y=0.7, horizontalalignment='center', fontsize=20), *args, **kwargs):
        ax, _ = self.graphic_record.plot(*args, **kwargs)

        # set plot span
        if auto_reverse and self.is_backward:
            ax.set_xlim(self.gene_location + self.span, self.gene_location - self.span)
        else:
            ax.set_xlim(self.gene_location - self.span, self.gene_location + self.span)

        # add title
        if add_title and self.title is not None:
            # ax.set_title(graphic_record.title)
            plt.text(s=self.title, transform=ax.transAxes, **title_kwargs)

        return ax, _

    def plot_to_string(self, tight_layout=True, *args, **kwargs):
        ax, _ = self.plot(*args, **kwargs)
        if tight_layout: ax.figure.tight_layout()
        f = BytesIO()
        plt.savefig(f, format="svg")
        return f.getvalue().decode('utf-8')

    def plot_gc(self, ax1=None, ax2=None, window_bp=100, *args, **kwargs):
        if ax1 is None or ax2 is None:
            fig, (ax1, ax2) = plt.subplots(
                2, 1, figsize=(12, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
            )

        # PLOT THE RECORD MAP
        self.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4, *args, **kwargs)

        # PLOT THE LOCAL GC CONTENT (we use 50bp windows)
        calc_gc = lambda s: 100.0 * len([c for c in s if c in "GC"]) / window_bp
        xx = np.arange(self.crop_window[0] - window_bp, self.crop_window[1] + window_bp)
        yy = [calc_gc(self.scaffold.seq[x: x + window_bp]) for x in xx]
        ax2.fill_between(xx + window_bp / 2, yy, alpha=0.3)
        ax2.set_ylim(bottom=0, top=100)
        ax2.set_ylabel("GC(%)")

        # ensure ax2 has same xlim as ax1
        ax2.set_xlim(ax1.get_xlim())

        return ax1, ax2

    def plot_bokeh(self, figure_width=12, figure_height='auto', viewspan=None, auto_reverse=True,
                   x_range=None):
        if not viewspan:
            viewspan = self.crop_window
        else:
            if auto_reverse and self.is_backward:
                viewspan = (self.gene_location + viewspan, self.gene_location - viewspan)
            else:
                viewspan = (self.gene_location - viewspan, self.gene_location + viewspan)

        bokeh = self.graphic_record.plot_with_bokeh(figure_width=figure_width, figure_height=figure_height)

        # autoscale plot
        bokeh.sizing_mode = 'scale_width'

        # overwrite plot tools. reason: remove hover tool
        tap = TapTool()

        tap.callback = CustomJS(code=JAVASCRIPT_TAP_CALLBACK)

        tool_objs, tool_map = process_tools_arg(bokeh, [tap, "xpan,xwheel_zoom,reset"])
        bokeh.tools = tool_objs

        if x_range:
            bokeh.x_range = x_range
        else:
            bokeh.x_range = Range1d(*viewspan)

        return bokeh

    def _crop_coordinates(self):
        assert self.scaffold_start <= self.gene_location <= self.scaffold_end

        crop_start = self.scaffold_start if self.gene_location - self.span <= self.scaffold_start else self.gene_location - self.span
        crop_end = self.scaffold_end if self.gene_location + self.span >= self.scaffold_end else self.gene_location + self.span

        return (crop_start, crop_end)

    def _add_start_and_end_feature(self, feature_span=30, color='#000000'):
        assert 0 <= self.scaffold_start and 0 < feature_span < self.scaffold_end
        add_start = add_end = False
        if self.gene_location - self.span <= self.scaffold_start:
            add_start = True
        if self.gene_location + self.span >= self.scaffold_end:
            add_end = True

        if not add_start and not add_end:
            return

        to_keep = []
        if add_start:
            start_feature = GraphicFeature(label=F'Start of contig\n{self.scaffold_id}',
                                           start=self.scaffold_start, end=self.scaffold_start + feature_span,
                                           strand=+1, color=color)
            start_feature.data['qualifiers'] = {'locus_tag': ['Start of contig']}
            to_keep.append(start_feature)

        for f in self.graphic_record.features:
            to_keep.append(f)

        if add_end:
            end_feature = GraphicFeature(label=F'End of contig\n{self.scaffold_id}',
                                         start=self.scaffold_end - feature_span, end=self.scaffold_end,
                                         strand=-1, color=color)
            end_feature.data['qualifiers'] = {'locus_tag': ['End of contig']}
            to_keep.append(end_feature)

        self.graphic_record = GraphicRecord(
            sequence=self.graphic_record.sequence,
            sequence_length=self.graphic_record.sequence_length,
            features=to_keep,
            feature_level_height=self.graphic_record.feature_level_height,
            first_index=self.graphic_record.first_index,
        )
