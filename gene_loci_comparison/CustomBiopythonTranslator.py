from Bio.SeqFeature import SeqFeature
from dna_features_viewer import GraphicFeature, BiopythonTranslator


class CustomBiopythonTranslator(BiopythonTranslator):
    label_fields = None

    def __init__(self, label_fields: [str], *args, **kwargs):
        self.label_fields = label_fields
        super().__init__(*args, **kwargs)

    def translate_feature(self, feature: SeqFeature):
        """Translate a Biopython feature into a Dna Features Viewer feature."""
        properties = dict(
            label=self.compute_feature_label(feature),
            color=self.compute_feature_color(feature),
            html=self.compute_feature_html(feature),
            fontdict=self.compute_feature_fontdict(feature),
            box_linewidth=self.compute_feature_box_linewidth(feature),
            box_color=self.compute_feature_box_color(feature),
            linewidth=self.compute_feature_linewidth(feature),
            label_link_color=self.compute_feature_label_link_color(feature),
            legend_text=self.compute_feature_legend_text(feature),
        )
        if self.features_properties is not None:
            other_properties = self.features_properties
            if hasattr(other_properties, "__call__"):
                other_properties = other_properties(feature)
            properties.update(other_properties)

        location = feature.location if feature.location_operator != 'join' else feature.location.parts[0]

        return GraphicFeature(
            start=location.start,
            end=location.end,
            strand=location.strand,
            **properties
        )
