import os
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from dna_features_viewer import GraphicFeature


def get_locus_tag(feature: GraphicFeature) -> str:
    return feature.data['qualifiers']['locus_tag'][0]


def get_scaffold_and_geneposition(gbk_file, locus) -> (SeqRecord, int):
    assert os.path.isfile(gbk_file)
    with open(gbk_file, "r") as input_handle:
        for scf in SeqIO.parse(input_handle, "genbank"):
            for f in scf.features:
                if f.type in ["gene", "CDS"] and "locus_tag" in f.qualifiers and f.qualifiers['locus_tag'][0] == locus:
                    location = f.location if f.location_operator != 'join' else f.location.parts[0]
                    f_start, f_end = location.start, location.end
                    gene_location = f_start + (f_end - f_start) // 2
                    return (scf, gene_location)
    raise KeyError(F'Gene {locus} was not found in file {gbk_file}')


JAVASCRIPT_TAP_CALLBACK = """\
// Get label of selected datapoint:
let label
if (typeof cb_data.source.data.hover_html != "undefined") {
    // clicked on gene box
    label = cb_data.source.data.hover_html[cb_data.source.selected.indices];
} else if (typeof cb_data.source.data.text != "undefined") {
    // clicked on gene text
    label = cb_data.source.data.text[cb_data.source.selected.indices];
}
if (typeof label == "undefined" ) {
    console.log('Something was clicked on, but no label could be extracted!');
} else if (typeof geneLabelClicked == "undefined" ) {
    console.log(label, 'was clicked, but function geneLabelClicked(label) is not implemented!');
} else {
    geneLabelClicked(label, cb_data);
}
cb_data.source.change.emit();\
"""


JAVASCRIPT_SYNC_SCROLL = """\
let start, end;
if (reverse) {
start = my_center - cb_obj.start + other_center;
end   = my_center - cb_obj.end   + other_center;
} else {
start = cb_obj.start - my_center + other_center;
end   = cb_obj.end   - my_center + other_center;
}
x_range.setv({start, end});\
"""
