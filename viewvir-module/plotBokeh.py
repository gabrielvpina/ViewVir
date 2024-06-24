from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord
)
from bokeh.resources import CDN
from bokeh.embed import file_html


features = [
    GraphicFeature(
        start=5, end=20, strand=+1, color="#ffd700", label="Small feature"
    ),
    GraphicFeature(
        start=20,
        end=500,
        strand=+1,
        color="#ffcccc",
        label="Gene 1 with a very long name",
    ),
    GraphicFeature(
        start=400, end=700, strand=-1, color="#cffccc", label="Gene 2"
    ),
    GraphicFeature(
        start=600, end=900, strand=+1, color="#ccccff", label="Gene 3"
    ),
]


# PLOT AND EXPORT A LINEAR VIEW OF THE CONSTRUCT
record = GraphicRecord(sequence_length=1000, features=features)
plot = record.plot_with_bokeh(figure_width=15)
#ax.figure.savefig("graphic_record_defined_by_hand.png") 
with open("plot_with_bokeh.html", "w+") as f:
    f.write(file_html(plot, CDN, "Example Sequence"))