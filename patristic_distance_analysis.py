# import library
import dendropy
import pandas as pd
import numpy as np
import plotly.express as px
from plotly.offline import iplot, plot

# get phylogenetic distance from nwk files
tree = dendropy.Tree.get(path="Phylogenetic_tree_of_6species_NRCs.nwk.nwktree.nwk", schema="newick", preserve_underscores=True)
pdm = tree.phylogenetic_distance_matrix()

labels = []
distances = []
for taxon1 in tree.taxon_namespace:
    labels.append(taxon1.label)
    each_rows = []
    for taxon2 in tree.taxon_namespace:
        weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
        each_rows.append(weighted_patristic_distance)
    distances.append(each_rows)
distances = pd.DataFrame(distances)
distances.index = labels
distances.columns = labels

# get pep IDs for each plant specie
TomatoNRCs = ['Solyc08g005500.4.1','Solyc04g015220.4.1','Solyc02g027080.2.1','Solyc12g009450.3.1','Solyc12g009460.1.1','Solyc05g008070.4.1','Solyc07g049700.1.1','Solyc12g017800.3.1','Solyc11g006630.1.1','Solyc11g020090.1.1','Solyc11g006530.3.1','Solyc11g020100.3.1','Solyc11g006640.3.1','Solyc11g006520.2.1','Solyc11g069020.3.1','Solyc10g008230.3.1','Solyc10g008240.3.1','Solyc07g005770.2.1','Solyc12g005970.1.1','Solyc01g087200.3.1','Solyc07g009180.1.1','Solyc08g074250.3.1','Solyc01g016960.1.1','Solyc04g011990.3.1','Solyc04g012010.3.1','Solyc04g012000.1.1','Solyc04g011980.1.1','Solyc04g011960.2.1','Solyc10g008220.4.1','Solyc03g005660.4.1','Solyc10g047320.2.1','Solyc01g090430.3.1','Solyc04g007050.4.1','Solyc04g015210.4.1','Solyc04g007070.3.1','Solyc04g007060.4.1','Solyc04g007030.4.1','Solyc02g070410.2.1','Solyc04g007490.3.1','Solyc05g013280.4.1','Solyc05g013260.3.1','Solyc05g005130.3.1','Solyc05g007170.3.1','Solyc05g012740.2.1','Solyc05g005330.3.1','Solyc05g012890.1.1','Solyc05g012910.4.1','Solyc05g007350.3.1','Solyc05g007640.3.1','Solyc05g007610.2.1','Solyc05g007630.3.1','Solyc04g008205.1.1','Solyc04g008130.3.1','Solyc05g042090.2.1','Solyc06g008785.1.1','Solyc06g008368.1.1','Solyc05g044490.3.1','Solyc05g043420.2.1','Solyc05g054340.4.1','Solyc05g050430.2.1','Solyc05g054010.4.1','Solyc09g064690.1.1','Solyc09g098100.4.1','Solyc09g098130.3.1','Solyc12g038900.1.1','Solyc05g009750.1.1','Solyc05g009740.1.1','Solyc06g064760.2.1','Solyc06g064690.2.1','Solyc06g064750.1.1','Solyc06g064680.1.1','Solyc06g064720.1.1']
NibenNRCs = [each_label for each_label in labels if each_label[:2] == "Nb"]
SweetNRCs = [each_label for each_label in labels if each_label[:2] == "it"]
CarrotNRCs = [each_label for each_label in labels if each_label[:2] == "DC"]
CoffeeNRCs = [each_label for each_label in labels if each_label[:2] == "Cc"]
MonkeyNRCs = [each_label for each_label in labels if each_label[:3] == "Mig"]

# prepare data for plot
data = pd.DataFrame()
hovertext = []
for i, specie in enumerate([NibenNRCs, SweetNRCs, CarrotNRCs, CoffeeNRCs, MonkeyNRCs]):
    data = pd.concat([data, distances.loc[TomatoNRCs, specie].min(axis=1)])
    Tomato_vs_specie = distances.loc[TomatoNRCs, specie]
    for i in range(Tomato_vs_specie.shape[0]):
        hovertext.append(Tomato_vs_specie.columns[Tomato_vs_specie.iloc[i, :] == Tomato_vs_specie.iloc[i, :].min()].values[0])

data["x"] = np.tile(range(len(TomatoNRCs)), 5)
data["specie"] = np.repeat(["N.bentham", "SweetPotato", "Carrot", "Coffee", "Monkey"], len(TomatoNRCs))
data["hovertext"] = hovertext
data.columns = ["y", "x", "specie", "gene"]
data.head()

# plot by plotly
fig = px.scatter(data, x="x", y="y", color="specie", hover_data=['gene'],
                 color_discrete_map={
                     "N.bentham":"#84b1aa",
                     "SweetPotato":"#BA77E8",
                     "Carrot":"#ce7e00",
                     "Coffee":"#1982EA",
                     "Monkey":"#3cb371"
                 })
fig.update_xaxes(
    ticktext=TomatoNRCs,
    tickvals=list(range(len(TomatoNRCs))),
    tickangle=-90,
    tickfont=dict(family='Arial', size=7),
)
fig.update_traces(marker_size=9)

def update_symbol(trace):
    if trace.name == "N.bentham":
        trace.update(marker_symbol="circle")
    elif trace.name == "SweetPotato":
        trace.update(marker_symbol="star-square")
    elif trace.name == "Carrot":
        trace.update(marker_symbol="square")
    elif trace.name == "Coffee":
        trace.update(marker_symbol="diamond")
    elif trace.name == "Monkey":
        trace.update(marker_symbol="cross")

fig.for_each_trace(
    update_symbol
)

fig.show()
plot(fig, filename='Distance_from_Tomato_NRCs.html')