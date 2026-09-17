import math
import os
import json

from pyecharts import options
from pyecharts.charts import Graph

from utils.chem_utils import ChemUtils


class MolDraw:

    def __init__(self, data_dp: str):
        self._mols_count = self._load_json(os.path.join(data_dp, 'mols_count.json'))
        self._mols_link = self._load_json(os.path.join(data_dp, 'mols_link.json'))
        self._graph_fp = os.path.join(data_dp, 'graph.html')

    def draw_graph(self):
        nodes = []
        cols = [2, 3, 4, 5, 6, 7]
        categories = [{"name": f"{col}th column"} for col in cols]
        for symbol, count in self._mols_count.items():
            nodes.append(options.GraphNode(name=symbol,
                                           symbol_size=math.pow(count, 1/2),
                                           category=cols.index(ChemUtils.get_row(ChemUtils.get_atomic_num(symbol)))))
        links = []
        for anum_pairs, count in self._mols_link.items():
            anum1, anum2 = anum_pairs.split('_')
            anum1 = int(anum1)
            anum2 = int(anum2)
            links.append(options.GraphLink(source=ChemUtils.get_symbol(anum1),
                                           target=ChemUtils.get_symbol(anum2),
                                           value=count,
                                           linestyle_opts=options.LineStyleOpts(width=2*math.log10(count))))
        c = (Graph(init_opts=options.InitOpts(width="800px", height="800px"))
             .add("",
                  nodes,
                  links,
                  categories=categories,
                  repulsion=300,
                  is_rotate_label=True,
                  layout="circular",
                  # layout="force",
                  linestyle_opts=options.LineStyleOpts(color="source", curve=0.3),
                  label_opts=options.LabelOpts(position="right"))
             .set_global_opts(title_opts=options.TitleOpts("Graph of Metals Relations"),
                              legend_opts=options.LegendOpts(is_show=False, orient="vertical", pos_left="2%",
                                                             pos_top="20%"))
             .render(self._graph_fp))

    @classmethod
    def _load_json(cls, fp: str):
        with open(fp, 'r', encoding='utf-8')as f:
            data_count = json.load(f)
        sorted_keys = sorted(data_count.keys(), key=lambda x: data_count[x])
        sorted_data_count = {k: data_count[k] for k in sorted_keys}
        return sorted_data_count


if __name__ == "__main__":
    md = MolDraw('../data/wos_nanozyme_mols')
    md.draw_graph()
