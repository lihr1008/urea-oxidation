import math
import os.path

import matplotlib
import pandas as pd
from matplotlib import pyplot as plt
from pyecharts import options as opts
from pyecharts.charts import Graph


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


class DrawUtils:

    @classmethod
    def draw_dict(cls, count_dict):
        plt.figure(figsize=(5, 12))
        keys = list(count_dict.keys())
        values = list(count_dict.values())
        plt.barh(keys, values)
        plt.show()

    @classmethod
    def draw_hbar(cls, names: [str], counts: [int], dp: str,
                  fig_name: str = None, ylabel: str = "Metals", xlabel: str = "Count"):
        fig_name = 'hbar.pdf' if fig_name is None else f"{fig_name}.pdf"
        names.reverse()
        counts.reverse()
        num_bar = len(names)
        height_bar = max(5.0, 12 * num_bar / 50)
        plt.figure(figsize=(4, height_bar), dpi=300)
        plt.barh(range(len(names)), counts, tick_label=names)
        for n, name in enumerate(names):
            c = counts[n]
            plt.text(c, n-0.3, c)

        max_count = max(counts)
        xlim = max_count * 1.2
        plt.xlim((0, xlim))
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.tight_layout()

        fp = os.path.join(dp, fig_name)
        plt.savefig(fp, format='pdf')
        plt.show()

    @classmethod
    def draw_heat(cls, data, names, dp: str, fig_name: str = None):

        fig_name = 'heat.pdf' if fig_name is None else f"{fig_name}.pdf"
        font_size = int((12 / len(names)) * 10) + 1
        plt.figure(figsize=(5, 5), dpi=600)
        plt.imshow(data,cmap='GnBu',vmin=0, vmax=data.max()+0.03,shape=5) #颜色
        idxes = list(range(len(names)))
        reverse_idxes = idxes.copy()
        reverse_idxes.reverse()
        reverse_names = names.copy()
        reverse_names.copy()
        plt.xticks(reverse_idxes, reverse_names, fontsize=font_size, rotation=-90)
        plt.yticks(range(len(names)), names, fontsize=font_size)
        cb=plt.colorbar(fraction=0.0463)
        cb.ax.tick_params(labelsize=15) #colorbar字号
        fp = os.path.join(dp, fig_name)
        plt.tight_layout()
        plt.savefig(fp, format='pdf')
        plt.show()

    @classmethod
    def draw_graph(cls, names: [str], counts: [int], names_link_df: pd.DataFrame,
                   idxes: [int], categories: [str], dp: str, fig_name: str = None,
                   line_width_ratio: float = 2):
        fig_name = 'graph.html' if fig_name is None else f"{fig_name}.html"
        nodes = []
        for name, count, idx in zip(names, counts, idxes):
            nodes.append(opts.GraphNode(name=name,
                                        symbol_size=math.pow(count, 1/2.5),
                                        category=idx))
        links = []
        for _, row in names_link_df.iterrows():
            n1 = row.k1
            n2 = row.k2
            c = row['count']
            if n1 not in names or n2 not in names:
                continue
            links.append(opts.GraphLink(source=n1,
                                        target=n2,
                                        value=c,
                                        linestyle_opts=opts.LineStyleOpts(width=line_width_ratio*math.log10(c))))
                                        # linestyle_opts=opts.LineStyleOpts(width=c/10)))
        fp = os.path.join(dp, fig_name)
        c = (Graph(init_opts=opts.InitOpts(width="800px", height="800px"))
             .add("",
                  nodes,
                  links,
                  categories=categories,
                  # repulsion=300,
                  is_rotate_label=True,
                  is_roam=True,
                  layout='circular',
                  # layout='force',
                  linestyle_opts=opts.LineStyleOpts(color='source', curve=0.3),
                  label_opts=opts.LabelOpts(position='right'))
             .set_global_opts(title_opts=opts.TitleOpts(""),
                              legend_opts=opts.LegendOpts(is_show=False, orient="vertical", pos_left="2%", pos_top="20%"))
             .render(fp))

    @classmethod
    def draw_graph_without_cate(cls, names: [str], counts: [int], names_link_df: pd.DataFrame,
                                dp: str, fig_name: str = None):
        fig_name = 'graph.html' if fig_name is None else f"{fig_name}.html"
        nodes = []
        for name, count in zip(names, counts):
            nodes.append(opts.GraphNode(name=name,
                                        symbol_size=math.pow(count, 1/2.5)))
        links = []
        for _, row in names_link_df.iterrows():
            n1 = row.k1
            n2 = row.k2
            c = row['count']
            if n1 not in names or n2 not in names:
                continue
            links.append(opts.GraphLink(source=n1,
                                        target=n2,
                                        value=c,
                                        linestyle_opts=opts.LineStyleOpts(width=2*math.log10(c))))
                                        # linestyle_opts=opts.LineStyleOpts(width=c/10)))
        fp = os.path.join(dp, fig_name)
        c = (Graph(init_opts=opts.InitOpts(width="800px", height="800px"))
             .add("",
                  nodes,
                  links,
                  is_rotate_label=True,
                  is_roam=True,
                  layout='circular',
                  # layout='force',
                  linestyle_opts=opts.LineStyleOpts(color='source', curve=0.3),
                  label_opts=opts.LabelOpts(position='right'))
             .set_global_opts(title_opts=opts.TitleOpts(""),
                              legend_opts=opts.LegendOpts(is_show=False, orient="vertical", pos_left="2%", pos_top="20%"))
             .render(fp))


if __name__ == "__main__":
    pass
