from typing import Iterator
import math
from collections import namedtuple
from enum import Enum
from xml.dom import minidom

import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
from fastode import FastXML

Vertex = namedtuple('Vertex', ['tag', 'type', 'smiles', 'value', 'idx'])


class GraphUtils:

    SIZE_VERTEX = 1

    @classmethod
    def get_leaf_vertices(cls, graph: nx.Graph) -> [Vertex]:
        leaf_vertices = []
        for node in graph.nodes:
            if node.value is not None:
                leaf_vertices.append(node)
        leaf_vertices = sorted(leaf_vertices, key=lambda x: x.idx)
        return leaf_vertices

    @classmethod
    def _get_children(cls, graph: nx.Graph, vertex: Vertex) -> [Vertex]:
        vertices = []
        for v in graph.neighbors(vertex):
            if v.idx > vertex.idx:
                vertices.append(v)
        return vertices

    @classmethod
    def _tree_pos(cls, graph) -> {Vertex, np.ndarray}:
        pos = {}
        leaf_vertices = cls.get_leaf_vertices(graph)
        for n, vertex in enumerate(leaf_vertices):
            pos[vertex] = [n*cls.SIZE_VERTEX, 0]
        while True:
            if len(pos) == graph.number_of_nodes():
                break
            for vertex in graph.nodes:
                if vertex in pos.keys():
                    continue
                nbrs = cls._get_children(graph, vertex)
                if any([nbr not in pos.keys() for nbr in nbrs]):
                    continue
                # print(len(pos))
                xs = [pos[nbr][0] for nbr in nbrs]
                ys = [pos[nbr][1] for nbr in nbrs]
                pos[vertex] = [sum(xs)/len(xs), max(ys)+cls.SIZE_VERTEX*2]
        return pos

    @classmethod
    def _draw_vertices(cls, pos, axes):
        xs, ys = list(zip(*list(pos.values())))
        axes.scatter(xs, ys, s=20, c='white', edgecolors='gray', linewidths=1, zorder=2)

    @classmethod
    def _draw_lines(cls, graph: nx.Graph, pos, axes):
        """
        1 / (x + g) + h = y
        :param graph:
        :param pos:
        :return:
        """
        for edge in graph.edges:
            v1, v2 = (edge[0], edge[1]) if edge[0].idx < edge[1].idx else (edge[1], edge[0])
            x1, y1 = pos[v1]
            x2, y2 = pos[v2]
            axes.plot([x2, x2], [y2, y1], color='gray', zorder=1)
            axes.plot([x2, x1], [y1, y1], color='gray', zorder=1)

    @classmethod
    def _get_vertex_by_value(cls, pos, value):
        for vertex in pos.keys():
            if vertex.value is not None and vertex.value.lower() == value:
                return vertex
        return None

    @classmethod
    def _get_path_distance(cls, path: [Vertex]) -> float:
        distance = len(path)
        for vertex in path:
            if vertex.tag == 'Document':
                return 10 + distance
            elif vertex.tag == 'Sentence':
                return 2 + distance
            elif 'Phrase' in vertex.tag:
                return 1 + distance
        return distance

    @classmethod
    def get_nearest_keyword_cner(cls, graph: nx.Graph, keywords: [str]):
        k_vertices = []
        for vertex in graph.nodes:
            if vertex.value is not None and vertex.value.lower() in keywords:
                k_vertices.append(vertex)
        if len(k_vertices) == 0:
            return None

        chem_vertices = []
        for vertex in graph.nodes:
            if vertex.value is not None and vertex.tag == 'OSCARCM':
                chem_vertices.append(vertex)

        for chem_vertex in chem_vertices:
            min_distance = 100
            best_keyword = None
            for k_vertex in k_vertices:
                path = nx.dijkstra_path(graph, k_vertex, chem_vertex)
                distance = cls._get_path_distance(path)
                if distance < min_distance:
                    min_distance = distance
                    best_keyword = k_vertex.value
            yield chem_vertex, best_keyword, min_distance

    @classmethod
    def get_keyword_to_cner_distance(cls, graph: nx.Graph, keyword: str) -> (Vertex, float):
        k_vertex = None
        for vertex in graph.nodes:
            if vertex.value is not None and vertex.value.lower() == keyword:
                k_vertex = vertex
                break
        if k_vertex is None:
            return None
        chem_vertices = []
        for vertex in graph.nodes:
            if vertex.value is not None and vertex.tag == 'OSCARCM':
                chem_vertices.append(vertex)

        for chem_vertex in chem_vertices:
            path = nx.dijkstra_path(graph, k_vertex, chem_vertex)
            distance = cls._get_path_distance(path)
            yield chem_vertex, distance

    @classmethod
    def _draw_path(cls, k_vertex: Vertex, n_vertex: Vertex, pos, graph: nx.Graph, axes):
        xk, yk = pos[k_vertex][0], pos[k_vertex][1]
        axes.scatter(xk, yk, s=50, marker='s', c='white', edgecolors='#0070C0', linewidths=2, zorder=3)
        path = nx.dijkstra_path(graph, k_vertex, n_vertex)
        if any([vertex.tag == 'Document' for vertex in path]):
            return
        for n, v1 in enumerate(path[:-1]):
            v2 = path[n+1]
            v1, v2 = (v1, v2) if v1.idx < v2.idx else (v2, v1)
            x1, y1 = pos[v1][0], pos[v1][1]
            x2, y2 = pos[v2][0], pos[v2][1]
            axes.plot([x2, x2], [y2, y1], color='#0070C0', zorder=1, linewidth=3)
            axes.plot([x2, x1], [y1, y1], color='#0070C0', zorder=1, linewidth=3)

    @classmethod
    def _draw_paths(cls, keywords: [str], name_vertices: [Vertex], pos, graph: nx.Graph, axes):
        for n_vertex in name_vertices:
            xn, yn = pos[n_vertex][0], pos[n_vertex][1]
            axes.scatter(xn, yn, s=50, marker='s', c='white', edgecolors='#C00000', linewidths=2, zorder=3)
        for keyword in keywords:
            k_vertex = cls._get_vertex_by_value(pos, keyword)
            if k_vertex is None:
                continue
            for n_vertex in name_vertices:
                cls._draw_path(k_vertex, n_vertex, pos, graph, axes)

    @classmethod
    def draw_tree(cls, graph: nx.Graph, keywords: [str]):
        figure, axes = plt.subplots(1, 1, figsize=(15, 5), dpi=300)
        pos = cls._tree_pos(graph)

        cls._draw_lines(graph, pos, axes)
        cls._draw_vertices(pos, axes)
        name_vertices = []
        for vertex in pos.keys():
            if vertex.value is not None and vertex.tag == 'OSCARCM':
                name_vertices.append(vertex)
        cls._draw_paths(keywords, name_vertices, pos, graph, axes)
        [axes.spines[loc_axis].set_visible(False) for loc_axis in ['top', 'right', 'bottom', 'left']]
        axes.set_xticks([])
        axes.get_yaxis().set_visible(False)
        axes.set_aspect(1)
        plt.tight_layout()
        plt.show()

    @classmethod
    def xml_to_graph(cls, xml_str: str) -> nx.Graph:
        xml = minidom.parseString(xml_str)
        return cls.xml_el_to_graph(xml)

    @classmethod
    def xml_el_to_graph(cls, xml: minidom.Element) -> nx.Graph:
        graph = nx.Graph()
        root = Vertex(tag='ROOT', type=None, smiles=None, value=None, idx=0)
        graph.add_node(root)
        graph, idx = cls._add_nodes_to_graph(graph, root, xml.firstChild)
        return graph

    @classmethod
    def _get_attrs(cls, element: minidom.Element):
        attrs = {}
        for key in ['type', 'smiles']:
            value = element.getAttribute(key)
            if len(value) != 0:
                attrs[key] = value
            else:
                attrs[key] = None
        return attrs

    @classmethod
    def _get_vertex(cls, element: minidom.Element, idx: int):
        tag = element.tagName
        attrs = cls._get_attrs(element)
        if isinstance(element.firstChild, minidom.Text):
            value = ' '.join(FastXML.el_to_texts(element))
        elif tag == 'OSCARCM':
            value = ' '.join(FastXML.el_to_texts(element))
        else:
            value = None
        return Vertex(tag=tag, type=attrs['type'], smiles=attrs['smiles'], value=value, idx=idx)

    @classmethod
    def _add_nodes_to_graph(cls, graph: nx.Graph, current_vertex: Vertex, element: minidom.Element,
                            idx: int = 0) -> (nx.Graph, int):
        idx += 1
        v = cls._get_vertex(element, idx)
        graph.add_node(v)
        graph.add_edge(current_vertex, v)
        if v.value is not None:
            return graph, idx
        else:
            for child_el in element.childNodes:
                graph, idx = cls._add_nodes_to_graph(graph, v, child_el, idx=idx)
        return graph, idx


if __name__ == "__main__":
    xml_s = '<?xml version="1.0"?><Document><Sentence><NounPhrase><CD>Three</CD><NN>novel</NN><MOLECULE><OSCARCM><OSCAR-CM>TADF</OSCAR-CM></OSCARCM><NN-CHEMENTITY>materials</NN-CHEMENTITY></MOLECULE><COMMA>,</COMMA><MOLECULE><OSCARCM><OSCAR-CM>DPAc-4PyPM</OSCAR-CM></OSCARCM></MOLECULE><COMMA>,</COMMA><MOLECULE><OSCARCM><OSCAR-CM>DPAc-6PyPM</OSCAR-CM></OSCARCM></MOLECULE><CC>and</CC><MOLECULE><OSCARCM><OSCAR-CM>DPAc-TPPM</OSCAR-CM></OSCARCM></MOLECULE></NounPhrase><VerbPhrase><VBP>have</VBP><VBN>been</VBN><VBN>designed</VBN></VerbPhrase><CC>and</CC><ActionPhrase type="Synthesize"><VerbPhrase><VB-SYNTHESIZE>synthesized</VB-SYNTHESIZE></VerbPhrase></ActionPhrase><STOP>.</STOP></Sentence><Sentence><ActionPhrase type="Add"><NounPhrase><MOLECULE><OSCARCM><OSCAR-CM>CH</OSCAR-CM></OSCARCM></MOLECULE><NN>center</NN><NN>dot</NN><NN>center</NN><NN>dot</NN><NN>center</NN><NN>dot</NN></NounPhrase><Unmatched><NN-MOLAR>N</NN-MOLAR></Unmatched><NounPhrase><NN>intramolecular</NN><MOLECULE><OSCARCM><OSCAR-CM>hydrogen-bonding</OSCAR-CM></OSCARCM></MOLECULE><NNS>interactions</NNS></NounPhrase><VerbPhrase><VBD>were</VBD><VB-ADD>introduced</VB-ADD><PrepPhrase><IN-INTO>into</IN-INTO><NounPhrase><MOLECULE><OSCARCM><OSCAR-CM>pyrimidine</OSCAR-CM><DASH>-</DASH><OSCAR-CM>pyridine</OSCAR-CM></OSCARCM></MOLECULE><NN>hybrid</NN><NNS>acceptors</NNS></NounPhrase></PrepPhrase></VerbPhrase></ActionPhrase><COMMA>,</COMMA><RB-CONJ>then</RB-CONJ><NounPhrase><DT-THE>the</DT-THE><MOLECULE><OSCARCM><OSCAR-CM>TADF</OSCAR-CM></OSCARCM></MOLECULE><NNS>properties</NNS></NounPhrase><VerbPhrase><MD>can</MD><VB>be</VB><VBN>tuned</VBN><PrepPhrase><IN-BY>by</IN-BY><NounPhrase><JJ>skillful</JJ><NN>manipulation</NN><PrepPhrase><IN-OF>of</IN-OF><NounPhrase><DT-THE>the</DT-THE><NNS>positions</NNS><PrepPhrase><IN-OF>of</IN-OF><NounPhrase><MOLECULE><OSCARCM><OSCAR-CM>CH</OSCAR-CM></OSCARCM></MOLECULE><NN>center</NN><NN>dot</NN><NN>center</NN><NN>dot</NN><NN>center</NN><NN>dot</NN></NounPhrase></PrepPhrase></NounPhrase></PrepPhrase></NounPhrase></PrepPhrase></VerbPhrase><Unmatched><NN-MOLAR>N</NN-MOLAR></Unmatched><PrepPhrase><IN-IN>in</IN-IN><NounPhrase><DT-THE>the</DT-THE><JJ-CHEM>acceptor</JJ-CHEM><NN>moiety</NN></NounPhrase></PrepPhrase><STOP>.</STOP></Sentence><Sentence><NounPhrase><DT>Both</DT><MOLECULE><OSCARCM><OSCAR-CM>DPAc-4PyPM</OSCAR-CM></OSCARCM></MOLECULE><CC>and</CC><MOLECULE><OSCARCM><OSCAR-CM>DPAc-6PyPM</OSCAR-CM></OSCARCM></MOLECULE></NounPhrase><VerbPhrase><VBP>exhibit</VBP></VerbPhrase><NounPhrase><DT>a</DT><JJR>higher</JJR><NNP>PLQY</NNP></NounPhrase><PrepPhrase><IN>than</IN><NounPhrase><MOLECULE><OSCARCM><OSCAR-CM>DPAc-TPPM</OSCAR-CM></OSCARCM></MOLECULE></NounPhrase></PrepPhrase><STOP>.</STOP></Sentence><Sentence><Unmatched><RB>Simultaneously</RB></Unmatched><Unmatched><COMMA>,</COMMA></Unmatched><NounPhrase><DT-THE>the</DT-THE><MOLECULE><OSCARCM><OSCAR-CM>k(RISC)</OSCAR-CM></OSCARCM></MOLECULE></NounPhrase><VerbPhrase><VBD>realized</VBD><PrepPhrase><IN-FOR>for</IN-FOR><NounPhrase><MOLECULE><OSCARCM><OSCAR-CM>DPAc-4PyPM</OSCAR-CM></OSCARCM></MOLECULE></NounPhrase></PrepPhrase></VerbPhrase><VerbPhrase><VBZ>is</VBZ><RB>nearly</RB></VerbPhrase><NounPhrase><CD>5</CD><CC>or</CC><CD>10</CD><NN>fold</NN></NounPhrase><JJR>higher</JJR><Unmatched><IN>than</IN></Unmatched><Unmatched><DT>that</DT></Unmatched><PrepPhrase><IN-OF>of</IN-OF><NounPhrase><PRP_>its</PRP_><NNS>analogues</NNS></NounPhrase></PrepPhrase><STOP>.</STOP></Sentence><Sentence><NounPhrase><DT-THE>The</DT-THE><JJ>optimized</JJ><JJ-CHEM>sky-blue</JJ-CHEM><NN>device</NN></NounPhrase><VerbPhrase><VB-USE>using</VB-USE></VerbPhrase><NounPhrase><MOLECULE><OSCARCM><OSCAR-CM>DPAc-4PyPM</OSCAR-CM></OSCARCM></MOLECULE></NounPhrase><PrepPhrase><IN-AS>as</IN-AS><NounPhrase><DT>a</DT><NN>dopant</NN></NounPhrase></PrepPhrase><VerbPhrase><VBD>demonstrated</VBD></VerbPhrase><NounPhrase><DT-THE>the</DT-THE><JJ>maximum</JJ><JJ>external</JJ><JJ-CHEM>quantum</JJ-CHEM><NN>efficiency</NN><COMMA>,</COMMA><JJ>current</JJ><NN>efficiency</NN><CC>and</CC><NN>power</NN><NN>efficiency</NN><PrepPhrase><IN-OF>of</IN-OF><NounPhrase><QUANTITY><PERCENT><CD>24.34</CD><NN-PERCENT>%</NN-PERCENT></PERCENT></QUANTITY><COMMA>,</COMMA><CD>53.89</CD><NN>cd</NN><MOLECULE><OSCARCM><OSCAR-CM>A(-1)</OSCAR-CM></OSCARCM></MOLECULE><CC>and</CC><CD>36.79</CD><NN>lm</NN><MOLECULE><OSCARCM><OSCAR-CM>W-1</OSCAR-CM></OSCARCM></MOLECULE></NounPhrase></PrepPhrase></NounPhrase><COMMA>,</COMMA><Unmatched><RB>respectively</RB></Unmatched><Unmatched><COMMA>,</COMMA></Unmatched><NounPhrase><WDT>which</WDT></NounPhrase><VerbPhrase><VBZ>is</VBZ><JJ>comparable</JJ></VerbPhrase><Unmatched><TO>to</TO></Unmatched><PrepPhrase><JJS>most</JJS><IN-OF>of</IN-OF><NounPhrase><DT-THE>the</DT-THE><RB>previously</RB><JJ-CHEM>reported</JJ-CHEM><JJ-CHEM>pyrimidine-based</JJ-CHEM><MOLECULE><OSCARCM><OSCAR-CM>TADF</OSCAR-CM></OSCARCM></MOLECULE><NNP>OLEDs</NNP></NounPhrase></PrepPhrase><STOP>.</STOP></Sentence></Document>'
    gr = GraphUtils.xml_to_graph(xml_s)
    GraphUtils.draw_tree(gr, ['tadf', 'device'])
