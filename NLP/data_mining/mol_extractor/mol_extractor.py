from collections import namedtuple
import xml
from xml.dom import minidom

from tqdm import tqdm
import pandas as pd
from rdkit import RDLogger
from fastode import FastXML, FastLog
from cuspy import ConfigUtils

from utils.formula_utils import FormulaUtils
from utils.tsv_utils import TSVUtils
from utils.graph_utils import GraphUtils


RDLogger.DisableLog("rdApp.*")
Mol = namedtuple('Mol', ['name', 'smiles', 'distance', 'keyword'])


class MolExtractor:

    @classmethod
    def extract_mols_from_xml(cls, xml_str: str, keywords: [str]) -> [Mol]:
        # doc = FastXML.parse_string(xml_str)
        graph = GraphUtils.xml_to_graph(xml_str)
        mols = []
        for chem_vertex, keyword, distance in GraphUtils.get_nearest_keyword_cner(graph, keywords):
            mols.append(Mol(name=chem_vertex.value,
                            smiles=chem_vertex.smiles,
                            distance=distance,
                            keyword=keyword))
        return mols

    def __init__(self, config):
        self._texts_tagged_fp = config.texts_tagged_fp
        self._unique_mols_fp = config.unique_mols_fp
        self._target_mol_type = config.target_mol_type
        self._logger = FastLog(config.log_fp, 'INFO').logger
        self._synonyms_df = pd.read_csv(config.keywords_synonyms_fp, sep='\t', encoding='utf-8')

    def _get_target_mols(self, sent_el):
        sent = ' '.join(FastXML.el_to_texts(sent_el))
        mols = []
        sent_graph = GraphUtils.xml_el_to_graph(sent_el)
        leaf_vertices = GraphUtils.get_leaf_vertices(sent_graph)
        words = [vertex.value for vertex in leaf_vertices]
        tags = [vertex.tag for vertex in leaf_vertices]
        for _, row in self._synonyms_df.iterrows():
            if row.word not in words or 'OSCARCM' not in tags:
                continue
            for chem_vertex, distance in GraphUtils.get_keyword_to_cner_distance(sent_graph, row.word):
                mols.append(Mol(name=chem_vertex.value,
                                smiles=chem_vertex.smiles,
                                distance=distance,
                                keyword=row.word))
        return sent, mols

    def _get_sent_and_mols(self, sent_el) -> (str, [Mol]):
        sent = ' '.join(FastXML.el_to_texts(sent_el))
        mols = []
        for mol_el in sent_el.getElementsByTagName('OSCARCM'):

            name = ' '.join(FastXML.el_to_texts(mol_el))
            smiles = mol_el.getAttribute('smiles')
            if len(smiles) == 0 and self._target_mol_type == 'ALLOY':
                smiles = FormulaUtils.formula_to_smiles(name)
            if smiles is None or len(smiles) == 0:
                smiles = None
            if smiles is None:
                continue
            mol = Mol(name=name, smiles=smiles, distance=1, keyword='*')
            mols.append(mol)
        return sent, mols

    @classmethod
    def get_first_sent_by_chem_name(cls, xml_el, chem_name: str):
        for sent_el in xml_el.getElementsByTagName('Sentence'):
            words = FastXML.el_to_texts(sent_el)
            sent = ' '.join(words)
            if chem_name in sent:
                return sent

    def _get_sents_and_mols(self, xml_el):
        graph = GraphUtils.xml_el_to_graph(xml_el)
        keywords = self._synonyms_df.word.values
        for chem_vertex, keyword, distance in GraphUtils.get_nearest_keyword_cner(graph, keywords):
            smiles = chem_vertex.smiles
            if (smiles is None or len(smiles) == 0) and self._target_mol_type == 'ALLOY':
                smiles = FormulaUtils.formula_to_smiles(chem_vertex.value)
            sent = self.get_first_sent_by_chem_name(xml_el, chem_vertex.value)
            yield sent, Mol(name=chem_vertex.value,
                            smiles=smiles,
                            distance=distance,
                            keyword=keyword)

    def process(self):
        texts_tagged_df = pd.read_csv(self._texts_tagged_fp, sep='\t', encoding='utf-8')
        mols_df = pd.DataFrame(columns=['pid', 'tid', 'name', 'smiles', 'distance', 'keyword', 'sentences'])
        pid = -1
        num_skip = 0
        num_wrong_xml = 0
        with tqdm(total=len(texts_tagged_df))as pbar:
            for _, row in texts_tagged_df.iterrows():
                pbar.update(1)
                pbar.set_postfix_str(f"skip: {num_skip}, add: {len(mols_df)}, num_wrong_xml: {num_wrong_xml}")
                try:
                    xml = FastXML.parse_string(row.xml)
                except Exception as e:
                    num_wrong_xml += 1
                    continue
                if row.pid != pid:
                    smiles_set = set()
                    pid = row.pid
                for sent, mol in self._get_sents_and_mols(xml):
                    if mol.smiles in smiles_set:
                        num_skip += 1
                        continue
                    smiles_set.add(mol.smiles)
                    mols_df = mols_df.append({'pid': row.pid, 'tid': row.tid,
                                              'name': mol.name, 'smiles': mol.smiles,
                                              'distance': mol.distance, 'keyword': mol.keyword,
                                              'sentences': sent}, ignore_index=True)
                    if len(mols_df) >= 1000:
                        TSVUtils.df_to_tsv(mols_df, self._unique_mols_fp, mode='a')
                        mols_df = pd.DataFrame(columns=['pid', 'tid', 'name', 'smiles',
                                                        'distance', 'keyword', 'sentences'])
        if len(mols_df) > 0:
            TSVUtils.df_to_tsv(mols_df, self._unique_mols_fp, mode='a')

    def _process(self):
        texts_tagged_df = pd.read_csv(self._texts_tagged_fp, sep='\t', encoding='utf-8')
        mols_df = pd.DataFrame(columns=['pid', 'tid', 'name', 'smiles', 'sentences'])
        pid = -1
        num_skip = 0
        num_wrong_xml = 0
        with tqdm(total=len(texts_tagged_df))as pbar:
            for _, row in texts_tagged_df.iterrows():
                pbar.update(1)
                pbar.set_postfix_str(f"skip: {num_skip}, add: {len(mols_df)}, num_wrong_xml: {num_wrong_xml}")
                try:
                    xml = FastXML.parse_string(row.xml)
                except Exception as e:
                    num_wrong_xml += 1
                    continue
                if row.pid != pid:
                    smiles_set = set()
                    pid = row.pid
                for sent_el in xml.getElementsByTagName('Sentence'):
                    sent, mols = self._get_sent_and_mols(sent_el)
                    for mol in mols:
                        if mol.smiles in smiles_set:
                            num_skip += 1
                            continue
                        smiles_set.add(mol.smiles)
                        mols_df = mols_df.append({'pid': row.pid, 'tid': row.tid,
                                                  'name': mol.name, 'smiles': mol.smiles,
                                                  'distance': mol.distance, 'keyword': mol.keyword,
                                                  'sentences': sent}, ignore_index=True)
                        if len(mols_df) >= 1000:
                            TSVUtils.df_to_tsv(mols_df, self._unique_mols_fp, mode='a')
                            mols_df = pd.DataFrame(columns=['pid', 'tid', 'name', 'smiles', 'sentences'])
        if len(mols_df) > 0:
            TSVUtils.df_to_tsv(mols_df, self._unique_mols_fp, mode='a')


if __name__ == "__main__":
    me = MolExtractor(ConfigUtils.load_config('../config.json').proj_config)
    me.process()
