import os.path

from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit import RDLogger
import pandas as pd
from cuspy import ConfigUtils

from utils.chem_utils import ChemUtils


RDLogger.DisableLog("rdApp.*")


class AlloyExtractor:

    def __init__(self, config):
        self._unique_mols_df = pd.read_csv(config.unique_mols_fp, sep='\t', encoding='utf-8')
        self._elements_count_fp = config.elements_count_fp
        self._elements_link_count_fp = config.elements_link_count_fp

        self._elements_count = {}
        self._elements_link_count = {}

    @classmethod
    def _get_elements(cls, smiles: str):
        if pd.isna(smiles):
            return []
        mol = AllChem.MolFromSmiles(smiles)
        if mol is None:
            return []
        elements = set()
        for atom in mol.GetAtoms():
            elements.add(atom.GetSymbol())
        return list(elements)

    def _count_elements(self, elements: [str]):
        for el in elements:
            atomic_num = ChemUtils.get_atomic_num(el)
            if not ChemUtils.is_metal(atomic_num):
                continue
            if el not in self._elements_count.keys():
                self._elements_count[el] = 1
            else:
                self._elements_count[el] += 1

    def _count_elements_link(self, elements: [str]):
        for n, el1 in enumerate(elements[:-1]):
            for el2 in elements[n+1:]:
                a1 = ChemUtils.get_atomic_num(el1)
                a2 = ChemUtils.get_atomic_num(el2)
                if not ChemUtils.is_metal(a1) or not ChemUtils.is_metal(a2):
                    continue
                pair = (el1, el2) if a1 > a2 else (el2, el1)
                if pair not in self._elements_link_count.keys():
                    self._elements_link_count[pair] = 1
                else:
                    self._elements_link_count[pair] += 1

    @classmethod
    def _sort_dict(cls, dic):
        sort_keys = sorted(dic.keys(), key=lambda k: dic[k], reverse=True)
        sort_dic = {}
        for k in sort_keys:
            sort_dic[k] = dic[k]
        return sort_dic

    def _save_count(self):
        self._elements_count = self._sort_dict(self._elements_count)
        self._elements_link_count = self._sort_dict(self._elements_link_count)
        els_count_dict = {'element': self._elements_count.keys(),
                          'count': self._elements_count.values()}
        els_count_df = pd.DataFrame(els_count_dict)
        els_count_df.to_csv(self._elements_count_fp, sep='\t', encoding='utf-8', index=False)

        el1s, el2s = zip(*self._elements_link_count.keys())
        els_link_count_dict = {'k1': el1s, 'k2': el2s,
                               'count': self._elements_link_count.values()}
        els_link_count_df = pd.DataFrame(els_link_count_dict)
        els_link_count_df.to_csv(self._elements_link_count_fp, sep='\t', encoding='utf-8', index=False)

    def process(self):
        if not os.path.exists(self._elements_count_fp) or not os.path.exists(self._elements_link_count_fp):
            with tqdm(total=len(self._unique_mols_df))as pbar:
                pbar.set_description("extract & count")
                for smiles in self._unique_mols_df.smiles:
                    pbar.update(1)
                    elements = self._get_elements(smiles)
                    self._count_elements(elements)
                    self._count_elements_link(elements)
            self._save_count()


if __name__ == "__main__":
    ae = AlloyExtractor(ConfigUtils.load_config('../config.json').proj_config)
    ae.process()
