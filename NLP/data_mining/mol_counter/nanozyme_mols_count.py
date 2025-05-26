
import json

import pandas as pd
from rdkit.Chem import AllChem
from rdkit import rdBase

from mol_counter.mols_unique import MolUnique
from mol_filter.nanozyme_mol_filter import NanoMolFilter as MolFilter
from utils.chem_utils import ChemUtils
from utils.dict_utils import DictUtils
from utils.draw_utils import DrawUtils
from utils._formula_utils import parse_formula

rdBase.DisableLog('rdApp.*')


class NanoMolsCount:

    @classmethod
    def filter(cls, unique_mols_fp: str, filter_mols_fp: str):
        # mols_table = pd.read_csv(mols_fp, sep='\t', encoding='utf-8')
        mols_table = pd.read_csv(unique_mols_fp, sep='\t', encoding='utf-8')
        mols_table = MolFilter.filter(mols_table)
        mols_table.to_csv(filter_mols_fp, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def get_metal_element_from_smiles(cls, smiles: str):
        if pd.isna(smiles):
            return []
        mol = AllChem.MolFromSmiles(smiles)
        if mol is None:
            return []
        metals = set()
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if not ChemUtils.is_metal(atomic_num):
                continue
            metals.add(atomic_num)
        return list(metals)

    @classmethod
    def get_metal_element_from_name(cls, name: str):
        try:
            symbol_count = parse_formula(name)
        except Exception as e:
            return []
        metals = set()
        all_is_upper = name.isupper()
        for key, value in symbol_count.items():
            atomic_num = ChemUtils.get_atomic_num(key)
            if atomic_num is None:
                if all_is_upper:
                    return []
                continue
            if not ChemUtils.is_metal(atomic_num):
                continue
            metals.add(atomic_num)
        return list(metals)

    @classmethod
    def count(cls, filter_mols_fp: str, mols_count_fp: str, mols_link_fp: str):
        mols_table = pd.read_csv(filter_mols_fp, sep='\t', encoding='utf-8')
        all_metals = []
        metal_link = {}
        for _, row in mols_table.iterrows():
            metals = cls.get_metal_element_from_smiles(row.smiles)
            if len(metals) == 0:
                metals = cls.get_metal_element_from_name(row['name'])

            # add link
            if len(metals) >= 2:
                for i, m1 in enumerate(metals[:-1]):
                    for m2 in metals[i+1:]:
                        link_symbol = f"{m1}_{m2}" if m1 < m2 else f"{m2}_{m1}"
                        if link_symbol not in metal_link.keys():
                            metal_link[link_symbol] = 0
                        metal_link[link_symbol] += 1
            all_metals.extend(metals)
        metal_count = {}
        for metal in set(all_metals):
            metal_count[ChemUtils.get_symbol(metal)] = all_metals.count(metal)
        sorted_metal_count = DictUtils.sort_count_dict(metal_count)
        with open(mols_count_fp, 'w', encoding='utf-8')as f:
            str_metal_count = {str(k): v for k, v in metal_count.items()}
            json.dump(str_metal_count, f)
        with open(mols_link_fp, 'w', encoding='utf-8')as f:
            json.dump(metal_link, f)
        # DrawUtils.draw_dict(sorted_metal_count)


if __name__ == "__main__":
    mfp = '../data/wos_nanozyme_mols/mol_table.tsv'
    umfp = '../data/wos_nanozyme_mols/unique_mol_table.tsv'
    fmfp = '../data/wos_nanozyme_mols/filter_mol_table.tsv'
    mcfp = '../data/wos_nanozyme_mols/mols_count.json'
    mlfp = '../data/wos_nanozyme_mols/mols_link.json'
    # NanoMolsCount.filter(umfp, fmfp)
    NanoMolsCount.count(fmfp, mcfp, mlfp)
