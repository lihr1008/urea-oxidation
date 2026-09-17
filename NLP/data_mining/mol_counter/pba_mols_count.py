

import pandas as pd
from rdkit.Chem import AllChem

from mol_counter.mols_unique import MolUnique
from mol_filter.pba_mol_filter import PBAMolFilter
from utils.chem_utils import ChemUtils
from utils.dict_utils import DictUtils
from utils.draw_utils import DrawUtils
from utils._formula_utils import parse_formula


class PBAMolsCount:

    @classmethod
    def filter(cls, mols_fp: str, filter_mols_fp: str):
        mols_table = pd.read_csv(mols_fp, sep='\t', encoding='utf-8')
        mols_table = MolUnique.unique_for_same_literature(mols_table)
        mols_table = PBAMolFilter.filter(mols_table)
        mols_table.to_csv(filter_mols_fp, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def get_other_element_from_smiles(cls, smiles: str):
        mol = AllChem.MolFromSmiles(smiles)
        metals = set()
        meet_fe = False
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if not ChemUtils.is_metal(atomic_num):
                continue
            if atomic_num == 26 and not meet_fe:
                meet_fe = True
                continue
            metals.add(atomic_num)
        return list(metals)

    @classmethod
    def get_other_element_from_name(cls, name: str):
        try:
            symbol_count = parse_formula(name)
        except Exception as e:
            return []
        metals = set()
        for key, value in symbol_count.items():
            atomic_num = ChemUtils.get_atomic_num(key)
            if atomic_num is None:
                continue
            if not ChemUtils.is_metal(atomic_num):
                continue
            if atomic_num == 26 and value == 1:
                continue
            metals.add(atomic_num)
        return list(metals)

    @classmethod
    def count(cls, filter_mols_fp: str):
        mols_table = pd.read_csv(filter_mols_fp, sep='\t', encoding='utf-8')
        all_metals = []
        for _, row in mols_table.iterrows():
            # metals = cls.get_other_element_from_smiles(row.smiles)
            metals = cls.get_other_element_from_name(row['name'])
            all_metals.extend(metals)
        metal_count = {}
        for metal in set(all_metals):
            metal_count[ChemUtils.get_symbol(metal)] = all_metals.count(metal)
        sorted_metal_count = DictUtils.sort_count_dict(metal_count)
        DrawUtils.draw_dict(sorted_metal_count)


if __name__ == "__main__":
    mfp = '../data/wos_pba_mols/mol_table.tsv'
    fmfp = '../data/wos_pba_mols/filter_mol_table.tsv'
    PBAMolsCount.filter(mfp, fmfp)
    PBAMolsCount.count(fmfp)
