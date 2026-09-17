import os

from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit.Chem import Fragments
from rdkit.Chem.rdchem import BondType
from rdkit.Chem import RDConfig, FragmentCatalog
import pandas as pd
from cuspy import ConfigUtils

from utils.dict_utils import DictUtils


class OrgExtractor:

    bt_to_symbol = {BondType.SINGLE: '-',
                    BondType.DOUBLE: '=',
                    BondType.TRIPLE: '#',
                    BondType.AROMATIC: ':'}
    fname = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
    fparams = FragmentCatalog.FragCatParams(1, 6, fname)

    def __init__(self, config):
        self._mol_fp = config.org_mols_fp
        self._atom_count_fp = config.atom_count_fp
        self._bond_count_fp = config.bond_count_fp
        self._frag_count_fp = config.frag_count_fp
        self._group_count_fp = config.group_count_fp
        self._atom_link_count_fp = config.atom_link_count_fp
        self._bond_link_count_fp = config.bond_link_count_fp
        self._frag_link_count_fp = config.frag_link_count_fp
        self._group_link_count_fp = config.group_link_count_fp

    @classmethod
    def _count_atom(cls, mol, atom_count):
        symbols = set([])
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol() if not atom.GetIsAromatic() else atom.GetSymbol().lower()
            if symbol == '*':
                continue
            symbols.add(symbol)
        for symbol in symbols:
            if symbol not in atom_count.keys():
                atom_count[symbol] = 1
            else:
                atom_count[symbol] += 1
        return atom_count

    @classmethod
    def _count_atom_link(cls, mol, atom_link_count, all_symbols):
        symbols = set([])
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol() if not atom.GetIsAromatic() else atom.GetSymbol().lower()
            if symbol == '*':
                continue
            symbols.add(symbol)
        symbols = list(symbols)
        for n, s1 in enumerate(symbols[:-1]):
            for s2 in symbols[n+1:]:
                i1 = all_symbols.index(s1)
                i2 = all_symbols.index(s2)
                pair = (s1, s2) if i1 < i2 else (s2, s1)
                if pair not in atom_link_count.keys():
                    atom_link_count[pair] = 1
                else:
                    atom_link_count[pair] += 1
        return atom_link_count

    # region ===== bond =====

    @classmethod
    def _get_bond_smarts(cls, bond):
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        begin_symbol = begin_atom.GetSymbol() if not begin_atom.GetIsAromatic() else begin_atom.GetSymbol().lower()
        end_symbol = end_atom.GetSymbol() if not end_atom.GetIsAromatic() else end_atom.GetSymbol().lower()
        if begin_symbol == '*' or end_symbol == '*':
            return None
        if begin_atom.GetAtomicNum() < end_atom.GetAtomicNum():
            atoms = [begin_symbol, end_symbol]
        elif begin_atom.GetAtomicNum() > end_atom.GetAtomicNum():
            atoms = [end_symbol, begin_symbol]
        elif begin_atom.GetIsAromatic():
            atoms = [end_symbol, begin_symbol]
        else:
            atoms = [begin_symbol, end_symbol]
        bond_symbol = cls.bt_to_symbol[bond.GetBondType()]
        return bond_symbol.join(atoms)

    @classmethod
    def _count_bond(cls, mol, bond_count):
        bonds = set([])
        for bond in mol.GetBonds():
            bond_smarts = cls._get_bond_smarts(bond)
            if bond_smarts is None:
                continue
            bonds.add(bond_smarts)
        for bond_smarts in bonds:
            if bond_smarts not in bond_count.keys():
                bond_count[bond_smarts] = 1
            else:
                bond_count[bond_smarts] += 1
        return bond_count

    @classmethod
    def _count_bond_link(cls,  mol, bond_link_count, all_symbols):
        bonds = set([])
        for bond in mol.GetBonds():
            bond_smarts = cls._get_bond_smarts(bond)
            if bond_smarts is None:
                continue
            bonds.add(bond_smarts)
        bonds = list(bonds)
        for n, s1 in enumerate(bonds[:-1]):
            for s2 in bonds[n+1:]:
                i1 = all_symbols.index(s1)
                i2 = all_symbols.index(s2)
                pair = (s1, s2) if i1 < i2 else (s2, s1)
                if pair not in bond_link_count.keys():
                    bond_link_count[pair] = 1
                else:
                    bond_link_count[pair] += 1
        return bond_link_count

    # endregion

    @classmethod
    def _count_frag(cls, mol, frag_count):
        frags = set([])
        for key, frag in Fragments.fns:
            if frag(mol) > 0:
                frags.add(key)
        for frag_key in frags:
            if frag_key not in frag_count.keys():
                frag_count[frag_key] = 1
            else:
                frag_count[frag_key] += 1
        return frag_count

    @classmethod
    def _count_group(cls, mol, group_count):
        group_names = set([])
        for n in range(cls.fparams.GetNumFuncGroups()):
            group = cls.fparams.GetFuncGroup(n)
            name = group.GetProp('_Name')
            if name in group_names:
                continue
            match = mol.GetSubstructMatch(group)
            if len(match) > 0:
                group_names.add(name)
                if name not in group_count.keys():
                    group_count[name] = 1
                else:
                    group_count[name] += 1
        return group_count

    @classmethod
    def _count_group_link(cls, mol, group_link_count, all_names):
        group_names = set([])
        for n in range(cls.fparams.GetNumFuncGroups()):
            group = cls.fparams.GetFuncGroup(n)
            name = group.GetProp('_Name')
            if name in group_names:
                continue
            match = mol.GetSubstructMatch(group)
            if len(match) > 0:
                group_names.add(name)
        group_names = list(group_names)
        for i, s1 in enumerate(group_names[:-1]):
            for s2 in group_names[i+1:]:
                i1 = all_names.index(s1)
                i2 = all_names.index(s2)
                pair = (s1, s2) if i1 < i2 else (s2, s1)
                if pair not in group_link_count.keys():
                    group_link_count[pair] = 1
                else:
                    group_link_count[pair] += 1
        return group_link_count

    @classmethod
    def _save_link_count(cls, link_count, fp: str):
        link_count = DictUtils.sort_count_dict(link_count, reverse=True)
        k1, k2 = zip(*link_count.keys())
        df = pd.DataFrame({'k1': k1, 'k2': k2, 'count': link_count.values()})
        df.to_csv(fp, sep='\t', encoding='utf-8', index=False)

    def process(self):
        mol_df = pd.read_csv(self._mol_fp, sep='\t', encoding='utf-8')
        atom_count = {}
        bond_count = {}
        frag_count = {}
        group_count = {}
        atom_link_count = {}
        bond_link_count = {}
        frag_link_count = {}
        group_link_count = {}
        with tqdm(total=len(mol_df))as pbar:
            for _, row in mol_df.iterrows():
                pbar.update(1)
                mol = AllChem.MolFromSmiles(row.smiles)
                if mol is None:
                    continue
                atom_count = self._count_atom(mol, atom_count)
                bond_count = self._count_bond(mol, bond_count)
                frag_count = self._count_frag(mol, frag_count)
                group_count = self._count_group(mol, group_count)

        atom_count = DictUtils.sort_count_dict(atom_count, reverse=True)
        bond_count = DictUtils.sort_count_dict(bond_count, reverse=True)
        frag_count = DictUtils.sort_count_dict(frag_count, reverse=True)
        group_count = DictUtils.sort_count_dict(group_count, reverse=True)
        atom_count_df = pd.DataFrame({'atom': atom_count.keys(), 'count': atom_count.values()})
        bond_count_df = pd.DataFrame({'bond': bond_count.keys(), 'count': bond_count.values()})
        frag_count_df = pd.DataFrame({'frag': frag_count.keys(), 'count': frag_count.values()})
        group_count_df = pd.DataFrame({'group': group_count.keys(), 'count': group_count.values()})
        atom_count_df.to_csv(self._atom_count_fp, sep='\t', encoding='utf-8', index=False)
        bond_count_df.to_csv(self._bond_count_fp, sep='\t', encoding='utf-8', index=False)
        frag_count_df.to_csv(self._frag_count_fp, sep='\t', encoding='utf-8', index=False)
        group_count_df.to_csv(self._group_count_fp, sep='\t', encoding='utf-8', index=False)

        with tqdm(total=len(mol_df))as pbar:
            for _, row in mol_df.iterrows():
                pbar.update(1)
                mol = AllChem.MolFromSmiles(row.smiles)
                if mol is None:
                    continue
                atom_link_count = self._count_atom_link(mol, atom_link_count, list(atom_count_df.atom.values))
                bond_link_count = self._count_bond_link(mol, bond_link_count, list(bond_count_df.bond.values))
                group_link_count = self._count_group_link(mol, group_link_count, list(group_count_df.group.values))
        self._save_link_count(atom_link_count, self._atom_link_count_fp)
        self._save_link_count(bond_link_count, self._bond_link_count_fp)
        self._save_link_count(group_link_count, self._group_link_count_fp)


if __name__ == "__main__":
    oe = OrgExtractor(ConfigUtils.load_config('../config.json').proj_config)
    oe.process()
