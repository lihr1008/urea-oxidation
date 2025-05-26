from rdkit.Chem import AllChem

from mol_filter.base_mol_filter import BaseMolFilter, pd
from utils.chem_utils import ChemUtils


class NanoMolFilter(BaseMolFilter):

    @classmethod
    def _is_right_smiles(cls, smiles: str) -> bool:
        if not isinstance(smiles, str):
            return False
        mol = AllChem.MolFromSmiles(smiles)
        num_metal = 0
        for atom in mol.GetAtoms():
            if ChemUtils.is_metal(atom.GetAtomicNum()):
                num_metal += 1
        return num_metal > 1

    @classmethod
    def _is_right_name(cls, name: str) -> bool:
        if not isinstance(name, str):
            return False
        return True
        # return 'Fe' in name and '(CN)' in name

    @classmethod
    def filter(cls, mols_table: pd.DataFrame):
        # filter_mol_table = mols_table.loc[mols_table.smiles.map(cls._is_right_smiles)]
        filter_mol_table = mols_table.loc[mols_table['name'].map(cls._is_right_name)]
        return filter_mol_table


if __name__ == "__main__":
    pass
