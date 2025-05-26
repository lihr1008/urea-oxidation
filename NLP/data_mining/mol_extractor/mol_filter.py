import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from cuspy import ConfigUtils
from fastode import FastLog

from utils.chem_utils import ChemUtils


RDLogger.DisableLog('rdApp.*')


class MolFilter:

    def __init__(self, config):
        self._target_mol_type = config.target_mol_type
        self._unique_mols_fp = config.unique_mols_fp
        self._filter_mols_fp = config.filter_mols_fp
        self._logger = FastLog(config.log_fp, 'INFO').logger

    def _filter(self, ser: pd.Series):
        if self._target_mol_type == 'ALLOY':
            return self._alloy_filter(ser)
        elif self._target_mol_type == 'ORG':
            return self._org_filter(ser)
        else:
            return True

    def _alloy_filter(self, ser):
        smiles = ser['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            self._logger.info(f"wrong smiles: {ser['smiles']}, name: {ser['name']}")
            return False


        for atom in mol.GetAtoms():
            if ChemUtils.is_metal(atom.GetAtomicNum()):
                return True
        return False

    def _org_filter(self, ser):
        smiles = ser['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            self._logger.info(f"wrong smiles: {ser['smiles']}, name: {ser['name']}")
            return False
        return mol.GetNumBonds() > 0

    def process(self):
        mols_df = pd.read_csv(self._unique_mols_fp, sep='\t', encoding='utf-8')
        filter_mol_df = mols_df.loc[mols_df.smiles.map(lambda x: (not pd.isna(x)) and len(x) > 0)]
        filter_mol_df = filter_mol_df.loc[filter_mol_df.apply(self._filter, axis=1)]
        filter_mol_df = filter_mol_df.sort_values(by='distance', axis=0)
        filter_mol_df.to_csv(self._filter_mols_fp, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    mf = MolFilter(ConfigUtils.load_config('../config.json').proj_config)
    mf.process()
