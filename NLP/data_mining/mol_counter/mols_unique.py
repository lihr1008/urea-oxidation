from tqdm import tqdm
import pandas as pd


class MolUnique:

    @classmethod
    def _unique_name_smiles_for_lid(cls, lid_mols_table: pd.DataFrame) -> pd.DataFrame:
        unique_lid_mols_table = pd.DataFrame(columns=lid_mols_table.columns)
        for _, row in lid_mols_table.iterrows():
            if pd.isna(row.smiles):
                if row['name'] in unique_lid_mols_table['name'].values:
                    continue
                else:
                    unique_lid_mols_table = unique_lid_mols_table.append(row, ignore_index=True)
            else:
                if row.smiles in unique_lid_mols_table.smiles.values:
                    continue
                else:
                    unique_lid_mols_table = unique_lid_mols_table.append(row, ignore_index=True)
        return unique_lid_mols_table

    @classmethod
    def _unique_name_smiles(cls, mols_table: pd.DataFrame):
        unique_lids = set(mols_table.lid.values)
        unique_mols_table = pd.DataFrame(columns=mols_table.columns)
        with tqdm(total=len(unique_lids))as pbar:
            for lid in unique_lids:
                pbar.update(1)
                lid_mols_table = mols_table.query(f"lid=={lid}")
                lid_mols_table = cls._unique_name_smiles_for_lid(lid_mols_table)
                unique_mols_table = unique_mols_table.append(lid_mols_table, ignore_index=True)
        return unique_mols_table

    @classmethod
    def unique_for_same_literature(cls, mols_table_fp: str, unique_mols_fp: str):
        mols_table = pd.read_csv(mols_table_fp, sep='\t', encoding='utf-8')
        unique_mols_table = cls._unique_name_smiles(mols_table)
        unique_mols_table.to_csv(unique_mols_fp, sep='\t', index=False, encoding='utf-8')
        return unique_mols_table


if __name__ == "__main__":
    df = pd.read_csv("../data/wos_nanozyme_mols/mol_table.tsv", sep='\t', encoding='utf-8')
    MolUnique.unique_for_same_literature("../data/wos_nanozyme_mols/mol_table.tsv",
                                         "../data/wos_nanozyme_mols/unique_mol_table.tsv")

