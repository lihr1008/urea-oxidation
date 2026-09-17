import os

from chem_nlpy import ChemicalTagger
import pandas as pd
from tqdm import tqdm
from table_utils import TableWriter

from data_utils.wos_utils import WOSUtils


class WOSMolExtractor:

    @classmethod
    def extract_from_paragraph(cls, para: str):
        try:
            mols = ChemicalTagger.text_to_mols(para)
        except Exception as e:
            print(e)
            return []
        return mols

    @classmethod
    def _load_extracted_lids(cls, lit_fp: str):
        if os.path.exists(lit_fp):
            df = pd.read_csv(lit_fp, sep='\t', encoding='utf-8')
            lids = set(df.lid)
        else:
            lids = set([])
        return lids

    @classmethod
    def extract(cls, wos_dp: str, result_dp: str):
        lit_table = pd.DataFrame(columns=['lid', 'doi', 'ut', 'head', 'abstract'])
        mol_table = pd.DataFrame(columns=['lid', 'mid', 'name', 'smiles'])
        lit_fp = os.path.join(result_dp, 'lit_table.tsv')
        mol_fp = os.path.join(result_dp, 'mol_table.tsv')
        mid = 0
        # parsed_lids = cls._load_extracted_lids(lit_fp)
        with tqdm(total=100000)as pbar:
            for i, doc in enumerate(WOSUtils.get_docs(wos_dp)):
                if pd.isna(doc['abstract']):
                    continue
                pbar.update(1)
                if i in [21544, 21545, 21546]:
                    continue

                doc['lid'] = i
                lit_table = lit_table.append(doc, ignore_index=True)
                doc_mols = cls.extract_from_paragraph(doc['abstract'])
                for mol in doc_mols:
                    mol['mid'] = mid
                    mol['lid'] = i
                    mid += 1
                    mol_table = mol_table.append(mol, ignore_index=True)
                if i % 1000 == 0:
                    TableWriter.write_tsv_by_df(lit_table, lit_fp)
                    TableWriter.write_tsv_by_df(mol_table, mol_fp)
                    lit_table = pd.DataFrame(columns=['lid', 'doi', 'ut', 'head', 'abstract'])
                    mol_table = pd.DataFrame(columns=['lid', 'mid', 'name', 'smiles'])
        TableWriter.write_tsv_by_df(lit_table, lit_fp)
        TableWriter.write_tsv_by_df(mol_table, mol_fp)
        # lit_table.to_csv(lit_fp, sep='\t', encoding='utf-8', index=False)
        # mol_table.to_csv(mol_fp, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    dp = "C:\\Users\\zhang\\OneDrive\\Projects\\CHEM-NLP\\Code\\chem_extractor\\data\\wos_eoc_mols\\wos"
    WOSMolExtractor.extract(dp, "C:\\Users\\zhang\\OneDrive\\Projects\\CHEM-NLP\\Code\\chem_extractor\\data\\wos_eoc_mols")
    # dp = "/home/zhangbc/Data/chem_extractor/wos_nanozyme"
    # WOSMolExtractor.extract(dp, '/home/zhangbc/Data/chem_extractor/wos_nanozyme_mols')
    # dp = "C:\\Users\\zhang\\OneDrive\\Projects\\Spider\\Data\\wos_nanozyme"
    # WOSMolExtractor.extract(dp, '../data/wos_nanozyme_mols')
    # dp = "C:\\Users\\zhang\\OneDrive\\Projects\\Spider\\Data\\wos_pba"
    # WOSMolExtractor.extract(dp, '../data/wos_pba_mols')

    # dp = "C:\\Users\\zhang\\OneDrive\\Projects\\Spider\\Data\\wos_urea_oxidation"
    # WOSMolExtractor.extract(dp, '../data/wos_urea_oxidation_mols')
