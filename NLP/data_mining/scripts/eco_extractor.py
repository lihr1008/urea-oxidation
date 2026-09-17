import os

from tqdm import tqdm
import pandas as pd
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from fastode import FastXML


class EOCExtractor:

    @classmethod
    def _contain_x_nbr_aromatic(cls, rdmol):
        for atom in rdmol.GetAtoms():
            if atom.GetAtomicNum() not in [9, 17, 35, 53]:
                continue
            for nbr_atom in atom.GetNeighbors():
                if nbr_atom.GetIsAromatic():
                    return True
        return False

    @classmethod
    def _contain_aromatic(cls, rdmol):
        for atom in rdmol.GetAtoms():
            if atom.GetIsAromatic():
                return True
        return False

    @classmethod
    def _is_right_smiles(cls, smiles: str):
        mol = AllChem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return cls._contain_aromatic(mol)
        # return cls._contain_x_nbr_aromatic(mol)

    @classmethod
    def _is_right_text(cls, text: str):
        if 'crystal' in text or 'flexible' in text or 'elastic' in text:
            return True
        return False

    @classmethod
    def process(cls):
        dp = "../data/wos_eoc_mols"
        # unique_mols_fp = os.path.join(dp, "unique_mols.tsv")
        # unique_mol_df = pd.read_csv(unique_mols_fp, sep='\t', encoding='utf-8')
        texts_tagged_fp = os.path.join(dp, "texts_tagged.tsv")
        texts_tagged_df = pd.read_csv(texts_tagged_fp, sep='\t', encoding='utf-8')

        mols_df = pd.DataFrame(columns=['pid', 'tid', 'mol', 'sentence'])
        added_smileses = set()
        with tqdm(total=len(texts_tagged_df))as pbar:
            for _, row in texts_tagged_df.iterrows():
                pbar.update(1)
                pbar.set_postfix_str(f"num: {len(mols_df)}")
                xml = FastXML.parse_string(row.xml)
                for sent_el in xml.getElementsByTagName('Sentence'):
                    sent_words = FastXML.el_to_texts(sent_el)
                    sent = ' '.join(sent_words)
                    if not cls._is_right_text(sent):
                        continue
                    for cm_el in sent_el.getElementsByTagName('OSCARCM'):
                        names = FastXML.el_to_texts(cm_el)
                        name = ' '.join(names)
                        smiles = cm_el.getAttribute('smiles')
                        if len(smiles) == 0:
                            continue
                        if not cls._is_right_smiles(smiles):
                            continue
                        if smiles in added_smileses:
                            continue
                        added_smileses.add(smiles)
                        mols_df = mols_df.append({'pid': row.pid, 'tid': row.tid, 'name': name,
                                                  'sentence': sent, 'smiles': smiles},
                                                 ignore_index=True)
        mols_df.to_csv('../data/wos_eoc_mols/eoc_mols.tsv', sep='\t', encoding='utf-8', index=False)
        PandasTools.AddMoleculeColumnToFrame(mols_df, 'smiles', 'mol')
        mols_df = mols_df.drop(['smiles'], axis=1)
        PandasTools.SaveXlsxFromFrame(mols_df, '../data/wos_eoc_mols/eoc_mols.xlsx', molCol='mol')


if __name__ == "__main__":
    EOCExtractor.process()
