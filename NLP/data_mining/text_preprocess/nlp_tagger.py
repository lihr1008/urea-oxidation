import os.path

from tqdm import tqdm
import pandas as pd
from cuspy import ConfigUtils
from fastode import FastLog
from chem_nlpy import ChemicalTagger

from utils.tsv_utils import TSVUtils


class NLPTagger:

    def __init__(self, config):
        self._text_fp = config.texts_fp
        self._text_tagged_fp = config.texts_tagged_fp
        self._logger = FastLog(config.log_fp, 'INFO')

    def get_max_tid(self):
        if os.path.exists(self._text_tagged_fp):
            _df = pd.read_csv(self._text_tagged_fp, sep='\t', encoding='utf-8')
            max_tid = _df.tid.max()
        else:
            max_tid = -1
        return max_tid

    def process(self):
        text_df = pd.read_csv(self._text_fp, sep='\t', encoding='utf-8')
        text_tagged_df = pd.DataFrame(columns=['pid', 'tid', 'text_type', 'xml'])
        max_tid = self.get_max_tid()
        with tqdm(total=len(text_df))as pbar:
            for idx, row in text_df.iterrows():
                pbar.update(1)
                if row.tid <= max_tid:
                    continue
                xml_str = ChemicalTagger.tag_text(row.text, host="http://localhost:8088/nlpj")

                xml_str = xml_str.replace('\n', '')
                text_tagged_df = text_tagged_df.append({'pid': row.pid,
                                                        'tid': row.tid,
                                                        'text_type': row.text_type,
                                                        'xml': xml_str}, ignore_index=True)
                if len(text_tagged_df) >= 2000:
                    TSVUtils.df_to_tsv(text_tagged_df, self._text_tagged_fp, mode='a')
                    text_tagged_df = pd.DataFrame(columns=['pid', 'tid', 'text_type', 'xml'])
        if len(text_tagged_df) > 0:
            TSVUtils.df_to_tsv(text_tagged_df, self._text_tagged_fp, mode='a')
        # text_tagged_df.to_csv(self._text_tagged_fp, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    nt = NLPTagger(ConfigUtils.load_config('../config.json').proj_config)
    nt.process()
