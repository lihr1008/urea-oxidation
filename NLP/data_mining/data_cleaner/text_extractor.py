import os

from tqdm import tqdm
import pandas as pd
from cuspy import ConfigUtils

from utils.tsv_utils import TSVUtils


class TextExtractor:

    def __init__(self, config):
        self._papers_fp = config.papers_fp
        self._texts_fp = config.texts_fp

    def process(self):
        papers_df = pd.read_csv(self._papers_fp, sep='\t', encoding='utf-8')
        texts_df = pd.DataFrame(columns=['tid', 'pid', 'text_type', 'text'])
        tid = 0
        if os.path.exists(self._texts_fp):
            os.remove(self._texts_fp)
        with tqdm(total=len(papers_df))as pbar:
            for _, row in papers_df.iterrows():
                pbar.update(1)
                pid = row['pid']
                title = row['title']
                abstract = row['abstract']
                if not pd.isna(title) and len(title) > 0:
                    texts_df = texts_df.append({'tid': tid,
                                                'pid': pid,
                                                'text_type': 'TITLE',
                                                'text': title}, ignore_index=True)
                    tid += 1
                if not pd.isna(abstract) and len(abstract) > 0:
                    texts_df = texts_df.append({'tid': tid,
                                                'pid': pid,
                                                'text_type': 'ABSTRACT',
                                                'text': abstract}, ignore_index=True)
                    tid += 1
                if 'contexts' in row.index:
                    for context in eval(row['contexts']):
                        if len(context) > 0:
                            texts_df = texts_df.append({'tid': tid,
                                                        'pid': pid,
                                                        'text_type': 'CONTEXT',
                                                        'text': context}, ignore_index=True)
                            tid += 1
                if len(texts_df) >= 2000:
                    TSVUtils.df_to_tsv(texts_df, self._texts_fp, mode='a')
                    texts_df = pd.DataFrame(columns=['tid', 'pid', 'text_type', 'text'])
        if len(texts_df) > 0:
            TSVUtils.df_to_tsv(texts_df, self._texts_fp, mode='a')

if __name__ == "__main__":
    te = TextExtractor(ConfigUtils.load_config('../config.json').proj_config)
    te.process()
