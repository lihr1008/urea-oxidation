from typing import Iterator
import os
import logging

from tqdm import tqdm
import pandas as pd
from cuspy import ConfigUtils

from utils.tsv_utils import TSVUtils


class PaperExtractor:

    def __init__(self, config):
        self._origin_papers_dp = config.origin_papers_dp
        self._papers_fp = config.papers_fp

    def _iter_origin_paper_df(self) -> Iterator[pd.DataFrame]:
        for fn in os.listdir(self._origin_papers_dp):
            fp = os.path.join(self._origin_papers_dp, fn)
            yield pd.read_csv(fp, sep='\t', encoding='utf-8')

    @classmethod
    def wos_to_df(cls, wos_fp: str):
        paper_df = pd.DataFrame(columns=['pid', 'source_fp', 'source_id', 'doi', 'title', 'abstract'])
        try:
            wos_df = pd.read_csv(wos_fp, sep='\t', encoding='utf-8')
        except Exception as e:
            return None
        for idx, row in wos_df.iterrows():
            paper_df = paper_df.append({'pid': idx,
                                        'source_fp': wos_fp,
                                        'source_id': idx,
                                        'doi': row.DI,
                                        'title': row.TI,
                                        'abstract': row.AB}, ignore_index=True)
        return paper_df

    def process(self):
        paper_df = pd.DataFrame(columns=['pid', 'source_fp', 'source_id', 'doi', 'title', 'abstract'])
        pid = 0 #计数
        tot = len(list(os.listdir(self._origin_papers_dp)))
        if os.path.exists(self._papers_fp):
            os.remove(self._papers_fp)
        with tqdm(total=tot)as pbar:
            for source_fn in os.listdir(self._origin_papers_dp):
                pbar.update(1)
                source_fp = os.path.join(self._origin_papers_dp, source_fn)
                try:
                    source_df = pd.read_csv(source_fp, sep='\t', encoding='utf-8')
                except Exception as e:
                    continue
                for idx, row in source_df.iterrows():
                    paper_df = paper_df.append({'pid': pid,
                                                'source_fp': source_fp,
                                                'source_idx': idx,
                                                'doi': row.DI,
                                                'title': row.TI,
                                                'abstract': row.AB}, ignore_index=True)
                    pid += 1
                    if len(paper_df) >= 1000: #如果文献数大于1000
                        TSVUtils.df_to_tsv(paper_df, self._papers_fp, mode='a')
                        paper_df = pd.DataFrame(columns=['pid', 'source_fp', 'source_id', 'doi', 'title', 'abstract'])
        if len(paper_df) > 0:
            TSVUtils.df_to_tsv(paper_df, self._papers_fp, mode='a')
        # paper_df.to_csv(self._papers_fp, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    pe = PaperExtractor(ConfigUtils.load_config('../config.json').proj_config)
    pe.process()
