from typing import Iterator
from collections import namedtuple
from enum import Enum
import os
import logging

import pandas as pd
from cuspy import ConfigUtils

logging.basicConfig(level=logging.INFO)


class TextType(Enum):

    TITLE = 0
    ABSTRACT = 1


Paper = namedtuple('Paper', ['source_fp', 'source_id', 'doi', 'title', 'abstract'])


class TextTagged:

    def __init__(self, config):
        self._papers_dp = config.eoc_config.papers_dp
        self._papers_tagged_fp = config.eoc_config.papers_tagged_fp

    def _papers_iter(self) -> Iterator[Paper]:
        for fn in os.listdir(self._papers_dp):
            fp = os.path.join(self._papers_dp, fn)
            logging.info(f"open file: {fp}")
            df = pd.read_csv(fp, sep='\t', encoding='utf-8')
            for idx, row in df.iterrows():
                paper = Paper(source_fp=fp,
                              source_id=idx,
                              doi=row.DI,
                              title=row.TI,
                              abstract=row.AB)
                yield paper

    def process(self):
        pass


if __name__ == "__main__":
    tt = TextTagged(ConfigUtils.load_config('../config.json'))
    tt.process()
