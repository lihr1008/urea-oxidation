
import os

import pandas as pd


class WOSUtils:

    @classmethod
    def _get_docs_from_tsv(cls, fp: str):
        df = pd.read_csv(fp, sep='\t', encoding='utf-8')
        for _, row in df.iterrows():
            head = row.TI
            abt = row.AB
            doi = row.DI
            ut = row.UT
            doc = {"head": head, "abstract": abt, "doi": doi, "ut": ut}
            yield doc

    @classmethod
    def get_docs(cls, dp: str):
        for fn in os.listdir(dp):
            fp = os.path.join(dp, fn)
            for doc in cls._get_docs_from_tsv(fp):
                yield doc


if __name__ == "__main__":
    pass
