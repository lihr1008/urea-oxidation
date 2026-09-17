import os

import pandas as pd


class TSVUtils:

    @classmethod
    def df_to_tsv(cls, df: pd.DataFrame, fp: str, mode: str):
        if mode == 'w' or not os.path.exists(fp):
            df.to_csv(fp, sep='\t', encoding='utf-8', index=False)
        else:
            df.to_csv(fp, sep='\t', encoding='utf-8', index=False, header=False, mode='a')


if __name__ == "__main__":
    pass
