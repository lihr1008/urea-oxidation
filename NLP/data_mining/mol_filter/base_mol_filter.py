from abc import abstractmethod

import pandas as pd


class BaseMolFilter:

    @classmethod
    def filter(cls, mols_table: pd.DataFrame) -> pd.DataFrame:
        raise NotImplementedError


if __name__ == "__main__":
    pass
