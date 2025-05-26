from rdkit.Chem.AllChem import GetPeriodicTable


class ChemUtils:

    @classmethod
    def get_row_by_symbol(cls, symbol: str):
        atomic_num = cls.get_atomic_num(symbol)
        return cls.get_row(atomic_num)

    @classmethod
    def get_row(cls, atomic_num: int):
        if atomic_num <= 2:
            return 1
        if 3 <= atomic_num <= 10:
            return 2
        elif 11 <= atomic_num <= 18:
            return 3
        elif 19 <= atomic_num <= 36:
            return 4
        elif 37 <= atomic_num <= 54:
            return 5
        elif 55 <= atomic_num <= 86:
            return 6
        elif 87 <= atomic_num:
            return 7

    @classmethod
    def is_trans(cls, atomic_num: int) -> bool:
        return 21 < atomic_num < 30 or 39 < atomic_num < 48 or 57 < atomic_num < 80

    @classmethod
    def is_metal(cls, atomic_num: int) -> bool:
        metal_ranges = [(3, 4), (11, 13), (19, 31), (37, 50), (55, 84), (87, 108)]
        # metal_ranges = [(21, 30), (39, 48), (57, 80)]
        # metal_ranges = [(19, 22), (24, 31)]
        for metal_range in metal_ranges:
            if metal_range[0] <= atomic_num <= metal_range[1]:
                return True
        return False

    @classmethod
    def get_symbol(cls, atomic_num: int) -> str:
        symbol = GetPeriodicTable().GetElementSymbol(atomic_num)
        return symbol

    @classmethod
    def get_atomic_num(cls, symbol: str) -> int:
        try:
            atomic_num = GetPeriodicTable().GetAtomicNumber(symbol)
            return atomic_num
        except Exception as e:
            print(f"Error Symbol: {symbol}")
            print(e)
            return None


if __name__ == "__main__":
    pass
