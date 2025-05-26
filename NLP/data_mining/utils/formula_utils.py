import re

from rdkit.Chem.AllChem import GetPeriodicTable


class FormulaUtils:

    @classmethod
    def _split_by_non_letter(cls, formula: str) -> [str]:
        return re.split(f'[^a-zA-Z]', formula)

    @classmethod
    def _split_by_upper(cls, tokens: [str]) -> [str]:
        new_tokens = []
        for token in tokens:
            new_tokens.extend(re.split('(?=[A-Z])', token))
        return new_tokens

    @classmethod
    def _is_elements(cls, symbol: str) -> bool:
        try:
            atomic_num = GetPeriodicTable().GetAtomicNumber(symbol)
            return True
        except Exception as e:
            return False

    @classmethod
    def _extract_elements(cls, tokens: [str]) -> [str]:
        symbols = []
        for token in tokens:
            if len(token) == 0 or token.islower():
                continue
            if len(token) > 3:
                return []
            is_element = False
            for l in range(min(len(token), 3), 0, -1):
                symbol = token[:l]
                if cls._is_elements(symbol):
                    symbols.append(symbol)
                    is_element = True
                    break
            if not is_element:
                return []
        return symbols

    @classmethod
    def formula_to_elements(cls, formula: str):
        tokens = cls._split_by_non_letter(formula)
        tokens = cls._split_by_upper(tokens)
        elements = cls._extract_elements(tokens)
        return elements

    @classmethod
    def formula_to_smiles(cls, formula: str):
        elements = cls.formula_to_elements(formula)
        if len(elements) == 0:
            return None
        return '.'.join([f"[{symbol}]" for symbol in elements])


if __name__ == "__main__":
    print(FormulaUtils.formula_to_smiles("ADVANTAGE"))
