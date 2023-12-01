from pathlib import Path


def explorer_file_path(path: str, filename: str, is_dir: bool):
    path: Path = Path(path)

    if is_dir:
        path = path / filename

    return path


def int_cell_id(explorer_cell_id: str) -> int:
    """Transforms an alphabetical cell id from the Xenium Explorer to an integer ID

    E.g., int_cell_id('aaaachba-1') = 10000"""
    code = explorer_cell_id[:-2] if explorer_cell_id[-2] == "-" else explorer_cell_id
    coefs = [ord(c) - 97 for c in code][::-1]
    return sum(value * 16**i for i, value in enumerate(coefs))


def str_cell_id(cell_id: int) -> str:
    """Transforms an integer cell ID into an Xenium Explorer alphabetical cell id

    E.g., str_cell_id(10000) = 'aaaachba-1'"""
    coefs = []
    for _ in range(8):
        cell_id, coef = divmod(cell_id, 16)
        coefs.append(coef)
    return "".join([chr(97 + coef) for coef in coefs][::-1]) + "-1"
