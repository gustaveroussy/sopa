from pathlib import Path


def explorer_file_path(path: str, filename: str, is_dir: bool):
    path: Path = Path(path)

    if is_dir:
        path = path / filename

    return path


def int_cell_id(explorer_cell_id: str) -> int:
    """Transforms an alphabetical cell id from the Xenium Explorer to an integer ID"""
    coefs = [ord(c) - 97 for c in explorer_cell_id][::-1]
    return sum(value * 26**i for i, value in enumerate(coefs))


def str_cell_id(cell_id: int) -> str:
    """Transforms an integer cell ID into an Xenium Explorer alphabetical cell id"""
    coefs = []
    for _ in range(8):
        cell_id, coef = divmod(cell_id, 26)
        coefs.append(coef)
    return "".join([chr(97 + coef) for coef in coefs][::-1])
