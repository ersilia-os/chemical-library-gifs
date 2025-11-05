import os
import tempfile
from rdkit import Chem

root = os.path.dirname(os.path.abspath(__file__))

def get_rgb_color(name):
    if name == "white":
        return (255, 255, 255)
    elif name == "mint":
        return (190, 230, 180)
    elif name == "pink":
        return (220, 160, 220)
    elif name == "purple":
        return (170, 150, 250)
    elif name == "orange":
        return (250, 160, 140)
    elif name == "yellow":
        return (250, 215, 130)
    elif name == "blue":
        return (140, 200, 250)
    elif name == "gray":
        return (210, 210, 210)
    else:
        raise Exception("Invalid color name")


def get_mol_svg(name, smiles, svg_file, color_name):
    mol2svg_exec = os.path.join(root, "..", "tools", "mol2svg")
    mol = Chem.MolFromSmiles(mol)
    if mol is None:
        return
    


if __name__ == "__main__":
    pass