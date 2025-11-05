import os
import tempfile
import subprocess
from rdkit import Chem

root = os.path.dirname(os.path.abspath(__file__))

def get_rgb_color(name):
    colors = {
        "white": (255, 255, 255),
        "mint": (190, 230, 180),
        "pink": (220, 160, 220),
        "purple": (170, 150, 250),
        "orange": (250, 160, 140),
        "yellow": (250, 215, 130),
        "blue": (140, 200, 250),
        "gray": (210, 210, 210),
    }
    if name in colors:
        return colors[name]
    else:
        raise ValueError(f"Invalid color name: {name}. Valid options are: {', '.join(colors.keys())}")


def get_mol_svg(name, smiles, output_dir, color_name):
    color = get_rgb_color(color_name)
    mol2svg_exec = os.path.join(root, "..", "tools", "mol2svg")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    sdf_file = os.path.join(output_dir, f"{name}.sdf")
    svg_file = os.path.join(output_dir, f"{name}.svg")
    with Chem.SDWriter(sdf_file) as writer:
        writer.write(mol)
    cmd = f"{mol2svg_exec} --output=svg --bgcolor={color[0]},{color[1]},{color[2]} {sdf_file} > {svg_file}"
    subprocess.run(cmd, shell=True, check=True)
    os.remove(sdf_file)



if __name__ == "__main__":
    get_mol_svg("benzene", "c1ccccc1", ".", "pink")