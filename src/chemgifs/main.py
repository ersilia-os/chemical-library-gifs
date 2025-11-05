import os
import csv
import shutil
import tempfile
import subprocess
import cairosvg

from tqdm import tqdm
from rdkit import Chem
from PIL import Image

import xml.etree.ElementTree as ET

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


def get_raw_mol_svg(name, smiles, output_dir, color):
    mol2svg_exec = os.path.join(root, "..", "tools", "mol2svg")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    sdf_file = os.path.join(output_dir, f"{name}.sdf")
    svg_file = os.path.join(output_dir, f"{name}.svg")
    with Chem.SDWriter(sdf_file) as writer:
        writer.write(mol)
    cmd = f"{mol2svg_exec} --autoscale=on --output=svg --bgcolor={color[0]},{color[1]},{color[2]} {sdf_file} > {svg_file}"
    subprocess.run(cmd, shell=True, check=True)
    os.remove(sdf_file)


def square_and_center_svg(input_svg, output_svg, background_rgb):
    tree = ET.parse(input_svg)
    root = tree.getroot()
    ns = ''
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0] + "}"
    vb = root.attrib.get("viewBox")
    if vb:
        x, y, w, h = map(float, vb.split())
    else:
        w = float(root.attrib.get("width", 0))
        h = float(root.attrib.get("height", 0))
        x, y = 0, 0
        root.set("viewBox", f"{x} {y} {w} {h}")
    size = max(w, h)
    dx = (size - w) / 2
    dy = (size - h) / 2
    root.set("viewBox", f"{x - dx} {y - dy} {size} {size}")
    root.set("width", str(size))
    root.set("height", str(size))
    fill_color = f"rgb({background_rgb[0]}, {background_rgb[1]}, {background_rgb[2]})"
    bg_rect = ET.Element(f"{ns}rect", {
        "x": str(x - dx),
        "y": str(y - dy),
        "width": str(size),
        "height": str(size),
        "fill": fill_color
    })
    root.insert(0, bg_rect)
    tree.write(output_svg)


def get_mol_svg(name, smiles, output_dir, color_name):
    color = get_rgb_color(color_name)
    svg_file = os.path.join(output_dir, f"{name}.svg")
    get_raw_mol_svg(name, smiles, output_dir, color)
    square_and_center_svg(svg_file, svg_file, color)


def read_smiles(input_csv):
    with open(input_csv, "r") as f:
        reader = csv.reader(f)
        smiles_list = []
        next(reader)
        for r in reader:
            smiles_list += [r[0]]
    return smiles_list


def svg_to_png(name, output_dir, size):
    input_svg = os.path.join(output_dir, name+".svg")
    output_png = os.path.join(output_dir, name+".png")
    cairosvg.svg2png(url=input_svg, write_to=output_png, output_width=size, output_height=size)
    os.remove(input_svg)


def pngs_to_gif(png_files, output_gif, duration, loop=0):
    images = [Image.open(p) for p in png_files]
    images = [img.convert("RGBA") for img in images]
    w, h = images[0].size
    images = [img.resize((w, h)) for img in images]
    images[0].save(
        output_gif,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=loop,
        optimize=False,
        disposal=2
    )


def main(input_csv, output_gif, color_name, size=512, duration_ms=200):
    tmp_dir = tempfile.mkdtemp("ersilia-")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    smiles_list = read_smiles(input_csv)
    for i, smiles in tqdm(enumerate(smiles_list)):
        name = "mol_{0}".format(str(i).zfill(6))
        get_mol_svg(name, smiles, tmp_dir, color_name)
        svg_to_png(name, tmp_dir, size=size)
    png_files = []
    for fn in os.listdir(tmp_dir):
        png_files += [os.path.join(tmp_dir, fn)]
    png_files = sorted(png_files)
    pngs_to_gif(png_files, output_gif, duration=duration_ms)
    shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    input_csv = os.path.join(root, "..", "..", "data", "example_compounds.csv")
    color_name = "white"
    main(input_csv, "test.gif", color_name)    