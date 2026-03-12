import os
import csv
import shutil
import tempfile
import subprocess
import cairosvg
import argparse

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from PIL import Image
import platform

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
        raise ValueError(
            f"Invalid color name: {name}. Valid options are: {', '.join(colors.keys())}"
        )


def resolve_mol2svg_exec():
    system = platform.system()
    machine = platform.machine()

    if system == "Linux" and machine == "x86_64":
        return os.path.abspath(
            os.path.join(root, "tools", "linux_x86", "mol2svg")
        )

    if system == "Darwin" and machine == "x86_64":
        return os.path.abspath(
            os.path.join(root, "tools", "macosx_x86", "mol2svg")
        )

    if system == "Darwin" and machine == "arm64":
        return os.path.abspath(
            os.path.join(root, "tools", "macosx_arm64", "mol2svg")
        )

    raise Exception("mol2svg binary not found for this platform")


def get_raw_mol_svg(name, smiles, output_dir, color):
    mol2svg_exec = resolve_mol2svg_exec()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    AllChem.Compute2DCoords(mol)
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
    ns = ""
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
    bg_rect = ET.Element(
        f"{ns}rect",
        {
            "x": str(x - dx),
            "y": str(y - dy),
            "width": str(size),
            "height": str(size),
            "fill": fill_color,
        },
    )
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
    input_svg = os.path.join(output_dir, name + ".svg")
    output_png = os.path.join(output_dir, name + ".png")
    cairosvg.svg2png(
        url=input_svg, write_to=output_png, output_width=size, output_height=size
    )
    os.remove(input_svg)


def make_grid_image(png_files, n_rows, n_cols, cell_size, background_rgb):
    grid_w = n_cols * cell_size
    grid_h = n_rows * cell_size
    grid = Image.new("RGBA", (grid_w, grid_h), background_rgb + (255,))
    for idx, p in enumerate(png_files):
        row = idx // n_cols
        col = idx % n_cols
        with Image.open(p) as img:
            img = img.convert("RGBA").resize((cell_size, cell_size))
        grid.paste(img, (col * cell_size, row * cell_size))
    return grid


def pngs_to_gif(png_files, output_gif, duration, n_rows, n_cols, cell_size, background_rgb, loop=0):
    cells_per_frame = n_rows * n_cols
    frames = []
    for i in range(0, len(png_files), cells_per_frame):
        chunk = png_files[i : i + cells_per_frame]
        frames.append(make_grid_image(chunk, n_rows, n_cols, cell_size, background_rgb))

    frames[0].save(
        output_gif,
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=loop,
        optimize=False,
        disposal=2,
    )


def run(input_csv, output_file, color_name, size, duration_ms, n_rows=1, n_cols=1, max_mols=None):
    background_rgb = get_rgb_color(color_name)
    is_png = output_file.lower().endswith(".png")
    tmp_dir = tempfile.mkdtemp("ersilia-")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    smiles_list = read_smiles(input_csv)
    if is_png:
        smiles_list = smiles_list[:n_rows * n_cols]
    elif max_mols is not None:
        smiles_list = smiles_list[:max_mols]
    for i, smiles in tqdm(enumerate(smiles_list)):
        name = "mol_{0}".format(str(i).zfill(6))
        get_mol_svg(name, smiles, tmp_dir, color_name)
        svg_to_png(name, tmp_dir, size=size)
    png_files = sorted(os.path.join(tmp_dir, fn) for fn in os.listdir(tmp_dir))
    if is_png:
        frame = make_grid_image(png_files, n_rows, n_cols, size, background_rgb)
        frame.save(output_file)
    else:
        pngs_to_gif(png_files, output_file, duration=duration_ms, n_rows=n_rows, n_cols=n_cols, cell_size=size, background_rgb=background_rgb)
    shutil.rmtree(tmp_dir)


def main():
    args = argparse.ArgumentParser(
        description="Generate animated GIF of molecules from SMILES."
    )

    args.add_argument("-i", "--input_csv", type=str, help="Input CSV file with SMILES.")
    args.add_argument("-o", "--output_gif", type=str, help="Output file (.gif or .png).")
    args.add_argument(
        "-c", "--color", type=str, default="white", help="Background color name."
    )
    args.add_argument(
        "-s", "--size", type=int, default=512, help="Size of each frame in pixels."
    )
    args.add_argument(
        "-d",
        "--duration",
        type=int,
        default=200,
        help="Duration of each frame in milliseconds.",
    )
    args.add_argument(
        "--n_rows",
        type=int,
        default=1,
        help="Number of rows in the grid per frame.",
    )
    args.add_argument(
        "--n_cols",
        type=int,
        default=1,
        help="Number of columns in the grid per frame.",
    )
    args.add_argument(
        "--max_mols",
        type=int,
        default=None,
        help="Maximum number of molecules to process.",
    )

    parsed_args = args.parse_args()
    run(
        input_csv=parsed_args.input_csv,
        output_file=parsed_args.output_gif,
        color_name=parsed_args.color,
        size=parsed_args.size,
        duration_ms=parsed_args.duration,
        n_rows=parsed_args.n_rows,
        n_cols=parsed_args.n_cols,
        max_mols=parsed_args.max_mols,
    )


if __name__ == "__main__":
    main()
