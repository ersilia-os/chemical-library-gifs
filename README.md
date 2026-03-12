# Chemical Library GIFs

Produce simple GIFs to quickly visualize chemical libraries.

This simple Python package produces a GIF file to display a collection of compounds. Below is a visualization of a sample chemical library of 1000 Enamine REAL compounds:

![Sample Library](assets/random_1000_1.gif)

Molecules are displayed using the excellent [mol2ps/mol2svg](https://homepage.univie.ac.at/norbert.haider/cheminf/mol2ps.html) tool developed by [Norbert Haider](https://homepage.univie.ac.at/norbert.haider/).

## Installation

First, create a Conda environment with Cairo:

```bash
conda create -n chemgifs python=3.12
conda activate chemgifs
conda install -c conda-forge cairo
```

Then install directly from GitHub:

```bash
pip install https://github.com/ersilia-os/chemical-library-gifs/archive/refs/heads/main.zip
```

Or clone and install in editable mode:

```bash
git clone https://github.com/ersilia-os/chemical-library-gifs.git
cd chemical-library-gifs
pip install -e .
```

## Usage

Simply pass a one-column CSV file with a header (`smiles`). The output can be a `.gif` (animated) or a `.png` (single static image).

```bash
chemgifs -i molecules.csv -o molecules.gif
```

### Arguments

| Argument | Short | Default | Description |
|---|---|---|---|
| `--input_csv` | `-i` | — | Input CSV file with a `smiles` column |
| `--output_gif` | `-o` | — | Output file. Use `.gif` for animation or `.png` for a static image |
| `--color` | `-c` | `white` | Background color name |
| `--size` | `-s` | `512` | Size of each molecule cell in pixels |
| `--duration` | `-d` | `200` | Duration of each GIF frame in milliseconds (ignored for PNG) |
| `--n_rows` | | `1` | Number of rows in the grid per frame |
| `--n_cols` | | `1` | Number of columns in the grid per frame |
| `--max_mols` | | all | Maximum number of molecules to process (GIF only) |

Accepted color names correspond to the official Ersilia [color palette](https://ersilia.gitbook.io/ersilia-book/styles/brand-guidelines): `white`, `mint`, `pink`, `purple`, `orange`, `yellow`, `blue`, `gray`.

### Examples

Animated GIF, one molecule per frame:
```bash
chemgifs -i molecules.csv -o molecules.gif -c yellow -s 512 -d 200
```

Animated GIF with a 4×8 grid per frame (32 molecules per frame):
```bash
chemgifs -i molecules.csv -o molecules.gif --n_rows 4 --n_cols 8
```

Static PNG with the first 6 molecules in a 2×3 grid:
```bash
chemgifs -i molecules.csv -o molecules.png --n_rows 2 --n_cols 3
```

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
