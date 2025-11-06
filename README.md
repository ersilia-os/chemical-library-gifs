# Chemical Library GIFs

Produce simple GIFs to quickly visualize chemical libraries.

This simple Python package produces a simple GIF file for a collection of compounds. Below is a visualization of a sample chemical library of 1000 Enamine REAL compounds:

![Sample Library](assets/random_1000_1.gif)

Molecules are displayed using the excellent [mol2ps/mol2svg](https://homepage.univie.ac.at/norbert.haider/cheminf/mol2ps.html) tool developed by [Norbert Haider](https://homepage.univie.ac.at/norbert.haider/).

## Installation

To get started, clone the repository:

```bash
git clone https://github.com/ersilia-os/chemical-library-gifs.git
cd chemical-library-gifs
```

To get started, create a Conda environment:

```bash
conda create -n chemgifs python=3.12
conda activate chemgifs
```

You will need Cairo installed:

```bash
conda install -c conda-forge cairo
```

Then install the package using pip:

```bash
pip install .
```

## Usage

Simply pass a one-column CSV file, with a header (`smiles`).

```bash
chemgifs -i molecules.csv -o molecules.gif -c yellow -s 512 -d 200
```

Accepted color names correspond to the official Ersilia [color palette](https://ersilia.gitbook.io/ersilia-book/styles/brand-guidelines). By default, the size of the image will be 512 pixels and the duration of each frame 200 ms.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
