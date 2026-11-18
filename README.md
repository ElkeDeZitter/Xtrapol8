<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/ElkeDeZitter/Xtrapol8">
    <img src="Xtrapol8_logo.png" alt="Logo" width="405" height="120">
  </a>

  <p align="center">
    Structural analysis of low-occupancy states
  </p>
</p>


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about">About</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
        <li><a href="#updating">Updating</a></li>
      </ul>
    </li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#command-line-usage">Command line usage</a></li>
        <li><a href="#graphical-user-interface-gui-usage">Graphical User Interface (GUI) usage</a></li>
        <li><a href="#tutorial">Tutorial</a></li>
      </ul>
    </li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
# About

<!-- [![Product Name Screen Shot][product-screenshot]](https://example.com) -->

Xtrapol8 is software for the calculation of Fourier difference maps and extrapolated structure factor amplitudes (ESFAs), and for the estimation of the occupancy of a triggered state that is present in the crystal at low occupancy. It is based on the cctbx toolbox and uses several CCP4 programs. To run Xtrapol8, you need CCP4 (including Coot-1) installed with a valid license. Phenix is optional but recommended (version 1.18 or higher), see [Prerequisites](#prerequisites).

An introduction to Fourier difference maps, structure factor extrapolation and the software itself can be found on the CCP4 YouTube channel: [CCP4 Study Weekend 2022](https://www.youtube.com/watch?v=H21L3-M_ge4) and [CCP4 Study Weekend 2024](https://www.youtube.com/watch?v=t6-7OBwtVts).

If Xtrapol8 was useful in your project, please cite us:<br />
De Zitter, E., Coquelle, N., Oeser, P., Barends, T. R. M., Colletier, J.-P., [Xtrapol8 enables automatic elucidation of low-occupancy intermediate-states in crystallographic studies](https://doi.org/10.1038/s42003-022-03575-7), Communications Biology, 5(1), 640 (2022).<br />

Don't forget to cite the software on which we build:
* Liebschner, D. et al., [Macromolecular structure determination using X-rays, neutrons and electrons: recent developments in Phenix](https://journals.iucr.org/d/issues/2019/10/00/di5033/index.html), Acta Crystallogr. D Struct. Biol. 75, 861–877 (2019).
* Agirre, J. et al., [The CCP4 suite: integrative software for macromolecular crystallography](https://journals.iucr.org/d/issues/2023/06/00/ai5011/index.html), Acta Crystallogr. D Struct. Biol. 79, 449–461 (2023).
* Emsley, P. et al., [Features and development of Coot](https://journals.iucr.org/d/issues/2010/04/00/ba5144/index.html), Acta Crystallogr. D Biol. Crystallogr. 66, 486–501 (2010).
* Grosse-Kunstleve, R. W. et al., [The Computational Crystallography Toolbox: crystallographic algorithms in a reusable software framework](https://journals.iucr.org/paper?ks0118), J. Appl. Crystallogr. 35, 126–136 (2002).

See the [release notes](https://github.com/ElkeDeZitter/Xtrapol8/releases) and the manuals for more information and new features.

<!-- GETTING STARTED -->
# Getting Started

This version of Xtrapol8 runs with the Python that comes with CCP4 (`ccp4-python`). The minimal requirement is CCP4 with Coot 1. We also recommend installing Phenix, so that refinement can be done with the Phenix programs. Once CCP4, Coot and (optionally) Phenix are installed, follow the steps below to make sure that Xtrapol8 can find everything it needs.

## Prerequisites

You need the [CCP4](https://www.ccp4.ac.uk) software suite with Coot-1, and optionally the [Phenix](https://www.phenix-online.org) software suite. Both have to be set up for use from the command line. Below is a detailed description of how to set up CCP4 and Phenix for Xtrapol8 in a terminal on macOS (on Linux, the steps are the same, only the installation paths differ):

1. Add Phenix to your PATH by sourcing the file `phenix_env.sh`, which is located in the Phenix installation folder. Replace the version number with the one you have installed (you can use Finder to check which Phenix version you have and where `phenix_env.sh` is located):
   ```bash
   source /Applications/phenix-2.2.1-6174/phenix_env.sh
   ```
2. Add CCP4 to your PATH by sourcing `ccp4.setup-sh` in the same terminal. Again, replace the version with the one you have installed:
   ```bash
   source /Applications/ccp4-9/bin/ccp4.setup-sh
   ```
3. Check that the Phenix, CCP4 and Coot-1 programs can be found (each command should return a path):
   ```bash
   which phenix.refine
   which scaleit
   which ccp4-python
   which coot-1
   ```

To avoid repeating steps 1 and 2 each time you open a new terminal, you can add them to your `~/.zprofile` (default shell on macOS), `~/.bash_profile` or `~/.bashrc` file. Open the file with your favorite editor (e.g. `nano`):
   ```bash
   nano ~/.zprofile
   ```
and add the following lines (with your own installation paths):
   ```bash
   # Setting up Phenix and CCP4
   source /Applications/phenix-2.2.1-6174/phenix_env.sh
   source /Applications/ccp4-9/bin/ccp4.setup-sh
   ```

Coot-1 is installed together with CCP4, so it is set up once CCP4 is set up. A Coot-1 installed separately also works, as long as `coot-1` can be found from the command line. This version of Xtrapol8 is not compatible with Coot-0.x versions.

If you do not have Phenix, use the CCP4 programs for refinement by setting `refinement.reciprocal_space=refmac5` and `refinement.real_space=coot`.

## Installation

Run all commands below in a terminal in which CCP4 has been set up (see [Prerequisites](#prerequisites)).

1. Clone the Xtrapol8 repository and go into the Xtrapol8 directory:
   ```bash
   git clone https://github.com/ElkeDeZitter/Xtrapol8.git
   cd Xtrapol8
   ```
   Alternatively, download the repository without git: on the GitHub page, select the `develop` branch, click on "Code" -> "Download ZIP" and unzip it in an easily accessible location (e.g. your home directory). Keep all files together in the same directory. Then go into the unzipped directory and skip step 2:
   ```bash
   cd </path/to>/Xtrapol8
   ```

2. Switch to the development branch (only needed if you used `git clone`):
   ```bash
   git checkout develop
   ```

3. Install or update the two dependencies (pypubsub and fonttools):
   ```bash
   ccp4-python -m pip install pypubsub
   ccp4-python -m pip install --upgrade fonttools
   ```

4. Install Xtrapol8 (the `.` refers to the current directory, so run this from inside the Xtrapol8 directory):
   ```bash
   ccp4-python -m pip install . --no-deps
   ```

5. Check the installation by printing the Xtrapol8 version:
   ```bash
   xtrapol8-version
   ```

## Updating

To update an existing installation, go into the Xtrapol8 directory, get the latest version and reinstall it:
   ```bash
   cd </path/to>/Xtrapol8
   git pull
   ccp4-python -m pip install . --no-deps
   ```
If you downloaded the ZIP file, download the new version and repeat steps 4 and 5 of the [Installation](#installation).

<!-- USAGE EXAMPLES -->
# Usage
## Command line usage

After installing Xtrapol8, you can run it with the command:
   ```bash
   xtrapol8
   ```

When Xtrapol8 is run without arguments, it prints the full list of parameters that can be set via an input file or the command line. An example input file with all parameters that can be changed is available as [Xtrapol8_example.phil](Xtrapol8_example.phil), and an input file with only the minimal parameters as [Xtrapol8_minimal.phil](Xtrapol8_minimal.phil).

* Example using an input file:
  1. Edit Xtrapol8.phil with your favorite text editor:
     ```bash
     nano Xtrapol8.phil
     ```
  2. Run Xtrapol8:
     ```bash
     xtrapol8 Xtrapol8.phil
     ```
* Example using command line arguments only:
     ```bash
     xtrapol8 input.reference_mtz=hiephiep.mtz input.triggered_mtz=hieperdepiep.mtz input.reference_pdb=hoera.pdb input.additional_files=jeej.cif input.additional_files=another.cif occupancies.list_occ=0.1,0.3,0.5 f_and_maps.f_extrapolated_and_maps=qfextr,qfgenick map_explorer.peak_integration_floor=3.5 map_explorer.peak_detection_threshold=4 output.outdir=fancy_party
     ```
* Example using an input file and command line arguments:
  1. Edit Xtrapol8.phil with your favorite text editor:
     ```bash
     nano Xtrapol8.phil
     ```
  2. Run Xtrapol8 with additional arguments. Arguments are read from left to right, so a parameter given later overrides the same parameter given earlier (here, the value in the input file):
     ```bash
     xtrapol8 Xtrapol8.phil refinement.phenix_keywords.refine.cycles=3
     ```

* Other command line programs:
  1. Show the Xtrapol8 version:
     ```bash
     xtrapol8-version
     ```
  2. Run the difference map analysis on an already completed Xtrapol8 run (run without arguments to get the help message):
     ```bash
     xtrapol8-difference-map-analysis
     ```
  3. Run the distance analysis on an already completed Xtrapol8 run (run without arguments to get the help message):
     ```bash
     xtrapol8-distance-analysis
     ```
  4. Calculate a distance difference matrix (run without arguments to get the help message):
     ```bash
     xtrapol8-ddm
     ```
  5. Run Refiner (run without arguments to get the help message). Refiner might still contain bugs in this version:
     ```bash
     xtrapol8-refiner
     ```
  6. Open the graphical user interface. The GUI might still contain bugs in this version:
     ```bash
     xtrapol8-gui
     ```

_For more information, please have a look at the [Command Line Manual](Xtrapol8_Command_Line_Manual.pdf)._

## Graphical User Interface (GUI) usage

Xtrapol8 can also be run via its GUI, which is launched from the command line:
  ```bash
  xtrapol8-gui
  ```

* The GUI can also use a pre-made input file:
  * The input file can be loaded directly when launching XtrapolG8:
    ```bash
    xtrapol8-gui Xtrapol8.phil
    ```
  * Or the input file can be loaded from within XtrapolG8 via the File menu:<br />
    "File" -> "Open Phil"

* XtrapolG8 can also be used to inspect the results of a previous run, regardless of whether it was run via the GUI or the command line:<br />
  "File" -> "Load Results"

_For more information, please have a look at the [GUI Manual](XtrapolG8_Manual.pdf)._

## Tutorial

A [basic tutorial](https://github.com/ElkeDeZitter/Xtrapol8-tutorial) is available, in which you can reproduce the results of one of the examples from the publication.

<!-- LICENSE -->
# License

Xtrapol8 is published under the [MIT License](LICENSE).<br />
Xtrapol8 makes use of adapted cctbx modules, for which we additionally refer to the [cctbx license](https://github.com/cctbx/cctbx_project/blob/master/LICENSE.txt).<br />
The user is responsible for legally obtaining and installing [CCP4](https://www.ccp4.ac.uk/download/), [Coot](https://www.ccp4.ac.uk/download/) and [Phenix](https://www.phenix-online.org/download/).

<!-- CONTACT -->
# Contact

Elke De Zitter - elke.de-zitter@ibs.fr

<!-- Nicolas Coquelle -  nicolas.coquelle@ibs.fr -->

<!-- Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de -->

Jacques-Philippe Colletier - jacques-Philippe.colletier@ibs.fr

[Institut de Biologie Structurale - DYNAMOP](https://www.ibs.fr/research/research-groups/dynamics-and-kinetics-of-molecular-processes-group-m-weik/)
