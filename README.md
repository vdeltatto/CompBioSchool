# Summer School on Theoretical and Computational Biophysics

Lecture slides and hands-on material for the school, organised by day.

## Structure

```
Day_1/   Introduction, statistical mechanics, and a hands-on notebook on the
         Verlet integrator and thermostats.
Day_2/   Force fields, hydrated systems, and the MD workflow. `MD_ED_exercise/`
         holds a full GROMACS tutorial (`MD_tutorial.md`) with its `input/`
         files, `output/` trajectories and analyses, and `notebooks/` for
         equilibration plots and essential-dynamics analysis.
Day_3/   Enhanced sampling. `US_Exercise/` is an umbrella-sampling exercise on
         the Müller-Brown potential, built on the local `molsim` package.
Day_4/   Density functional theory.
Day_5/   Overview of unsupervised methods; clustering and Markov state modeling
```

Day 2 has extra requirements beyond the Python environment — see
[GROMACS, VMD and exercise data](#gromacs-vmd-and-exercise-data).

## Setting up the Python environment

Everything in the notebooks runs on **Python 3.10+**. Pick either route.

### Option A — conda / miniconda

If you don't have conda, install
[Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) first. Then:

```bash
conda create -n compbio python=3.12
conda activate compbio
pip install -r requirements.txt
```

### Option B — venv

First check your Python version with `python3 --version`:

- **macOS**: the built-in `python3` is often an old system version (e.g. 3.9)
  that will **not** work. If it's too old, install a newer one with Homebrew
  (`brew install python@3.12`) and use `python3.12` below.
- **Linux**: the system `python3` is usually fine.
- **Windows**: get Python from [python.org](https://www.python.org/downloads/),
  ticking "Add to PATH" during install.

Then, from this directory:

```bash
python3 -m venv venv
source venv/bin/activate     # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

### Day 3: the `molsim` package

The Day 3 notebook imports `molsim`, a local package that is not on PyPI.
With your environment active, install it from this directory:

```bash
pip install -e Day_3/US_Exercise
```

### Running the notebooks

Register the environment as a Jupyter kernel:

```bash
python -m ipykernel install --user --name compbio --display-name "Python (compbio)"
```

Then open the notebooks in VS Code or with `jupyter lab`, selecting the
**"Python (compbio)"** kernel.

## GROMACS, VMD and exercise data

The Day 2 exercise runs the simulations with **GROMACS** and visualises the
results with **VMD**. Both are separate from the Python environment above and
must be installed on their own.

### GROMACS

Install `gmx` with whichever fits your system:

```bash
conda install -c conda-forge gromacs   # any platform with conda
brew install gromacs                   # macOS
sudo apt install gromacs               # Ubuntu / Debian
```

There is no native Windows build — use WSL, or run the simulations on the
cluster. To compile from source (needed for GPU support), follow the
[GROMACS installation guide](https://manual.gromacs.org/current/install-guide/index.html).

Check it worked:

```bash
gmx --version
```

Note that `notebooks/extract_structures.sh` calls `gmx trjconv` directly from
the shell, so it needs GROMACS but *not* the Python environment.

### VMD

Download from the [VMD site](https://www.ks.uiuc.edu/Research/vmd/) (free
registration required) and follow `MD_ED_exercise/vmd-tutorial.pdf`.

### Exercise data

The trajectories and reference outputs for the exercise of Day 2 are too large for this
repository, so you can download them from Google Drive at this [link](https://drive.google.com/drive/folders/1E-5biEIZL0NDnWPQt-0zna1TtQSTFfiQ?usp=drive_link).

Download the folder and place it under `Day_2/`, keeping its name, so that the
paths used in `MD_tutorial.md` resolve correctly.
