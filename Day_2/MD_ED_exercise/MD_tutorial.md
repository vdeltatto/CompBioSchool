# Molecular Dynamics of β-GAG⁺ with GROMACS

This tutorial describes the workflow used to prepare and simulate the β-GAG⁺ peptide in explicit water with GROMACS. The complete workflow will cover:

1. System preparation
2. Energy minimization
3. NVT equilibration
4. NPT equilibration
5. NVT production dynamics
6. Trajectory preparation
7. Essential Dynamics (ED) analysis

The tutorial covers the complete workflow from **system preparation** and **equilibration** to **production molecular dynamics**, **Essential Dynamics analysis**, construction of the **2D free-energy landscape**, and structural interpretation of its main conformational basins.

---

# References

This exercise is based primarily on the molecular-dynamics and Essential Dynamics analysis of aqueous cationic GAG⁺ described in:

1. M. Monti, M. Stener, and M. Aschi, **“A computational approach for modeling electronic circular dichroism of solvated chromophores,”** *Journal of Computational Chemistry* **43** (2022), 2023–2036.  
   https://doi.org/10.1002/jcc.27001

For additional explanations of the general GROMACS workflow, system preparation, equilibration, production dynamics, and trajectory analysis, the following tutorial is also useful:

2. J. A. Lemkul, **“Lysozyme in Water – GROMACS Tutorial.”**  
   http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html

For an introduction to molecular visualization and trajectory inspection with VMD:

3. **VMD Tutorial – Theoretical and Computational Biophysics Group, University of Illinois.**  
   https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/
   A PDF of the VMD Tutorial is also available in the folder of the exercise
   
For additional information about the MDAnalysis python package:

4. **Documentation MDAnalysis**  
   https://docs.mdanalysis.org/1.1.0/documentation_pages/overview.html

The GAG⁺ system and the Essential Dynamics interpretation used in this exercise follow Reference 1, while several aspects of the step-by-step GROMACS organization are inspired by Reference 2. Reference 3 provides additional guidance for using VMD.

---

# 1. System preparation

Before running a molecular dynamics simulation, the molecular system must be prepared. In this first stage we will:

1. generate the peptide topology;
2. define the simulation box;
3. solvate the system;
4. neutralize the total charge by adding a counterion.

The starting structure is provided in the `input` directory:

`gag-no-h.pdb`

The structure does not contain hydrogen atoms. These will be added automatically by GROMACS according to the selected force field and protonation states.

---

## 1.1 Prepare the system-generation directory

To keep all files generated during system preparation together, create a dedicated directory:

**mkdir syst_gen**

Move into the new directory:

**cd syst_gen**

Copy the starting PDB structure from the `input` directory:

**cp ../input/gag-no-h.pdb .**

All commands in the remainder of this section are run from inside `syst_gen`.

---

## 1.2 Generate the topology

The first step is to convert the PDB structure into a GROMACS-compatible structure and generate the corresponding topology.

Run:

**gmx pdb2gmx -f gag-no-h.pdb -ter**

The option `-ter` allows us to choose the protonation state of the N- and C-termini interactively.

During the execution, select:

```
Force field:
15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

Water model:
3: TIP3P  TIP 3-point

N-terminus:
0: GLY-NH3+

C-terminus:
2: GLY-COOH
```

The most important output files are:

conf.gro
topol.top

where:

- `conf.gro` contains the coordinates of the peptide after processing by `pdb2gmx`;
- `topol.top` contains the molecular topology and force-field information required for the simulation.

At this stage, it is good practice to inspect the generated structure and verify that the termini and hydrogen atoms are consistent with the intended protonation state.

---

## 1.3 Construct the simulation box

The peptide is placed at the center of a cubic simulation box with dimensions

2.7 × 2.7 × 2.7 nm³

using:

**gmx editconf -f conf.gro -box 2.7 2.7 2.7 -c -o conf-c-box.gro**

The relevant options are:

-f      input structure
-box    simulation-box dimensions in nm
-c      center the molecule in the box
-o      output structure

The resulting file is:

conf-c-box.gro

This structure contains the peptide centered inside the empty simulation box.

---

## 1.4 Solvate the system

The simulation box is filled with explicit water molecules using:

**gmx solvate -cp conf-c-box.gro -cs spc216.gro -o conf-c-box-solv.gro -p topol.top**

The options are:

-cp    structure to be solvated
-cs    solvent configuration
-o     output solvated structure
-p     topology file to update


For this system, solvation produces **639 water molecules**.

At the end of `topol.top`, the `[ molecules ]` section should therefore contain:

```
[ molecules ]
; Compound        #mols
Protein             1
SOL                639
```

The solvated structure is:

conf-c-box-solv.gro


---

## 1.5 Neutralize the simulation box

β-GAG⁺ has a total charge of **+1**. We therefore add one negatively charged counterion so that the total charge of the simulation box is zero.

In this case, one water molecule is replaced with one chloride ion, Cl⁻.

### 1.5.1 Generate the temporary input for ion insertion

First, generate a temporary `.tpr` file:

**gmx grompp -f ../input/minim.mdp -p topol.top -c conf-c-box-solv.gro -o ions.tpr -maxwarn 10**

`gmx grompp` combines:

- the simulation parameters from `minim.mdp`;
- the topology in `topol.top`;
- the coordinates in `conf-c-box-solv.gro`;

into the binary run-input file:

ions.tpr

This temporary file is required by `gmx genion`.

> **Note**
>
> `-maxwarn 10` allows `grompp` to proceed despite warnings. Warnings should still be read carefully before continuing, since not every warning is harmless.

---

### 1.5.2 Replace one water molecule with Cl⁻

Run:

**gmx genion -s ions.tpr -o conf-c-box-solv-ion.gro -p topol.top -nname CL -nn 1**

The options used here are:

-s       input .tpr file
-o       output coordinate file
-p       topology file to update
-nname   name of the negative ion
-nn      number of negative ions to add

GROMACS will ask which group should be used for ion substitution.

Select the water group:

Group    12 (          Water) has  1917 elements

One water molecule is then removed and replaced by one chloride ion.

Because the `-p topol.top` option was included, the topology is updated automatically.

The final `[ molecules ]` section should be:

[ molecules ]
; Compound        #mols
Protein             1
SOL                638
CL                  1

The final neutralized structure is:

conf-c-box-solv-ion.gro


---

## 1.6 Final system

At the end of the system-preparation stage, the simulation box contains:

1 β-GAG⁺ peptide
638 water molecules
1 Cl⁻ counterion

The system is electrically neutral and ready for the next stage:

Energy minimization

The coordinate file that will be used as the starting structure for energy minimization is:

conf-c-box-solv-ion.gro

---

# 2. Energy minimization

The solvated and neutralized system is now assembled. Before starting molecular dynamics, we first relax the structure to remove possible steric clashes or unfavorable local geometries. This step is called **energy minimization (EM)**.

As in the ion-addition step, `grompp` is used to combine the structure, topology, and simulation parameters into a binary `.tpr` file. This time, however, the `.tpr` file is passed to `mdrun`, which performs the energy minimization.

---

## 2.1 Prepare the working directory

At the end of the system-preparation stage, we are still inside the `syst_gen` directory. First, return to the main working directory:

**cd ../**

To keep the different stages of the simulation organized, create a directory for the equilibration workflow and a subdirectory for the energy minimization:

**mkdir equilibration**

**cd equilibration**

**mkdir energy_minimization**

**cd energy_minimization**

Copy the files required for the minimization into this directory:

**cp ../../input/minim.mdp ../../syst_gen/conf-c-box-solv-ion.gro ../../syst_gen/topol.top ../../input/run_gromacs.sh .**


The three main input files are:

- `minim.mdp`: simulation parameters for the energy minimization;
- `conf-c-box-solv-ion.gro`: solvated and neutralized starting structure;
- `topol.top`: topology of the complete system.

---

## 2.2 Generate the energy-minimization input

Generate the binary GROMACS input file using:

**gmx grompp -f minim.mdp -p topol.top -c conf-c-box-solv-ion.gro -o en_min.tpr**

The resulting file,

`en_min.tpr`

contains the information required by `mdrun` to perform the minimization.

---

## 2.3 Run the energy minimization

The energy minimization can be run using:

**gmx mdrun -deffnm en_min -v**

For the school, the same command is also provided inside the submission script `run_gromacs.sh`. It should be possible to run the minimization interactively (so just run the command on the terminal)

The most relevant options are:

- `-deffnm en_min`: uses `en_min` as the default name for the input and output files;
- `-ntmpi 1`: runs with one thread-MPI rank;
- `-ntomp 2`: uses two OpenMP threads;
- `-v`: prints the progress of the minimization to the terminal.

Because the input file is called `en_min.tpr`, `-deffnm en_min` is sufficient. If the `.tpr` file had a different name, it could instead be specified explicitly using the `-s` option.

After the calculation, the most important output files are:

```
en_min.log    ASCII-text log file of the minimization
en_min.edr    binary energy file
en_min.trr    full-precision trajectory
en_min.gro    energy-minimized structure
```

The structure in `en_min.gro` will be used as the starting point for the equilibration dynamics.

---

## 2.4 Check whether the minimization converged

At the end of the minimization, `mdrun` reports the final potential energy and forces. A successful run may end with output similar to:

```
Steepest Descents converged to Fmax < 1000 in 95 steps

Potential Energy  = -2.6801924e+04
Maximum force     =  9.9622791e+02 on atom 8
Norm of force     =  8.8724764e+01
```

The most important quantity to check is the **maximum force**, `Fmax`.

In the provided `minim.mdp` file, the convergence criterion is:

`emtol = 1000.0`

Therefore, the target is:

**Fmax < 1000 kJ mol⁻¹ nm⁻¹**

If GROMACS reports that the steepest-descent minimization has converged below this threshold, the minimization has reached the requested force tolerance.

The potential energy should also be inspected. It will normally decrease during minimization and is commonly negative for a solvated biomolecular system, but its absolute value depends strongly on the size and composition of the system. Therefore, the magnitude of the potential energy alone should **not** be used as a universal convergence criterion.

If the minimization stops with `Fmax` significantly larger than `emtol`, inspect the output carefully. Large residual forces can indicate unfavorable contacts or an inadequately minimized starting structure, and the minimization settings or initial geometry may need to be reconsidered.

---

## 2.5 Analyze the potential energy

The `.edr` file contains the energy terms collected during the minimization. These quantities can be extracted using `gmx energy`.

Run:

**gmx energy -f en_min.edr -o potential.xvg**

GROMACS will display a list of available energy terms. Select **Potential** and then enter `0` to finish the selection.

For the present setup, this can be entered as:

**10 0**

> **Note**
>
> The numerical index assigned to `Potential` can depend on the GROMACS version and on the system. Always check the list printed on the screen rather than assuming that it will be `10`.

The resulting file is:

`potential.xvg`

Plot `potential.xvg` using **cell 1 of the `equilibration_plots.ipynb` notebook located in: ../../notebooks/**.

**NOTE: the notebook assumes that the *.xvg is in the same folder, so copy the potential.xvg file there**:

cp potential.xvg ../../notebooks

During a successful energy minimization, the potential energy should decrease and approach a stable value as the structure relaxes.

---

---

# 3. NVT equilibration

Energy minimization gives us a reasonable starting structure by removing unfavorable contacts and local geometries. Before beginning the production dynamics, however, the system must be equilibrated.

Equilibration is commonly performed in two stages. First, we equilibrate the system in the **NVT ensemble**, where the Number of particles, Volume, and Temperature are kept constant. During this stage, the system is brought to the target temperature and the solvent is allowed to reorganize around the peptide.

To avoid large structural changes in the peptide while the solvent is still adapting to its surroundings, we apply **position restraints** to the heavy atoms of the peptide. These restraints are defined in the `posre.itp` file generated by `pdb2gmx`. The peptide atoms are still allowed to move, but deviations from their reference positions are energetically penalized.

The reference coordinates for the position restraints are provided to `grompp` through the `-r` option.

For this exercise, the NVT equilibration is run for **1 ns at 303 K**. The temperature should fluctuate around the target value without showing a systematic drift before proceeding to the next equilibration stage.

The file `nvt.mdp` is already provided and contains the parameters required for this stage.

---

## 3.1 Prepare the NVT equilibration directory

Leave the `energy_minimization` directory:

**cd ../**

Create a new directory for the NVT equilibration:

**mkdir nvt_equil**

**cd nvt_equil**

Copy the required files into the new directory:

**cp ../../input/nvt.mdp ../../syst_gen/topol.top ../../syst_gen/posre.itp ../energy_minimization/en_min.gro ../../input/run_gromacs.sh .**

The main files required for this stage are:

- `nvt.mdp`: parameters for the NVT equilibration;
- `topol.top`: topology of the complete system;
- `posre.itp`: position restraints for the peptide;
- `en_min.gro`: energy-minimized structure from the previous step;
- `run_gromacs.sh`: SLURM submission script used to run GROMACS.

---

## 3.2 Generate the NVT input file

Generate the binary input file using:

**gmx grompp -f nvt.mdp -c en_min.gro -p topol.top -r en_min.gro -o nvt.tpr -maxwarn 10**

The `-r en_min.gro` option provides the reference coordinates used by the position restraints.

The resulting file is:

`nvt.tpr`

---

## 3.3 Run the NVT equilibration

Open `run_gromacs.sh`. Comment out the command used for the previous stage and activate the NVT command:

**gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 2 -v**

For this exercise, the job requests:

```
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus=1
```

Submit the calculation with:

**sbatch run_gromacs.sh**

With these resources, the **1 ns NVT equilibration** should take approximately **1–2 minutes**.

After the calculation finishes, the final structure is written to:

`nvt.gro`

and the energy information is stored in:

`nvt.edr`

---

## 3.4 Check the temperature

After the NVT simulation is complete, extract the temperature from the energy file using:

**gmx energy -f nvt.edr -o temperature.xvg**

When GROMACS displays the list of available quantities, select **Temperature** and then terminate the selection with `0`.

For the present setup:

**15 0**

> **Note**
>
> The numerical index associated with `Temperature` may depend on the GROMACS version and on the system. Always check the list printed on the screen before making the selection.

The resulting file is:

`temperature.xvg`

Plot this file using **cell 2 of the `equilibration_plots.ipynb` notebook**.

**NOTE: the notebook assumes that the *.xvg is in the same folder, so copy the temperature.xvg file there**:

cp temperature.xvg ../../notebooks

During a properly equilibrated NVT simulation, the instantaneous temperature will fluctuate, but it should remain centered around the target value of **303 K** without a systematic drift.

---

---

# 4. NPT equilibration

The previous **NVT equilibration** brought the system to the target temperature. Before starting the production dynamics, we also need to equilibrate the **pressure and density** of the system.

This second equilibration stage is performed in the **NPT ensemble**, where the Number of particles, Pressure, and Temperature are controlled. The NPT ensemble is also called the **isothermal-isobaric ensemble**.

Compared with the NVT stage, the main addition in `npt.mdp` is the pressure-coupling section. In this exercise, pressure is controlled using the **C-rescale barostat**.

Pressure generally relaxes more slowly than temperature, so the NPT equilibration is run for **2 ns**.

Unlike the previous NVT equilibration, **position restraints are not applied during this stage**. This allows the complete solvated system, including the peptide, to relax while the box dimensions and density adjust to the target pressure.

The file `npt.mdp` is already provided and contains the parameters required for this stage.

---

## 4.1 Prepare the NPT equilibration directory

Leave the `nvt_equil` directory:

**cd ../**

Create a new directory for the NPT equilibration:

**mkdir npt_equil**

**cd npt_equil**

Copy the required files into the new directory:

**cp ../../input/npt.mdp ../../syst_gen/topol.top ../nvt_equil/nvt.gro ../nvt_equil/nvt.cpt ../../input/run_gromacs.sh .**

The main files required for this stage are:

- `npt.mdp`: parameters for the NPT equilibration;
- `topol.top`: topology of the complete system;
- `nvt.gro`: final structure obtained from the NVT equilibration;
- `nvt.cpt`: checkpoint containing the state of the system at the end of the NVT equilibration;
- `run_gromacs.sh`: SLURM submission script used to run GROMACS.

---

## 4.2 Generate the NPT input file

Generate the binary input file using:

**gmx grompp -f npt.mdp -c nvt.gro -p topol.top -t nvt.cpt -o npt.tpr -maxwarn 10**

The `-t nvt.cpt` option allows the NPT stage to continue from the state reached at the end of the NVT equilibration, including the velocities stored in the checkpoint.

The resulting file is:

`npt.tpr`

---

## 4.3 Run the NPT equilibration

Open `run_gromacs.sh`. Comment out the command used for the NVT stage and activate the NPT command:

**gmx mdrun -deffnm npt -ntmpi 1 -ntomp 2 -v**

For this exercise, the job requests:

```
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus=1
```

Submit the calculation with:

**sbatch run_gromacs.sh**

With these resources, the **2 ns NPT equilibration** should take approximately **2–3 minutes**.

After the calculation finishes, the final structure is written to:

`npt.gro`

the checkpoint is written to:

`npt.cpt`

and the energy information is stored in:

`npt.edr`

---

## 4.4 Check the pressure

Extract the pressure from the energy file using:

**gmx energy -f npt.edr -o pressure.xvg**

When GROMACS displays the list of available quantities, select **Pressure** and then terminate the selection with `0`.

For the present setup:

**17 0**

> **Note**
>
> The numerical index associated with `Pressure` may depend on the GROMACS version and on the system. Always check the list printed on the screen before making the selection.

The resulting file is:

`pressure.xvg`

Plot this file using **cell 3 of the `equilibration_plots.ipynb` notebook**.

**NOTE: the notebook assumes that the *.xvg is in the same folder, so copy the pressure.xvg file there**:

cp pressure.xvg ../../notebooks

---

## 4.5 Check the density

Extract the density using:

**gmx energy -f npt.edr -o density.xvg**

Select **Density** and then terminate the selection with `0`.

For the present setup:

**23 0**

> **Note**
>
> The numerical index associated with `Density` may depend on the GROMACS version and on the system. Always check the list printed on the screen before making the selection.

The resulting file is:

`density.xvg`

Plot this file using **cell 4 of the `equilibration_plots.ipynb` notebook**.

**NOTE: the notebook assumes that the *.xvg is in the same folder, so copy the density.xvg file there**:

cp density.xvg ../../notebooks

The pressure and density plots will be used to assess whether the system is sufficiently equilibrated before moving to the production dynamics.

---

---

# 5. NVT production dynamics

After completion of the NVT and NPT equilibration stages, the system has reached the desired temperature and an equilibrated density. We can now begin the **production molecular dynamics**, during which the trajectory used for subsequent analysis is collected.

In our workflow, the position restraints were already removed during the NPT equilibration. The production simulation therefore continues with the fully unrestrained system.

For the production run, we return to the **NVT ensemble**. The simulation box obtained at the end of the NPT equilibration is kept fixed, while the temperature is maintained at **303 K**.

As in the transition from NVT to NPT equilibration, the checkpoint file from the previous stage is used when generating the new `.tpr` file. This allows the production simulation to start from the full-precision coordinates and velocities reached at the end of NPT equilibration.

The file `prod.mdp` is already provided and contains the parameters required for the production simulation. For this exercise, the production trajectory has a total length of **50 ns**.

---

## 5.1 Prepare the production directory

At the end of the NPT equilibration, we are inside `equilibration/npt_equil`. Return to the main working directory:

**cd ../../**

Create a directory for the production simulation:

**mkdir production**

**cd production**

Copy the required files:

**cp ../input/prod.mdp ../equilibration/npt_equil/npt.gro ../equilibration/npt_equil/npt.cpt ../syst_gen/topol.top ../input/run_gromacs.sh .**

The main files required for this stage are:

- `prod.mdp`: parameters for the NVT production simulation;
- `npt.gro`: final structure obtained from the NPT equilibration;
- `npt.cpt`: checkpoint containing the state of the system at the end of NPT equilibration;
- `topol.top`: topology of the complete system;
- `run_gromacs.sh`: SLURM submission script used to run GROMACS.

---

## 5.2 Generate the production input file

Generate the binary input file using:

**gmx grompp -f prod.mdp -p topol.top -c npt.gro -t npt.cpt -o prod.tpr -maxwarn 10**

The `-t npt.cpt` option transfers the full-precision coordinates and velocities from the end of the NPT equilibration into the new production input file.

The resulting file is:

`prod.tpr`

---

## 5.3 Run the production dynamics

Open `run_gromacs.sh` and activate the final GROMACS command:

**gmx mdrun -deffnm prod -ntmpi 1 -ntomp 2 -v -nb gpu -pme gpu -pin on**

The production simulation uses the same computational resources:

```
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus=1
```

The command uses one MPI rank and two OpenMP threads. The short-range non-bonded interactions and PME electrostatics are offloaded to the GPU.

**IMPORTANT**
A **50 ns** production simulation takes approximately **30 minutes** with the resources used for this exercise. To avoid spending this time during the school, the completed production trajectory is already provided in the `output` directory, together with the reference outputs from the equilibration stages.

You therefore do **not** need to submit the 50 ns production calculation during the exercise.

Copy the provided trajectory into the production directory:

**cp ../output/production/prod.xtc .**

The file `prod.xtc` contains the coordinate trajectory that will be used for the subsequent analysis.

---

---

# 6. Trajectory preparation

Before performing the Essential Dynamics analysis, the production trajectory must be prepared. Because molecular dynamics simulations employ **periodic boundary conditions (PBC)**, the peptide can appear to cross the simulation-box boundaries during the trajectory. This is physically harmless, but it makes visualization and structural analysis inconvenient.

We therefore generate a new trajectory in which the peptide is kept at the center of the simulation box and molecules are reconstructed consistently across periodic boundaries.

---

## 6.1 Center the production trajectory

Run:

**gmx trjconv -f prod.xtc -s prod.tpr -center -boxcenter rect -pbc mol -ur compact -o prod_c.xtc**

GROMACS will ask for two selections.

For centering, select:

**1: Protein**

For the output trajectory, select:

**0: System**

The resulting trajectory is:

`prod_c.xtc`

The `-center` option centers the selected group, while `-pbc mol` keeps molecules whole with respect to the periodic boundaries. The `-ur compact` option writes the system using a compact representation of the periodic unit cell.


---

---

# 7. Essential Dynamics analysis

Essential Dynamics (ED), also commonly referred to as principal component analysis (PCA) of a molecular-dynamics trajectory, is used to identify the dominant collective motions sampled during the simulation.

Instead of analyzing the motion of every Cartesian coordinate independently, ED describes the trajectory in terms of a new set of coordinates called **eigenvectors** or **principal components**. The corresponding **eigenvalues** indicate how much of the structural fluctuation is associated with each eigenvector.

For this analysis, we focus on the **backbone atoms of the peptide**, since they provide a compact description of the motions associated with changes in the peptide conformation and secondary structure.

---

## 7.1 Calculate the covariance matrix

Run:

**gmx covar -f prod_c.xtc -s prod.tpr**

`gmx covar` calculates the covariance matrix of the atomic positional fluctuations and diagonalizes it to obtain the eigenvalues and eigenvectors.

Before calculating the covariance matrix, GROMACS performs a **least-squares fit** of the trajectory to a reference structure. This fitting step removes the overall translation and rotation of the peptide, allowing the analysis to focus on its **internal motions** rather than on the motion of the entire molecule through space.

GROMACS will ask for two selections.

First, select the group to use for the least-squares fit:

**4: Backbone**

Then select the group for which the covariance matrix should be calculated:

**4: Backbone**

Using the same backbone atoms for both operations ensures that the covariance analysis describes the internal fluctuations of the peptide backbone after removal of its overall roto-translational motion.

Among the files generated by `gmx covar`, the two most important for the following analysis are:

`eigenval.xvg`

`eigenvec.trr`

The file `eigenval.xvg` contains the eigenvalues of the covariance matrix, while `eigenvec.trr` contains the corresponding eigenvectors.

---

## 7.2 Inspect the eigenvalue spectrum

Plot `eigenval.xvg` using **cell 1 of the `ED_ana_backbone.ipynb` notebook that can be found in the ../notebooks folder**.

**NOTE: the notebook expects the files to be in the same folder. You can copy the notebook in the current production folder**

Each eigenvalue measures the amount of structural fluctuation associated with the corresponding eigenvector. The eigenvectors are ordered from the largest eigenvalue to the smallest, so the first few components describe the largest-amplitude collective motions sampled during the trajectory.

For a two-dimensional representation of the Essential Dynamics space to be informative, the first two eigenvectors should account for a substantial fraction of the total fluctuation.

This contribution can be expressed as:

**(λ₁ + λ₂) / Σᵢ λᵢ**

where `λ₁` and `λ₂` are the first two eigenvalues and the denominator is the sum of all eigenvalues obtained from the covariance analysis.

The first cell of the notebook calculates this quantity automatically and prints:

`First 2 eigenvalues cover XX.XX% of the total variance`

This provides a quantitative measure of how much of the backbone fluctuation is retained when the trajectory is reduced to the first two essential coordinates.

There is no universal percentage that the first two components must reproduce. The important point is to check whether they capture a sufficiently large fraction of the total fluctuation for the two-dimensional projection to provide a meaningful representation of the dominant conformational motions.

---

## 7.3 Project the trajectory onto the first two eigenvectors

We now project every frame of the production trajectory onto the first two eigenvectors:

**gmx anaeig -f prod_c.xtc -s prod.tpr -first 1 -last 2 -2d**

The eigenvectors calculated in the previous step are read from the default file:

`eigenvec.trr`

When prompted for the fitting/reference and analysis groups, select:

**4: Backbone**

for both selections, consistently with the covariance analysis.

The resulting file is:

`2dproj.xvg`

For each trajectory frame, this file contains its projection onto the first and second eigenvectors. These two coordinates define the **two-dimensional essential plane** that will be used to visualize the conformational space sampled by the peptide.

---

## 7.4 Construct the 2D free-energy landscape

The file `2dproj.xvg` contains the projection of every trajectory frame onto the first two eigenvectors. Each frame is therefore represented by a point in the two-dimensional essential plane:

**Proj-1, Proj-2**

We can use the density of points in this plane to identify the regions of conformational space that are sampled most frequently during the simulation.

The corresponding analysis is performed using the cells under **“2. 2D free-energy histogram + low-energy structures”** in the `ED_ana_backbone.ipynb` notebook.

---

### 7.4.1 Build the 2D histogram

The first step is to divide the essential plane into a two-dimensional grid and count how many trajectory frames fall into each bin.

In the provided notebook, the plane is divided into:

```
bin_x = 50
bin_y = 50
```

giving a **50 × 50 histogram**.

For the 50 ns trajectory used in this exercise, coordinates were saved every 2 ps, giving **25,001 trajectory frames**, including the initial frame. Therefore, the notebook uses:

```
dime = 25001
```

The histogram population of each bin provides an estimate of how frequently that region of the essential plane is sampled.

Highly populated bins correspond to frequently visited conformations, whereas sparsely populated regions correspond to less frequently sampled conformations.

---

### 7.4.2 Convert the histogram into a free-energy landscape

The probability distribution in the essential plane can be converted into a **relative free-energy landscape** using the Boltzmann relation:

**ΔGᵢ = -RT ln(Pᵢ / Pmax)**

where:

- `Pᵢ` is the population of a given histogram bin;
- `Pmax` is the population of the most populated bin;
- `R` is the gas constant;
- `T` is the temperature of the production simulation.

Before running this part of the notebook, make sure that the value of `temp` corresponds to the temperature used in `prod.mdp`.

Because the probabilities are normalized to the most populated bin, the most populated region is assigned:

**ΔG = 0 kJ mol⁻¹**

and all other bins are expressed relative to this minimum.

The notebook writes the resulting free-energy surface to:

`2d_histo.txt`

This file contains three columns:

```
Proj-1    Proj-2    ΔG
```

The free-energy values are calculated in J mol⁻¹ and converted to kJ mol⁻¹ when the landscape is plotted.

> **Note**
>
> This is a **projected relative free-energy landscape** obtained from the population of the first two essential coordinates. It should therefore be interpreted as a representation of the conformational sampling in this two-dimensional space, rather than as an absolute free energy of the system.

---

### 7.4.3 Identify the most populated regions

We are also interested in identifying representative configurations belonging to the most populated regions of the landscape.

In the notebook, bins with a relative free energy of:

**ΔG ≤ 2.5 kJ mol⁻¹**

are selected as low-energy, highly populated regions.

For each selected histogram bin, the notebook searches `2dproj.xvg` for an actual trajectory frame whose Proj-1 and Proj-2 coordinates lie sufficiently close to the center of that bin.

The corresponding trajectory time, projections, and free energy are written to:

`low_en_structs.txt`

with the columns:

```
Time    Proj-1    Proj-2    ΔG
```

Because the trajectory contains one frame every **2 ps**, the projection index can be directly converted into the corresponding simulation time.

These points are useful because they correspond to **real configurations sampled during the MD trajectory**, rather than to artificial coordinates located only at the centers of the histogram bins.

Strictly speaking, the single bin with `ΔG = 0` is the **most probable** region. The configurations selected with the 2.5 kJ mol⁻¹ cutoff represent a broader set of **low-free-energy, highly populated regions** of the landscape.

---

### 7.4.4 Plot the essential-dynamics landscape

Finally, run the plotting cell under **“3. Plot histogram + low-energy configurations”** in the `ED_ana_backbone.ipynb` notebook.

The notebook plots the free-energy surface obtained from `2d_histo.txt` as a contour map and superimposes the low-energy projections from `low_en_structs.txt`.

The axes correspond to:

- **Proj-1**: projection onto the first eigenvector;
- **Proj-2**: projection onto the second eigenvector.

The color scale represents the relative free energy in **kJ mol⁻¹**, while the black points indicate trajectory configurations belonging to the low-energy regions selected above.

The resulting figure is saved as:

`2dhi.pdf`

> **Note: your landscape may appear mirrored**
>
> When comparing your free-energy landscape with the reference figure, you may find that one of the projection axes appears inverted, so that the landscape looks mirrored horizontally or vertically.
>
> This is completely normal. The sign of an eigenvector obtained from diagonalization of the covariance matrix is arbitrary: an eigenvector **v** and the corresponding vector **−v** describe exactly the same collective motion. Therefore, Proj-1 and/or Proj-2 may have the opposite sign in an independent analysis, even when the underlying conformational dynamics are equivalent.
>
> A mirrored landscape is therefore **not an error**. The positions of the basins may be reflected with respect to one axis, but their relative populations, free energies, and structural interpretation remain unchanged. If desired for visual comparison only, a projection coordinate can be multiplied by `-1`, but this is not necessary for the analysis.

---

## 7.5 Interpret the 2D landscape

The two-dimensional landscape provides a compact representation of the conformational space explored by the peptide during the simulation.

Regions with **low free energy** correspond to highly populated conformational states. These appear as basins or minima in the landscape. If several distinct minima are present, the trajectory has sampled multiple preferred backbone conformations.

Regions of higher free energy are less frequently populated and separate the major conformational basins in the projected space.

The position of a minimum along Proj-1 and Proj-2 does not by itself provide a direct structural interpretation. To understand what distinguishes the different regions of the landscape, representative structures from the corresponding basins must be extracted from the trajectory and compared.

> **Important**
>
> The landscape is based only on the first two eigenvectors. Even when these components capture a substantial fraction of the total fluctuation, additional motions are present in the remaining dimensions. The 2D landscape should therefore be viewed as a reduced representation of the dominant conformational dynamics.

---

## 7.6 Select representative structures from the basins

The final part of the `ED_ana_backbone.ipynb` notebook groups the low-energy configurations contained in `low_en_structs.txt` into spatial clusters in the essential plane. These clusters provide an operational way of identifying the different low-free-energy basins sampled by the trajectory.

The clustering is performed using **DBSCAN**, with the parameters:

```
eps = 0.03
min_samples = 3
```

The notebook prints the number of clusters that were identified. The value of `eps` can be adjusted if the clustering does not reproduce the basins that are visible in the 2D free-energy landscape.

For each cluster, the notebook then selects a number of configurations spanning the range of free energies sampled within that basin.

The number of structures selected per basin is controlled by:

```
n_reps = 5
```

You can increase or decrease `n_reps` depending on how many structures you want to inspect.

The selected structures are written to:

`selected_structures.txt`

with the columns:

```
cluster    time    proj1    proj2    dG
```

The `time` column gives the simulation time in **ps** corresponding to the selected trajectory frame.

The notebook also plots the clustered low-energy configurations and highlights the selected representative structures. This allows you to check visually that the selected structures are distributed across the different basins.

---

## 7.7 Extract the representative structures from the trajectory

The configurations listed in `selected_structures.txt` correspond to actual frames of the centered production trajectory.

For each selected time, extract the corresponding structure using:

**gmx trjconv -f prod_c.xtc -s prod.tpr -dump time -o conf_time.gro**

Replace `time` with the value, in ps, reported in `selected_structures.txt`.

For example, if a selected structure corresponds to 19226 ps, the command would be:

**gmx trjconv -f prod_c.xtc -s prod.tpr -dump 19226 -o conf_19226.gro**

When prompted for the output group, select:

**1: Protein**

if you only want to extract the peptide structure.

Repeat this procedure for the representative configurations that you want to inspect.

The resulting `.gro` files can then be opened in **VMD** and compared to determine the structural differences between the conformational basins identified in the Essential Dynamics landscape.

To extract the representative structures automatically, use the script:

`extract_structures.sh`

available in:

`../notebooks`

The script reads the simulation times of the representative low-energy conformations and extracts the corresponding frames from the centered trajectory using `gmx trjconv`.

Internally, it runs commands of the form:

**gmx trjconv -f prod_c.xtc -s prod.tpr -dump "$time" -o "$out"**

The extracted `.gro` structures are automatically organized into separate folders according to the free-energy basin to which they belong.

---

## 7.8 Relate the free-energy basins to peptide conformations

Open the representative `.gro` structures extracted from the different basins in **VMD** and compare the peptide backbone conformations.

A ready-to-use VMD visualization state is available in:

`../output/production/ED_analysis_backbone`

The file is called:

`selected_low_en_confs.vmd`

Copy it into your current folder, since it needs to find the `low_en_structures` folder generated in Section 7.7:

**cp ../output/production/ED_analysis_backbone/selected_low_en_confs.vmd .**

**NOTE:**
The file contains **absolute paths** pointing to the folder where it was originally created, which will not match your own directory. Before loading it in VMD, these paths must be updated to point to your own `low_en_structures` folder.

First, find your current working directory:

**pwd**

Then open the file with a text editor, for example `vi`:

**vi selected_low_en_confs.vmd**

Inside the file there are three lines starting with `mol new`, one per basin. Search for them by typing:

**/mol new**

and pressing **Enter**, then press **n** to jump to the next match. Each of these lines ends with a path of the form:

`.../low_en_structures/basin_N/selected_basin_N.gro type gro first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all`

Replace everything **before** `low_en_structures` (not included) with the output of `pwd`, so that each line reads:

`mol new <your pwd output>/low_en_structures/basin_N/selected_basin_N.gro type gro first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all`

Repeat this for all three `mol new` lines (`basin_0`, `basin_1`, `basin_2`), then save and quit (`:wq`).

> **Note**
>
> If you are comfortable with `vi`, a single substitution command updates all three lines at once instead of editing them one by one — replace `<old path prefix>` with the path shown in the file (everything up to, but not including, `low_en_structures`):
>
> `:%s#<old path prefix>#<your pwd output>#g`

Open **VMD**, then go to:

**File → Load Visualization State**

and load `selected_low_en_confs.vmd`.

Three `.gro` files will be loaded:

`selected_basin_0.gro`  
`selected_basin_1.gro`  
`selected_basin_2.gro`

Each file contains representative structures from one of the identified free-energy basins.

For clarity, only one basin is displayed initially, while the other two are hidden. In the provided visualization state, `selected_basin_2.gro` is displayed first.

To look at the different representative configurations belonging to the **same basin**, use the trajectory controls in the **VMD Main** window and click the **step-forward button** (`|>`) to advance one configuration at a time.

To display a **different basin**:

1. Go to **Graphics → Representations**.
2. In the **Selected Molecule** menu, choose the `.gro` file corresponding to the basin you want to inspect.
3. In the representation list, locate the row corresponding to the **CPK** representation.
4. **Double-click the CPK row** to make that basin visible.
5. Select the basin that was previously displayed and **double-click its CPK row** to hide it.

You can then return to the **VMD Main** window and use the **step-forward button** (`|>`) again to browse through the representative configurations of the newly selected basin.

For additional guidance on using VMD, see the official VMD tutorial:

https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/

If you have issues, you can also try loading the `.gro` files in the `low_en_structures` folder following the VMD tutorial (using File -> New Molecule).

You should observe a predominance of **pPII-like** conformations, followed by **β-strand-like** conformations and a smaller population of **right-handed α (αR)-like** conformations. In this small peptide, pPII- and αR-like backbone geometries can give rise to more globally compact structures, whereas β-like conformations are generally more open and extended.

These conformational classes do **not map one-to-one onto individual Essential Dynamics basins**. Both of the main low-energy basins contain a mixture of backbone conformations, although their relative populations differ. This reflects the highly dynamic nature of the peptide: local backbone conformations continuously interconvert while the molecule explores nearby regions of the global conformational landscape.

The two main basins are connected by a relatively small free-energy barrier, approximately **3–4.5 kJ/mol** from visual inspection of the calculated landscape. This is consistent with **rapid interconversion between the sampled conformational states**, rather than with two rigid and well-separated structures.

For the trajectory analyzed in this exercise, the global minimum lies within the more compact region of the landscape, while the more extended β-like configurations are somewhat less populated. The important picture is therefore that GAG⁺ exists as a **dynamic conformational ensemble**, with several backbone geometries accessible within closely connected low-free-energy regions.

This interpretation is consistent with the reference study. At lower temperature, GAG⁺ shows a higher propensity for the **pPII-like** state, which is enthalpically favored, whereas increasing the temperature shifts the conformational population toward the more **β-strand-like** state, which is entropically favored. This temperature-dependent redistribution is supported by computational studies and experimental spectroscopic data.

For a complete discussion of the GAG⁺ conformational landscape, its temperature dependence, and its relationship with the experimental circular-dichroism spectra, see:

M. Monti, M. Stener, and M. Aschi, *Journal of Computational Chemistry* **43** (2022), 2023–2036.  
https://doi.org/10.1002/jcc.27001

---

# 8. Additional structural analysis

The Essential Dynamics landscape identifies the main conformational states sampled by GAG⁺. We can now complement this picture with additional structural descriptors. We first examine the **backbone dihedral angles** using a Ramachandran plot, and then ask how the different clusters identified in Section 7.6 differ in **compactness** and in their **interaction with the surrounding water**. We also track two complementary diagnostics, RMSD and RMSF, that describe how the peptide's structure evolves and fluctuates over the trajectory as a whole, independently of any specific cluster.

For these analyses we will use **MDAnalysis**, a Python library for the analysis of molecular-dynamics trajectories.

All calculations are performed on the centered production trajectory:

`prod_c.xtc`

using the production topology:

`prod.tpr`

The trajectory can be loaded in MDAnalysis with:

```
import MDAnalysis as mda

u = mda.Universe("prod.tpr", "prod_c.xtc")
```

The same trajectory frames used for the Essential Dynamics analysis are therefore used here, which will later allow us to relate each structural descriptor to the position of the peptide in the 2D free-energy landscape.

---

## 8.1 Ramachandran plot

The backbone dihedral angles **φ** (`C(i-1)-N-CA-C`) and **ψ** (`N-CA-C-N(i+1)`) define the local backbone conformation of a peptide. Plotting φ against ψ therefore allows us to relate the structures sampled during the trajectory to familiar regions of Ramachandran space, including **pPII-like**, **β-strand-like**, and **right-handed α (αR)-like** conformations.

For a peptide containing only three residues, the terminal residues cannot each provide a complete `(φ, ψ)` pair: φ requires the carbonyl carbon of the previous residue, while ψ requires the nitrogen of the following residue. Therefore, **GLY1** and **GLY3** each lack one of the two angles, and only the central residue, **ALA2**, has a complete `(φ, ψ)` pair.

The Ramachandran analysis is performed with `MDAnalysis.analysis.dihedrals.Ramachandran`:

```
from MDAnalysis.analysis.dihedrals import Ramachandran

protein = u.select_atoms("protein")
rama = Ramachandran(protein).run()
angles = rama.results.angles
```

When the full `protein` selection is passed to `Ramachandran`, MDAnalysis automatically excludes residues for which a complete `(φ, ψ)` pair cannot be defined. For this system, the returned angles therefore correspond to **ALA2** only; MDAnalysis may print a warning indicating that terminal residues have been removed from the analysis.

The notebook plots the `(φ, ψ)` values sampled during the production trajectory. The low-energy representative structures identified in Section 7.6 are also marked on the plot and colored according to their cluster.

For orientation, the plot includes approximate reference points for the canonical regions:

- **pPII:** φ ≈ -75°, ψ ≈ 145°;
- **β-strand / extended:** φ ≈ -130°, ψ ≈ 130°;
- **right-handed α (αR):** φ ≈ -60°, ψ ≈ -45°.

These values are **approximate visual references**, not rigorous boundaries between conformational states. In particular, sampling the αR region means that the local backbone dihedral angles are helix-like; for a three-residue peptide, this does **not** imply formation of a complete α-helix.

---

## 8.2 Root mean square deviation (RMSD)

The **root mean square deviation**, RMSD, measures how far the peptide's structure has drifted, atom-by-atom after an optimal least-squares superposition, from a reference structure — here the first production frame (t = 0 ns).

RMSD is computed for the whole peptide using the `RMSD` class from `MDAnalysis.analysis.rms`:

```
from MDAnalysis.analysis import rms

rmsd_protein = rms.RMSD(u, u, select="protein").run().results.rmsd[:, 2]
```

This class performs the fit itself — both translation and rotation are removed before the RMSD is computed, using the standard Theobald QCP algorithm — so no separate alignment step is needed.

Unlike Rg, RMSD does not describe an intrinsic property of a given structure: it describes deviation from one specific, somewhat arbitrary reference frame. For this reason RMSD is tracked here as a standalone diagnostic of the peptide's overall conformational drift and stability over the trajectory, and is **not** included in the final per-cluster comparison table in Section 8.8.

The low-energy representative structures from Section 7.6 are marked on the RMSD time series (colored by their cluster), the same convention used throughout the rest of this section. Time is reported in **ns**.

---

## 8.3 Root mean square fluctuation (RMSF)

While RMSD tracks how the whole structure drifts over time relative to one reference frame, the **root mean square fluctuation**, RMSF, asks a complementary, per-atom question: how much does each individual atom move around its own average position over the trajectory? This gives a flexibility profile of the peptide, not a time series.

Computing RMSF properly requires first aligning every frame to the trajectory's own **average structure** (not just a single reference frame), following [MDAnalysis's standard recipe](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rms.html):

```
from MDAnalysis.analysis import align, rms

average = align.AverageStructure(u2, u2, select="protein and not name H*", ref_frame=0).run()
align.AlignTraj(u2, average.results.universe, select="protein and not name H*", in_memory=True).run()

heavy = u2.select_atoms("protein and not name H*")
rmsf = rms.RMSF(heavy).run().results.rmsf
```

> **Note**
>
> `AlignTraj(..., in_memory=True)` modifies the coordinates of the universe it is given, in place. Run this on a separate copy of the universe (`u2` above, not the `u` used everywhere else in this section) so it doesn't affect the other analyses, which rely on the original, unaligned trajectory.

RMSF is computed **per heavy atom** (`protein and not name H*`), not per residue: for this 3-residue peptide there are only 14 heavy atoms, so a per-atom profile (labeled by residue and atom name, e.g. `GLY1-N`, `ALA2-CB`) is more informative than collapsing it to 3 per-residue values.

Like RMSD, RMSF is an ensemble/per-atom quantity rather than a per-frame one — there is no meaningful "RMSF of one structure" — so it is not part of the final per-cluster comparison table in Section 8.8 either.

---

## 8.4 Radius of gyration

The **radius of gyration**, Rg, provides a simple measure of the compactness of a molecular structure. A larger Rg generally corresponds to a more extended conformation, while a smaller Rg indicates a more compact structure.

For GAG⁺, we will calculate Rg for two atom selections:

- the **whole peptide**;
- the **peptide backbone**.

In MDAnalysis these groups can be selected as:

```
protein = u.select_atoms("protein")
backbone = u.select_atoms("protein and backbone")
```

For every trajectory frame, calculate:

```
protein.radius_of_gyration()
backbone.radius_of_gyration()
```

and store the corresponding simulation time (in **ns**).

The notebook plots the two Rg time series (as two stacked panels, whole peptide and backbone). MDAnalysis reports distances in **Å**, so the radius of gyration will also be obtained in Å.

The comparison between the two definitions is useful because the full-peptide Rg includes all atoms, while the backbone Rg focuses more directly on the overall peptide conformation.

The low-energy representative structures identified in Section 7.6 are marked directly on both panels (colored by their cluster), so their Rg values can be read off and compared against the overall mean ± standard deviation, rather than splitting the whole trajectory into basins.

This provides a direct question to explore:

**Do the different conformational clusters differ systematically in Rg — is one more extended than the others?**

Section 8.8's table lets you check this cluster by cluster, using the representative structures identified in Section 7.6.

---

## 8.5 Peptide-water radial distribution function

The **radial distribution function**, RDF, describes how the density of one group of atoms varies as a function of the distance from another group.

Here, we use the RDF to characterize the organization of water around the peptide.

A useful general peptide-water RDF can be obtained between:

- **peptide heavy atoms**;
- **water oxygen atoms**.

For example:

```
protein_heavy = u.select_atoms("protein and not name H*")
water_O = u.select_atoms("resname SOL and name OW")
```

> **Note**
>
> Check the atom names in your topology before running the analysis. In the present GROMACS system the water residue is expected to be called `SOL`, but the oxygen atom name should always be verified.

The RDF can be calculated using the `InterRDF` class from `MDAnalysis.analysis.rdf`.

The resulting function, **g(r)**, provides information about the hydration structure around the peptide. Peaks in g(r) indicate distances at which water molecules are preferentially found relative to peptide atoms, while g(r) approaching 1 corresponds to a bulk-like distribution.

For the first analysis, calculate the RDF using **all frames of the production trajectory**.

The analysis can then be repeated separately for each cluster's low-energy representative frames identified in Section 7.6 (a modest but physically meaningful sample specific to that conformational family, rather than a lone structure or the whole unclustered trajectory).

Comparing these cluster-resolved RDFs allows us to ask whether different conformational families produce different hydration environments. Differences may appear, for example, in the position or intensity of the first hydration-shell peak.

Because each cluster contributes a different number of representative frames, each cluster-specific RDF should be **normalized independently** (one `InterRDF` instance per cluster) rather than compared using raw pair counts.

---

## 8.6 Peptide-water hydrogen bonds

Water does not interact with the peptide only through its spatial distribution. We can also quantify the more specific directional interactions formed through **hydrogen bonds**.

For this analysis we are interested only in hydrogen bonds formed **between the peptide and water**, excluding peptide-peptide and water-water hydrogen bonds.

MDAnalysis provides the `HydrogenBondAnalysis` class in:

```
MDAnalysis.analysis.hydrogenbonds.hbond_analysis
```

The analysis can be restricted to peptide-water interactions using the two selections:

```
protein
resname SOL
```

through the `between` option.

The standard MDAnalysis hydrogen-bond definition uses geometric criteria based on the donor-acceptor distance and donor-hydrogen-acceptor angle. The default criteria are a donor-acceptor distance of **3.0 Å** and a donor-hydrogen-acceptor angle of at least **150°** — the exact MDAnalysis defaults, used here explicitly rather than left implicit.

For every trajectory frame, calculate the **number of peptide-water hydrogen bonds** and plot this quantity as a function of time (in **ns**).

This is done **two ways**: **All** (every donor/acceptor on the peptide, termini included) and **Backbone** (restricted to the mainchain amide N-H and carbonyl C=O only). Do not rely on `guess_hydrogens()`/`guess_acceptors()` for the "All" selections without checking them first — for this topology the auto-guesser silently drops the ALA backbone carbonyl O and the C-terminal OT oxygen from the acceptor list, undercounting real hydrogen bonds. Explicit wildcard selections (`name O*`, `name H*`) avoid this. For "Backbone", `protein and name H` isolates exactly the two standard backbone amide hydrogens for this system (excluding the N-terminal `H1/H2/H3` and the C-terminal `HO`), and `backbone and name O*` gives the three backbone carbonyl oxygens (excluding the C-terminal `OT`) — always verify this against your own topology's atom names rather than assuming it.

As with Rg, the low-energy representative structures identified in Section 7.6 are marked on both time series (colored by their cluster).

This allows us to investigate whether the different conformational clusters differ in the number of hydrogen bonds they form with the solvent, and whether that difference is driven by the termini or by the backbone itself.

The interpretation should be statistical: individual hydrogen bonds continuously form and break during the trajectory, so we are interested in the **distribution and average number of peptide-water hydrogen bonds** rather than in one particular hydrogen bond observed in a single structure.

---

## 8.7 Solvent-accessible surface area (SASA)

The **solvent-accessible surface area**, SASA, gives a geometric measure of how much of the peptide's surface is exposed to solvent overall — complementing the hydrogen-bond count, which only counts specific directional interactions. A more compact conformation generally buries more surface and has a lower SASA.

MDAnalysis has no built-in SASA calculator, so this uses the external `freesasa` package:

```
pip install freesasa
```

For each frame, the current peptide coordinates are written to a temporary PDB file and passed to `freesasa.calc`:

```
import freesasa, tempfile, os

freesasa.setVerbosity(freesasa.silent)

def frame_sasa(frame_index):
    u.trajectory[frame_index]
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        protein.write(tmp.name)
        tmp_path = tmp.name
    try:
        return freesasa.calc(freesasa.Structure(tmp_path)).totalArea()
    finally:
        os.remove(tmp_path)
```

> **Note**
>
> `freesasa.setVerbosity(freesasa.silent)` suppresses a per-atom warning about the C-terminal `OT` atom name (from `pdb2gmx`'s neutral-COOH naming), which is not a standard PDB atom name — `freesasa` still correctly identifies it as an oxygen with a standard radius, so the warning can be safely ignored.

Computing SASA this way for every frame takes about a minute over the full 50 ns trajectory. As with the other descriptors, the low-energy representative structures are marked on the time series (in **ns**), colored by their cluster.

---

## 8.8 Connect the structural descriptors with the free-energy landscape

The most informative part of these analyses comes from combining them with the low-energy representative structures identified in Section 7.6.

For each representative structure, the notebook builds one table combining its position in the landscape (proj1, proj2, ΔG) with its Rg (whole peptide and backbone), H-bond count (all and backbone-only), and SASA — grouped by cluster and sorted by increasing ΔG within each cluster.

For each cluster, the notebook also reports:

- the mean ± standard deviation of each quantity, computed over that cluster's representative structures;
- the minimum and maximum value observed within the cluster.

These are printed alongside the overall mean ± standard deviation (over the whole trajectory) reported earlier for each quantity, so a cluster's typical values can be compared directly against the general trend.

**RMSD and RMSF are not part of this table.** RMSD measures deviation from one specific reference frame rather than an intrinsic property of a structure, and RMSF is a per-atom, whole-trajectory quantity with no well-defined value for a single representative structure. Both remain useful as their own standalone diagnostics (Sections 8.2–8.3).

----------------------

For the trajectory analyzed in this exercise, **basin 0 is strongly enriched in pPII-like conformations**, whereas **basin 1 contains a more heterogeneous mixture of pPII-like (predominant) and αR-like (small %) conformations**. The **β-like population** is relatively small but present in both basins.

> **Note**
>
> A minor αR-like population is also acknowledged in the literature for small, unfolded GxG-type peptides, but it is comparatively less explored than the pPII/β equilibrium. This is largely because varying the temperature mainly redistributes population between pPII and the predominant β state, so most studies — including the reference work this exercise is based on — focus on the two-state pPII↔β picture rather than resolving finer sub-populations like the αR-like conformations identified here.

The clearest difference between the two basins is their **compactness**. Basin 0 has a larger average radius of gyration than basin 1, both for the whole peptide and for the backbone. Basin 1 therefore contains, on average, more compact structures. This is consistent with the presence of both pPII and αR-like conformations, which can produce a more compact global geometry in this small peptide.

In contrast, the differences in **SASA** and in the average number of **peptide-water hydrogen bonds** are small. The two basins therefore do not differ strongly in their overall solvent exposure or simply in the number of hydrogen bonds they form with water.

This illustrates an important point when analyzing MD trajectories: **different conformational states are not necessarily distinguished by every structural descriptor**. In this case, the main difference between the basins is geometric, whereas simple global measures of hydration are much more similar.

The descriptors used here (RDF, total hydrogen-bond count, SASA) are all **global, ensemble-averaged** measures of hydration — they summarize the whole peptide's interaction with water into a single number per frame, but do not resolve *where* or *how* that hydration happens. More detailed, finer-grained solvent properties would be needed to probe subtler differences between the backbone conformations, for example:

- **which specific peptide groups are hydrated** — e.g., whether water preferentially solvates the charged termini versus the backbone carbonyls/amides differently in one conformation than another;
- **hydrogen-bond geometries and lifetimes** — not just how many H-bonds form, but how linear/strong they are and how long they persist, since two conformations could form the same *number* of H-bonds with very different stability;
- **the orientation and organization of nearby water molecules** — e.g., whether one conformation induces a more ordered hydration shell, or bridges specific peptide groups through additional water-mediated contacts, that a bulk RDF or total H-bond count would not reveal on its own.

These would require frame-by-frame or atom/residue-resolved analysis of the hydration shell, going beyond the global averages computed in this exercise.

---


