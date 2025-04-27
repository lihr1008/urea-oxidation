# urea-oxidation
Description of folders:
1. AutoCal:
(DFT) We used the atomate2 workflow to call VASP software to perform structure optimization and static calculations on adsorption substrates and adsorption structures with different doping structures to obtain their energy and charge information.
(MD) Firstly, we use GROMACS software to generate initial configurations with different metal compositions.Then we use the lammps_interface package to generate a structure file of type lammps and an input parameter file for the UFF force field which add energy minimization and NVT simulation. After that, we call LAMMPS software based on subprocess to automate the simulation calculations.
2. MD Preprocess: We process 20,000 MD trajectories with different metal compositions, starting from step 10,000, and uniformly sample 100 conformations for each trajectory. Centering on each metal, we obtain multi-metallic local structural models, resulting in a total of 6,000.
3. MD Properties: Based on DFT calculation results, we perform statistical averaging on the multi-metallic local structural models obtained from conformation sampling to acquire the adsorption energy and charge transfer corresponding to each metal composition.
4. Model:
`pretrain.py`: The premodel where the input is the metal composition from MD, and the output is the corresponding catalytic properties (adsorption energy and charge transfer).
`train.py`: Theory-practice integration model. The input is the experimental metal composition, which first goes through the pre-model in pretrain to obtain its catalytic properties, then passes through fully connected layers to get its overpotential. Additionally, the SHAP algorithm analyzes the importance of the proportions of five metals in the test set on the prediction results, excluding non-high-entropy data for plotting (paint_shap.py).
`grid_search.py`: We perform a grid search with a step size of 1% over the chemical space for each metal composition ranging from 5% to 35%, using the theory-practice integration model to predict its overpotential, followed by stratified sampling for experimental validation.
