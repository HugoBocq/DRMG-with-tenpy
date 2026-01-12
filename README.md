# DRMG-with-tenpy
Perform DRMG simulations on different semi-infinite geometries to define the effect of three-body scalar chirality term on a spin Hamiltonian.

### Prerequirements 
- Install the TeNPy Python library which is a tensor-network package to simulate strongly correlated quantum systems. (https://tenpy.readthedocs.io/en/latest/) 

### Structure of the repository
The example file 'simulation.yml' defines the algorithm parameters, the initial conditions/state, the values measured/computed.
The lattice geometry and the interactions are defined within the python file "model.py" in a class that is imported and called in 'simulation.yml'. 
Result folder contains .pkl files with simulation results for given setup/geometry. There are plotted and analyzed in jupyter-notebook 'plot_XXX.ipynb'.

### Running the code
Run the code from the command line in the appropriate environment where TeNPy is installed by entering: 'tenpy-run simulation.yml'

### Example data plotted and analyzed
- triangular ladder with exchange and scalar chirality Hamiltonian for spin s=1/2. 
