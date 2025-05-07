from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from lammps_interface.lammps_main import LammpsSimulation
from lammps_interface.structure_data import from_CIF, write_CIF
import sys
from pathlib import Path
import shutil

parser = CifParser(f'box.cif')
structure = parser.parse_structures(primitive=False)

class Parameters:
    def __init__(self, cif):
        # File options
        self.cif_file = cif
        self.output_cif = False
        self.output_raspa = False

        # Force field options
        self.force_field = 'UFF'
        self.mol_ff = None
        self.h_bonding = False
        self.dreid_bond_type = 'harmonic'
        self.fix_metal = False

        # Simulation options
        self.minimize = False
        self.bulk_moduli = False
        self.thermal_scaling = False
        self.npt = False
        self.nvt = False
        self.cutoff = 12.5
        self.replication = None
        self.orthogonalize = False
        self.random_vel = False
        self.dump_dcd = 0
        self.dump_xyz = 0
        self.dump_lammpstrj = 0
        self.restart = False

        # Parameter options
        self.tol = 0.4
        self.neighbour_size = 5
        self.iter_count = 10
        self.max_dev = 0.01
        self.temp = 298.0
        self.pressure = 1.0
        self.nprodstp = 200000
        self.neqstp = 200000

        # Molecule insertion options
        self.insert_molecule = ""
        self.deposit = 0

    def show(self):
        for v in vars(self):
            print('%-15s: %s' % (v, getattr(self, v)))

cif_path = Path('box.cif').resolve()
run_dir = Path(f"lammps_run_{cif_path.stem}").resolve()
run_dir.mkdir(parents=True, exist_ok=True)

target_cif = run_dir / cif_path.name
if not target_cif.exists():
    shutil.copy(cif_path, target_cif)

par = Parameters(str(target_cif))
par.show()
sim = LammpsSimulation(par)
cell, graph = from_CIF(par.cif_file)
sim.set_cell(cell)
sim.set_graph(graph)
sim.split_graph()
sim.assign_force_fields()
sim.compute_simulation_size()
sim.merge_graphs()
if par.output_cif:
    print("CIF file requested. Exiting...")
    write_CIF(graph, cell)
    sys.exit()

sim.write_lammps_files(wd=run_dir)

with open(run_dir/f"in.box", 'a') as f:
    f.write("""
neighbor        2.0 bin
neigh_modify    delay 5
minimize 1.0e-4 1.0e-6 100000 100000

velocity all create 2273.0 1814 dist gaussian
fix 3 all npt temp 2000.0 2000.0 100.0 iso 0.0 0.0 1000
thermo 1000

neighbor        6.0 bin
neigh_modify    delay 0 every 1 check yes
timestep        0.2
run 5000

neighbor        2.0 bin
neigh_modify    delay 5
timestep        1.0

dump test1 all custom 10000 movie.data id type x y z
dump test all xyz 10000 movie.xyz
dump_modify test element Mn Fe Co Ni Zn N C Co
run            1000000
    """)

import subprocess
lammps_cmd = "mpirun -np 64 lmp_cpu -in in.box"

print(f"Running LAMMPS in: {run_dir}")
result = subprocess.run(
    lammps_cmd,
    cwd=str(run_dir),
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

print("=== LAMMPS stdout ===")
print(result.stdout)

print("=== LAMMPS stderr ===")
print(result.stderr)
