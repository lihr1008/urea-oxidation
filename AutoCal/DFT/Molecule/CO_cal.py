from atomate2.vasp.jobs.adsorption import MolRelaxMaker,MolStaticMaker
from atomate2.vasp.powerups import update_user_kpoints_settings,update_user_incar_settings
from jobflow import run_locally, JobStore
from maggma.stores.mongolike import MemoryStore
from monty.serialization import dumpfn
import pandas as pd
import copy
import os
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.inputs import KpointsSupportedModes
from jobflow import Flow


monkhorst_kpoints = Kpoints(
    comment="Automatic mesh",
    num_kpts=0,
    style=KpointsSupportedModes.Monkhorst,
    kpts=[(1, 1, 1)],
    kpts_shift=(0, 0, 0),
)
energies={
    'relax_energy':[],
    'static_energy':[]
}


parser = CifParser(f'CO.cif')
structure = parser.parse_structures(primitive=False)

# print(structure)
#generate the flow and reduce the k-point mesh for the relaxation jobs
relax_maker = MolRelaxMaker()
relax_job = relax_maker.make(structure[0])

magmom_dict = {element: 1.0 for element in structure[0].chemical_system_set}
incar_relax_dict = {
    # Control
    "SYSTEM": "CO",
    "ISTART": 0,
    "ICHARG": 2,
    "LCHARG": False,
    "LWAVE": False,
    "LREAL": "Auto",
    "NPAR": 4,
    "LASPH": False,

    # SCF
    "GGA": "PE",
    "ALGO": "Normal",
    "EDIFF": 1E-5,
    "ISMEAR": 0,
    "SIGMA": 0.1,
    "PREC": "Accurate",
    "NELM": 60,
    "NELMIN": 6,
    "ENCUT": 400,
    "IVDW": 12,

    # OPT
    "EDIFFG": -0.01,
    "NSW": 700,
    "IBRION": 2,
    "ISIF": 2,
    "LMAXMIX": 2,

    # Spin-polar
    "ISPIN": 2,
    "MAGMOM": magmom_dict

}

incar_static_dict = {
    # Control
    "SYSTEM": "CO",
    "ISTART": 0,
    "ICHARG": 2,
    "LCHARG": True,
    "LWAVE": False,
    "LREAL": "Auto",
    "NPAR": 4,
    "LASPH": False,

    # SCF
    "GGA": "PE",
    "ALGO": "Normal",
    "EDIFF": 1E-5,
    "ISMEAR": 0,
    "SIGMA": 0.1,
    "PREC": "Accurate",
    "NELM": 300,
    "NELMIN": 6,
    "ENCUT": 400,
    "IVDW": 12,

    # OPT
    "EDIFFG": -0.01,
    "NSW": 0,
    "IBRION": -1,
    "ISIF": 2,
    "LMAXMIX": 2,

    # bader
    "LAECHG":True,

    # Spin-polar
    "ISPIN": 2,
    "MAGMOM": magmom_dict

}

relax_job = update_user_incar_settings(relax_job,incar_relax_dict)
relax_job = update_user_kpoints_settings(
    flow=relax_job,
    kpoints_updates=monkhorst_kpoints,
    name_filter="relax"
)

static_maker = MolStaticMaker()
static_job = static_maker.make(relax_job.output.structure)
static_job = update_user_incar_settings(static_job, incar_static_dict)
static_job = update_user_kpoints_settings(
    flow=static_job,
    kpoints_updates=monkhorst_kpoints,
    name_filter="static"
)

jobs = [relax_job, static_job]
flow = Flow(jobs=jobs)


# run the workflow using a custom store so that we can easily compile test data
store = JobStore(MemoryStore(), additional_stores={"data": MemoryStore()})
dir_path = f'cal/CO'
os.makedirs(dir_path, exist_ok=True)
run_locally(
    flow,
    store=store,
    create_folders=True,
    ensure_success=True,
    root_dir=dir_path
)

# dump all of the job outputs to the outputs.json file in the current directory
outputs = list(store.query(load=True))
dumpfn(outputs,f"outputs_CO.json")
energies["relax_energy"].append(outputs[0]['output']['output']['energy'])
energies["static_energy"].append(outputs[1]['output']['output']['energy'])

df = pd.DataFrame(energies)
df.to_excel("energies.xlsx", index=False)
