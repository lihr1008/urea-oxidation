from atomate2.vasp.jobs.adsorption import SlabRelaxMaker,SlabStaticMaker
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
def apply_selective_constraints(struct,
                                fix_elements=None,
                                free_z_threshold=None,
                                fixed_indices=None):

    for site in struct:
        site.properties["selective_dynamics"] = [False, False, False]


    if fix_elements:
        for i, site in enumerate(struct):
            if site.specie.symbol in fix_elements:
                struct[i].properties["selective_dynamics"] = [True, True, True]


    if free_z_threshold:
        z_min, z_max = free_z_threshold
        for i, site in enumerate(struct):
            if z_min <= site.frac_coords[2] <= z_max:
                struct[i].properties["selective_dynamics"] = [True, True, True]
            else:
                struct[i].properties["selective_dynamics"] = [False, False, False]


    if fixed_indices:
        for idx in fixed_indices:
            struct[idx].properties["selective_dynamics"] =  [True, True, True]

    return struct
excel = r'metal_sites.xlsx'
data=pd.read_excel(excel)
metals=data.values
energies={
    'metals':[],
    'relax_energy':[],
    'static_energy':[]
}

monkhorst_kpoints = Kpoints(
    comment="Automatic mesh",
    num_kpts=0,
    style=KpointsSupportedModes.Monkhorst,
    kpts=[(1, 1, 1)],
    kpts_shift=(0, 0, 0),
)

# 读取POSCAR_urea文件
with open('slab_CO_temp.cif', 'r') as file:
    content = file.read()

new_content=copy.deepcopy(content)

for i in range(len(metals)):
    new_content = copy.deepcopy(content)
    modified_content = new_content.replace("Aa", metals[i][0])
    modified_content = modified_content.replace("Bb", metals[i][1])
    modified_content = modified_content.replace("Cc", metals[i][2])
    modified_content = modified_content.replace("Dd", metals[i][3])
    modified_content = modified_content.replace("Ee", metals[i][4])
    modified_content = modified_content.replace("Ff", metals[i][5])

    with open(f'cal/slab_CO_{metals[i][0]}_{metals[i][1]}_{metals[i][2]}_{metals[i][3]}_{metals[i][4]}_{metals[i][5]}.cif', 'w') as file:
        file.write(modified_content)


    parser = CifParser(f'cal/slab_CO_{metals[i][0]}_{metals[i][1]}_{metals[i][2]}_{metals[i][3]}_{metals[i][4]}_{metals[i][5]}.cif')
    structure = parser.parse_structures(primitive=False)
    modified_structure = apply_selective_constraints(
        structure[0],
        free_z_threshold=(0.3, 0.8),
    )
    #generate the flow and reduce the k-point mesh for the relaxation jobs
    relax_maker = SlabRelaxMaker()
    relax_job = relax_maker.make(modified_structure)

    magmom_dict = {element: 1.0 for element in modified_structure.chemical_system_set}
    incar_relax_dict = {
        # Control
        "SYSTEM": "slab_CO",
        "ISTART": 0,
        "ICHARG": 2,
        "LCHARG": False,
        "LWAVE": False,
        "LREAL": "Auto",
        "NPAR": 4,
        "LASPH": False,

        # SCF
        "GGA": "PE",
        "ALGO": "Fast",
        "EDIFF": 1E-4,
        "ISMEAR": 0,
        "SIGMA": 0.1,
        "PREC": "Normal",
        "NELM": 60,
        "NELMIN": 6,
        "ENCUT": 400,
        "IVDW": 12,

        # OPT
        "EDIFFG": -0.03,
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
        "SYSTEM": "slab_CO",
        "ISTART": 0,
        "ICHARG": 2,
        "LCHARG": True,
        "LWAVE": False,
        "LREAL": "Auto",
        "NPAR": 4,
        "LASPH": False,

        # SCF
        "GGA": "PE",
        "ALGO": "Fast",
        "EDIFF": 1E-4,
        "ISMEAR": 0,
        "SIGMA": 0.1,
        "PREC": "Normal",
        "NELM": 300,
        "NELMIN": 6,
        "ENCUT": 400,
        "IVDW": 12,

        # OPT
        "EDIFFG": -0.03,
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

    static_maker = SlabStaticMaker()
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
    dir_path = f'cal/{metals[i][0]}_{metals[i][1]}_{metals[i][2]}_{metals[i][3]}_{metals[i][4]}_{metals[i][5]}'
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
    dumpfn(outputs,
           f"outputs_{metals[i][0]}_{metals[i][1]}_{metals[i][2]}_{metals[i][3]}_{metals[i][4]}_{metals[i][5]}.json")
    energies["metals"].append(f'{metals[i][0]}_{metals[i][1]}_{metals[i][2]}_{metals[i][3]}_{metals[i][4]}_{metals[i][5]}')
    energies["relax_energy"].append(outputs[0]['output']['output']['energy'])
    energies["static_energy"].append(outputs[1]['output']['output']['energy'])

df = pd.DataFrame(energies)
df.to_excel("energies.xlsx", index=False)
