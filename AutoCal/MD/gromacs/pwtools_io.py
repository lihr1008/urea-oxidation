"""High level Structure and Trajectory I/O. HDF5 convenience IO functions."""

import warnings, os
try:
    import h5py
except ImportError:
    pass

import pickle
import numpy as np
from pwtools.common import frepr, cpickle_load
from pwtools.constants import Ha, eV
from pwtools import parse, atomic_data, lammps
from pwtools import crys
from pwtools import common
from pwtools import pwscf
##warnings.simplefilter('always')

try:
    import CifFile as pycifrw_CifFile
except ImportError:
    pass


def wien_sgroup_input(struct, lat_symbol='P'):
    """Generate input for WIEN2K's ``sgroup`` symmetry analysis tool.

    length: can be any

    Parameters
    ----------
    struct : Structure instance
    lat_symbol : str, e.g. 'P'

    Notes
    -----
    From sgroup's README:

    / ------------------------------------------------------------
    / in input file symbol "/" means a comment
    / and trailing characters are ignored by the program

    / empty lines are allowed

    P  /type of lattice; choices are P,F,I,C,A

    /  parameters of cell:
    /  lengths of the basis vectors and
    /  angles (degree unit is used)  alpha=b^c  beta=a^c  gamma=a^b
    /   |a|  |b|   |c|               alpha  beta  gamma

       1.0   1.1   1.2                90.   91.    92.

    /Number of atoms in the cell
    4

    /List of atoms
    0.1 0.2 0.3  / <-- Atom positions in units of the vectors a b c
    Al           / <-- name of this atom

    0.1 0.2 0.4  /....
    Al1

    0.2 0.2 0.3
    Fe

    0.1 0.3 0.3
    Fe

    / ------------------------------------------------------------------
    """
    empty = '\n\n'
    txt = "/ lattice type symbol\n%s" %lat_symbol
    txt += empty
    txt += "/ a b c alpha beta gamma\n"
    txt += " ".join(["%.13e"]*6) % tuple(struct.cryst_const)
    txt += empty
    txt += "/ number of atoms\n%i" %struct.natoms
    txt += empty
    txt += "/ atom list (crystal cooords)\n"
    fmt = ' '.join(['%.13e']*3)
    for sym, coord in zip(struct.symbols, struct.coords_frac):
        txt += fmt % tuple(coord) + '\n' + sym + '\n'
    return txt


def write_wien_sgroup(filename, struct, **kwds):
    """
    Write `struct` to input file for WIEN2K's ``sgroup`` symmetry analysis
    tool.

    Parameters
    ----------
    filename : str
        name of the output file
    struct : Structure
    **kwds : see wien_sgroup_input()
    """
    txt = wien_sgroup_input(struct, **kwds)
    common.file_write(filename, txt)


def write_cif(filename, struct):
    """Q'n'D Cif writer. Uses PyCifRW.

    length: Angstrom

    Parameters
    ----------
    filename : str
        name of output .cif file
    struct : Structure, length units Angstrom assumed
    """
    elements = list({}.fromkeys(struct.symbols).keys())
    with open("element.txt", 'w') as f:
      for s in elements:
       f.write(s + ' ')

    ffmt = "%.13e"
    cf = pycifrw_CifFile.CifFile()
    block = pycifrw_CifFile.CifBlock()

    block['_cell_length_a'] = frepr(struct.cryst_const[0], ffmt=ffmt)
    block['_cell_length_b'] = frepr(struct.cryst_const[1], ffmt=ffmt)
    block['_cell_length_c'] = frepr(struct.cryst_const[2], ffmt=ffmt)
    block['_cell_angle_alpha'] = frepr(struct.cryst_const[3], ffmt=ffmt)
    block['_cell_angle_beta'] = frepr(struct.cryst_const[4], ffmt=ffmt)
    block['_cell_angle_gamma'] = frepr(struct.cryst_const[5], ffmt=ffmt)
    block['_symmetry_space_group_name_H-M'] = 'P 1'
    block['_symmetry_Int_Tables_number'] = 1
    # assigning a list produces a "loop_"
    block['_symmetry_equiv_pos_as_xyz'] = ['x,y,z']

    # atoms
    #
    # _atom_site_label: We just use symbols, which is then =
    #   _atom_site_type_symbol, but we *could* use that to number atoms of each
    #   specie, e.g. Si1, Si2, ..., Al1, Al2, ...
    data_names = ['_atom_site_label',
                  '_atom_site_fract_x',
                  '_atom_site_fract_y',
                  '_atom_site_fract_z',
                  '_atom_site_type_symbol',
                  '_atom_site_charge']
    _xyz2str = lambda arr: [ffmt %x for x in arr]
    labels = []
    for i in range(struct.natoms):
      labels.append(struct.symbols[i]+str(i+1))
    atom_to_charge1 = {25:'2.000',
                       26:'2.000',
                       27:'2.000',
                       28:'2.000',
                       30:'2.000'}
    atom_to_charge2 = {27:'3.000',
                        6:'-0.510',
                        7:'-0.490'}

    charge = []
    for i in range(60):
      charge.append(atom_to_charge1[struct.znucl[i]])
    for i in range(520):
      charge.append(atom_to_charge2[struct.znucl[i+60]])

    data = [labels,
            _xyz2str(struct.coords_frac[:,0]),
            _xyz2str(struct.coords_frac[:,1]),
            _xyz2str(struct.coords_frac[:,2]),
            struct.symbols,
            charge]
    bond_name = ['_geom_bond_atom_site_label_1',
                 '_geom_bond_atom_site_label_2',
                 '_geom_bond_distance',
                 '_ccdc_geom_bond_type']
    
    # FeCoNiZnCrMn-OH
    # bond = [labels[60:struct.natoms:2],
    #         labels[60+1:struct.natoms:2],
    #         [str(1.0)]*int((struct.natoms-60)/2),
    #         ['S']*int((struct.natoms-60)/2)]

    # MgZnNiCuCo-TPA
    bond_1 = []
    bond_2 = []
    bond_3 = []
    bond_4 = []
    for i in range(40):
      bond_1.append(labels[59+13*i+1])
      bond_1.append(labels[59+13*i+2])
      bond_1.append(labels[59+13*i+3])
      bond_1.append(labels[59+13*i+4])
      bond_1.append(labels[59+13*i+5])
      bond_1.append(labels[59+13*i+6])
      bond_1.append(labels[59+13*i+7])
      bond_1.append(labels[59+13*i+8])
      bond_1.append(labels[59+13*i+9])
      bond_1.append(labels[59+13*i+10])
      bond_1.append(labels[59+13*i+11])
      bond_1.append(labels[59+13*i+12])
      bond_2.append(labels[59+13*i+7])
      bond_2.append(labels[59+13*i+9])
      bond_2.append(labels[59+13*i+10])
      bond_2.append(labels[59+13*i+8])
      bond_2.append(labels[59+13*i+11])
      bond_2.append(labels[59+13*i+12])
      bond_2.append(labels[59+13*i+13])
      bond_2.append(labels[59+13*i+13])
      bond_2.append(labels[59+13*i+13])
      bond_2.append(labels[59+13*i+13])
      bond_2.append(labels[59+13*i+13])
      bond_2.append(labels[59+13*i+13])
      bond_3.append('1.173')
      bond_3.append('1.173')
      bond_3.append('1.173')
      bond_3.append('1.173')
      bond_3.append('1.173')
      bond_3.append('1.173')
      bond_3.append('1.923')
      bond_3.append('1.923')
      bond_3.append('1.923')
      bond_3.append('1.923')
      bond_3.append('1.923')
      bond_3.append('1.923')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
      bond_4.append('S')
    bond = [bond_1,bond_2,bond_3,bond_4]

    # "loop_" with multiple columns
    block.AddCifItem([[data_names], [data]])
    block.AddCifItem([[bond_name], [bond]])
    cf['pwtools'] = block
    # maxoutlength = 2048 is default for cif 1.1 standard (which is default in
    # pycifrw 3.x). Reset default wraplength=80 b/c ASE's cif reader cannot
    # handle wrapped lines.
    common.file_write(filename, cf.WriteOut(wraplength=2048))


def write_xyz(filename, obj, name='pwtools_dummy_mol_name'):
    """Write VMD-style [VMD] XYZ file.

    length: Angstrom

    Parameters
    ----------
    filename : target file name
    obj : Trajectory or Structure
    name : str, optional
        Molecule name.

    References
    ----------
    [VMD] http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html
    """
    traj = crys.struct2traj(obj)
    xyz_str = ""
    for istep in range(traj.nstep):
        xyz_str += "%i\n%s\n%s" %(traj.natoms,
                                  name + '.%i' %(istep + 1),
                                  pwscf.atpos_str_fast(traj.symbols,
                                                       traj.coords[istep,...]),
                                  )
    common.file_write(filename, xyz_str)


def write_axsf(filename, obj):
    """Write animated XSF file for Structure (only 1 step) or Trajectory.

    Note that forces are converted eV / Ang -> Ha / Ang.

    length: Angstrom
    forces: Ha / Angstrom

    Parameters
    ----------
    filename : target file name
    obj : Structure or Trajectory

    References
    ----------
    [XSF] http://www.xcrysden.org/doc/XSF.html
    """
    # Notes
    # -----
    # XSF: The XSF spec [XSF] is a little fuzzy about what PRIMCOORD actually
    #     is (fractional or cartesian Angstrom). Only the latter case results
    #     in a correctly displayed structure in xcrsyden. So we use that.
    #
    # Speed: The only time-consuming step is calling atpos_str*() in the loop
    #     b/c that transforms *every* single float to a string, which
    #     effectively is a double loop over `ccf`. No way to get faster w/ pure
    #     Python.
    #
    traj = crys.struct2traj(obj)
    # ccf = cartesian coords + forces (6 columns)
    if traj.is_set_attr('forces'):
        ccf = np.concatenate((traj.coords, traj.forces*eV/Ha), axis=-1)
    else:
        ccf = traj.coords
    axsf_str = "ANIMSTEPS %i\nCRYSTAL" %traj.nstep
    for istep in range(traj.nstep):
        axsf_str += "\nPRIMVEC %i\n%s" %(istep+1,
                                         common.str_arr(traj.cell[istep,...]))
        axsf_str += "\nPRIMCOORD %i\n%i 1\n%s" %(istep+1,
                                                 traj.natoms,
                                                 pwscf.atpos_str_fast(traj.symbols,
                                                                      ccf[istep,...]))
    common.file_write(filename, axsf_str)


def write_lammps(filename, struct, symbolsbasename='lmp.struct.symbols'):
    """Write Structure object to lammps format. That file can be read in a
    lammps input file by ``read_data``. Write file ``lmp.struct.symbols`` with
    atom symbols.

    Parameters
    ----------
    filename : str
        name of file to write
    symbolsbasename : str, optional
        file for atom symbols
    struct : Structure

    References
    ----------
    ase.calculators.lammpsrun (ASE 3.8).
    """
    dr = os.path.dirname(filename)
    fn = os.path.join(dr, symbolsbasename)
    common.file_write(fn, '\n'.join(struct.symbols))
    common.file_write(filename, lammps.struct_str(struct))


def write_h5(fn, dct, **kwds):
    """Write dictionary with arrays (or whatever HDF5 handles) to h5 file.

    Dict keys are supposed to be HDF group + dataset names like `/a/b/c/dset`.
    The leading slash can be skipped. Default file mode is 'a' (see below for
    details).

    Parameters
    ----------
    fn : str
        filename (e.g. 'foo.h5', 'bar.hdf')
    dct : dict
    **kwds :
        keywords to ``h5py.File`` (e.g. ``mode='w'``)

    Notes
    -----
    The default file opening mode is the (old?) ``h5py.File`` default value,
    which is ``mode='a'``, i.e. read+append mode. In this mode, existing keys
    cannot be reused (overwritten), only new ones can be appended. The file is
    created if nonexistent. To overwrite, use ``mode='w'``, but this is the
    same as deleting the file and writing a new one! If you want to overwrite
    some or all existing keys and add new ones, use smth like::

    >>> old = read_h5('file.h5')
    >>> old.update({'/old/key': new_value, '/new/key': some_more_data})
    >>> write_h5('file.h5', old, mode='w')
    """
    mode = kwds.pop('mode') if 'mode' in kwds else 'a'
    fh = h5py.File(fn, mode=mode, **kwds)
    for key,val in dct.items():
        fh[key] = val
    fh.close()


def read_h5(fn):
    """Read h5 file into dict.

    Dict keys are the group + dataset names, e.g. '/a/b/c/dset'. All keys start
    with a leading slash even if written without (see :func:`write_h5`).

    Parameters
    ----------
    fn : str
        filename

    Examples
    --------
    >>> read_h5('foo.h5').keys()
    ['/a/b/d1', '/a/b/d2', '/a/c/d3', '/x/y/z']
    """
    fh = h5py.File(fn, mode='r')
    dct = {}
    def get(name, obj):
        if isinstance(obj, h5py.Dataset):
            _name = name if name.startswith('/') else '/'+name
            val = obj[()]
            dct[_name] = obj.asstr()[()] if isinstance(val, bytes) else val
    fh.visititems(get)
    fh.close()
    return dct


def load_h5(*args, **kwds):
    """Alias for :func:`read_h5`. Deprecated."""
    warnings.warn("load_h5() is deprcated, use read_h5() instead",
                   DeprecationWarning)
    return read_h5(*args, **kwds)


def read_pickle(filename):
    """Load object written by ``pickle.dump()``, e.g. files written by
    :meth:`~pwtools.base.FlexibleGetters.dump()`."""
    return pickle.load(open(filename, 'rb'))


class ReadFactory:
    """Factory class to construct callables to parse files."""
    def __init__(self, parser=None, struct_or_traj=None, doc=''):
        """
        Parameters
        ----------
        parser : one of parse.*File parsing classes
        struct_or_traj : str
            {'struct','traj'}
            Whether the callables should return parser.get_{struct,traj}()'s
            return value.
        doc : str
            docstring header
        """
        self.parser = parser
        self.struct_or_traj = struct_or_traj
        # Overwrite class docstring. That shows up in sphinx autodoc of the
        # class instances (all read_* "functions"). The more natural thing is
        # to tell sphinx to use the __call__.__doc__ since the instance is used
        # as a callable. Now, they are treated as normal class instance, which
        # is perfectly right from sphinx' point of view. We would need to add
        # `doc` to __call__.__doc__, but how? Fancy decorator magic!
        self.__doc__ = doc + """

        Parameters
        ----------
        filename : str
            Name of the file to parse.
        **kwds : keywords args
            passed to the parser class (e.g. units=...)

        Returns
        -------
        ret : :class:`~pwtools.crys.Structure` (SCF runs) or \
              :class:`~pwtools.crys.Trajectory` (MD-like runs)
        """

    def __call__(self, filename, **kwds):
        """
        Parameters
        ----------
        filename : str
            Name of the file to parse.
        **kwds : keywords args
            passed to the parser class (e.g. units=...)
        """
        if self.struct_or_traj == 'struct':
            return self.parser(filename, **kwds).get_struct()
        elif self.struct_or_traj == 'traj':
            return self.parser(filename, **kwds).get_traj()
        else:
            raise Exception("unknown struct_or_traj: %s" %struct_or_traj)


read_cif = ReadFactory(parser=parse.CifFile,
                       struct_or_traj='struct',
                       doc="Read Cif files."
                       )
read_pdb = ReadFactory(parser=parse.PDBFile,
                       struct_or_traj='struct',
                       doc="Read PDB files."
                       )
read_pw_scf = ReadFactory(parser=parse.PwSCFOutputFile,
                          struct_or_traj='struct',
                          doc="Read Pwscf SCF run ouput."
                          )
read_pw_md = ReadFactory(parser=parse.PwMDOutputFile,
                         struct_or_traj='traj',
                         doc="Read Pwscf md/relax/vc-relax run ouput."
                         )
read_pw_vcmd = ReadFactory(parser=parse.PwVCMDOutputFile,
                           struct_or_traj='traj',
                           doc="Read Pwscf vc-md run ouput."
                           )
read_cpmd_scf = ReadFactory(parser=parse.CpmdSCFOutputFile,
                            struct_or_traj='struct',
                            doc="Read CPMD SCF (single point) run ouput."
                            )
read_cpmd_md = ReadFactory(parser=parse.CpmdMDOutputFile,
                           struct_or_traj='traj',
                           doc="Read CPMD MD (fixed and variable cell, BO and CP) run ouput."
                           )
read_cp2k_scf = ReadFactory(parser=parse.Cp2kSCFOutputFile,
                           struct_or_traj='struct',
                           doc="Read CP2K SCF (single point) run ouput."
                           )
read_cp2k_md = ReadFactory(parser=parse.Cp2kMDOutputFile,
                           struct_or_traj='traj',
                           doc="Read CP2K MD run ouput (all text)."
                           )
read_cp2k_md_dcd = ReadFactory(parser=parse.Cp2kDcdMDOutputFile,
                               struct_or_traj='traj',
                               doc="Read CP2K MD run ouput (coordinates in dcd binary format)."
                               )
read_cp2k_relax = ReadFactory(parser=parse.Cp2kRelaxOutputFile,
                              struct_or_traj='traj',
                              doc="Read CP2K relaxation run ouput (all text)."
                              )
read_lammps_md_txt = ReadFactory(parser=parse.LammpsTextMDOutputFile,
                                 struct_or_traj='traj',
                                 doc="Read LAMMPS MD run ouput (all text)."
                                 )
read_lammps_md_dcd = ReadFactory(parser=parse.LammpsDcdMDOutputFile,
                                 struct_or_traj='traj',
                                 doc="Read LAMMPS MD run ouput (coordinates in dcd format)."
                                 )

