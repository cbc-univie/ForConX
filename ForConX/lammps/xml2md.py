import numpy as np
#import xml.etree.ElementTree as ET
from ..md_xml import helpers, input_output
from ..md_xml import molecule
from ..md_xml import bonds, angles, dihedrals, impropers
from ..md_xml import nonbonded

def read_pdb(pdbfile):
    atoms = []
    box = [0.,0.,0.]
    pdb = open(pdbfile, 'r')
    for line in pdb:
        if line.startswith('CRYST1'):
            #  COLUMNS       DATA  TYPE    FIELD          DEFINITION
            #  -------------------------------------------------------------
            #   1 -  6       Record name   "CRYST1"
            #   7 - 15       Real(9.3)     a              a (Angstroms).
            #  16 - 24       Real(9.3)     b              b (Angstroms).
            #  25 - 33       Real(9.3)     c              c (Angstroms).
            #  34 - 40       Real(7.2)     alpha          alpha (degrees).
            #  41 - 47       Real(7.2)     beta           beta (degrees).
            #  48 - 54       Real(7.2)     gamma          gamma (degrees).
            #  56 - 66       LString       sGroup         Space  group.
            #  67 - 70       Integer       z              Z value.
            a = float(line[7:16])
            b = float(line[16:25])
            c = float(line[25:34])
            alpha = np.radians(float(line[34:41]))
            beta = np.radians(float(line[41:48]))
            gamma = np.radians(float(line[48:55]))
            yz = c*(np.cos(alpha)-np.cos(gamma)*np.cos(beta))/np.sin(gamma)
            box = [a, b*np.sin(gamma), c*np.sqrt(np.sin(beta)**2 - yz**2)] # lz
            if float(line[34:41]) != 90. or float(line[41:48]) != 90. or float(line[48:55]) != 90.:
                box.append(b*np.cos(gamma)) # xy
                box.append(c*np.cos(beta)) # xz
                box.append(yz) # yz
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            #  COLUMNS        DATA  TYPE    FIELD        DEFINITION
            #  -------------------------------------------------------------------------------------
            #   1 -  6        Record name   "ATOM  "
            #   7 - 11        Integer       serial       Atom  serial number.
            #  13 - 16        Atom          name         Atom name.
            #  17             Character     altLoc       Alternate location indicator.
            #  18 - 20        Residue name  resName      Residue name.
            #  22             Character     chainID      Chain identifier.
            #  23 - 26        Integer       resSeq       Residue sequence number.
            #  27             AChar         iCode        Code for insertion of residues.
            #  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
            #  39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
            #  47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
            #  55 - 60        Real(6.2)     occupancy    Occupancy.
            #  61 - 66        Real(6.2)     tempFactor   Temperature  factor.
            #  77 - 78        LString(2)    element      Element symbol, right-justified.
            #  79 - 80        LString(2)    charge       Charge  on the atom.
            iatom = int(line[7:12])
            aname = line[12:17].strip()
            mname = line[17:22].strip().split()[0]
            imol = int(line[23:27])
            x = float(line[31:39])
            y = float(line[39:47])
            z = float(line[47:54])
            atoms.append((iatom, aname, imol, mname, x, y, z))
        else: continue
    return atoms, box

def write(root):
    inlmpelt = root.find('output/command')
    if inlmpelt is None:
        inlmp = 'in.lmp'
    else:
        inlmp = inlmpelt.get("file")
    datalmpelt = root.find('output/data')
    if datalmpelt is None:
        datalmp = 'data.lmp'
    else:
        datalmp = datalmpelt.get("file")
    molextelt = root.find('output/molecule')
    if molextelt is None:
        mol_suffix = '.mol'
    else:
        mol_suffix = molextelt.get("ext")
    pairlmpelt = root.find('output/pair')
    if pairlmpelt is None:
        pairlmp = 'pair_coeffs.lmp'
    else:
        pairlmp = pairlmpelt.get("file")
    nonbondedelt = nonbonded.nonbondedElement(root)
    mixing = nonbondedelt.mixing_sigma
    if nonbondedelt.mixing_epsilon != "geometric":
        nonbonded.atom2vdw()

    co = input_output.mdElement(root,'output')
    units = (co.convert_distance(output="ANGSTROEM"), co.convert_energy(output="KCAL"))

    out_formats = {}
    natoms = 0
    nbonds = 0
    nangles = 0
    ndihedrals = 0
    nimpropers = 0
    molecules = []
    for item in root.findall('./molecule'):
        mname = item.get('name')
        molecules.append(mname)
        mol = molecule.moleculeElement(item,mname)
        nmol = mol.nmol
        natoms += nmol * len(mol.list('atom'))
        nbonds += nmol * len(mol.list('bond'))
        nangles += nmol * len(mol.list('angle'))
        ndihedrals += nmol * len(mol.list('dihedral'))
        nimpropers += nmol * len(mol.list('improper'))

    atom_types = set()
    for item in root.findall('./molecule'):
        mname = item.get('name')
        mol = molecule.moleculeElement(item,mname)
        names = mol.list('atom') #molecule.namelist_atom(item)
        for name in names:
            obj = molecule.atomClass(item, name)
            tp = obj.type
            m = obj.mass
            q = obj.charge
            atom_types.add((tp, m, q))
    atom_types = list(atom_types)
    atom_style = ['charge', 'template mymols']
    natom_types = len(atom_types)

    pair_types = set()
    pair_style = set()
    nonb = nonbonded.nonbondedElement(root)
    nonb.atom2vdw()
    factor_coul = nonb.elec14_total()
    factor_lj = nonb.vdw14_total()
    for item in nonb.list('vdw'):
        style = "lj/charmm/coul/long"
        options = "14. 16."
        obj = nonbonded.vdwClass(root,item)
        tp = obj.type
        sigma = obj.sigma * units[0]
        epsilon = obj.epsilon * units[1]
        sigma14 = sigma * units[0]
        if obj.vdw14 is None:
            epsilon14 = epsilon
        else:
            epsilon14 = epsilon * obj.vdw14
        pair_types.add((style, tp, epsilon, sigma, epsilon14, sigma14))
        if ('pair', style) not in out_formats:
            out_formats['pair',style] = "%16.8f %16.8f %16.8f %16.8f"
        pair_style.add((style, options))
    pair_types = list(pair_types)
    pair_style = list(pair_style)
    npair_types = len(pair_types)

    bond_types = set()
    bond_style = set()
    shaked_bonds = set()
    for item in root.findall('./bonds/*'):
        tag = item.tag.lower()
        tp = item.get('type')
        style = tag.upper()
        if tag == 'harm':
            style = 'harmonic'
            obj = bonds.harmClass(root, tp)
            r0 = obj.r0 * units[0]
            try:
                k = float(obj.k) * units[1]/units[0]**2
            except TypeError:
                k = 750. # kcal/mol (arbitrary value for shaked bonds)
                shaked_bonds.add((style, tp, k, r0))
            bond_types.add((style, tp, k, r0))
            if ('bond', style) not in out_formats:
                out_formats['bond',style] = "%16.8f %16.8f"
        bond_style.add(style)
    bond_types = list(bond_types)
    bond_style = list(bond_style)
    nbond_types = len(bond_types)

    angle_types = set()
    angle_style = set()
    angle_map = {}
    for item in root.findall('./angles/*'): # same angle can have multiple entries if urey style
        style = item.tag.lower()
        angle_style.add(style)
        tp = item.get('type')
        if tp in angle_map: angle_map[tp].add(style)
        else: angle_map[tp] = set([style])
    for tp, styles in angle_map.items():
        if 'urey' in angle_style:
            k, theta0, kr, r0 = 0., 0., 0., 0.
            if 'harm' in styles:
                obj = angles.harmClass(root, tp)
                k = float(obj.k) * units[1]
                theta0 = obj.theta0
            if 'urey' in styles:
                obj = angles.ureyClass(root, tp)
                kr = float(obj.k) * units[1]/units[0]**2
                r0 = float(obj.r0) * units[0]
            angle_types.add(('charmm', tp, k, theta0, kr, r0))
        else:
            if styles == set(['harm']):
                obj = angles.harmClass(root, tp)
                k = float(obj.k) * units[1]
                theta0 = obj.theta0
                angle_types.add(('harmonic', tp, k, theta0))
    if 'urey' in angle_style:
        angle_style.discard('urey')
        angle_style.discard('harm')
        angle_style.add('charmm')
    if 'harm' in angle_style:
        angle_style.discard('harm')
        angle_style.add('harmonic')
    out_formats['angle','harmonic'] = "%16.8f %16.8f"
    out_formats['angle','charmm'] = "%16.8f %16.8f %16.8f %16.8f"
    angle_types = list(angle_types)
    angle_style = list(angle_style)
    nangle_types = len(angle_types)

    dihedral_types = set()
    dihedral_style = set()
    for item in root.findall('./molecule'):
        mname = item.get('name')
        mol = molecule.moleculeElement(item, mname)
        names = mol.list('dihedral') #molecule.namelist_dihedral(item)
        for name in names:
            obj = molecule.dihedralClass(item, name)
            tp = obj.type
            item2 = root.find('./dihedrals/*[@type="'+tp+'"]')
            tag = item2.tag.lower()
            style = tag.upper()
            if tag == 'cos':
                style = 'charmm'
                obj2 = dihedrals.cosClass(root, tp)
                ks = obj2.k
                ns = obj2.n
                deltas = obj2.delta
                for i in xrange(len(ns)):
                    k = ks[i] * units[1]
                    # if k == 0: continue
                    n = ns[i]
                    delta = deltas[i]
                    factor14 = 0.
                    dihedral_types.add((style, tp, k, int(n), int(delta), factor14))
                if ('dihedral', style) not in out_formats:
                    out_formats['dihedral',style] = "%16.8f %6d %6d %16.8f"
            elif tag == 'ryck':
                # WARNING style has to be charmm if pair style is lj/charmm/coul/...
                #style = 'multi/harmonic'
                #obj2 = dihedrals.ryckClass(root, tp)
                #dihedral_types.add((style, tp) + tuple([k * units[1] for k in obj2.k[:5]]))
                #if ('dihedral', style) not in out_formats:
                #    out_formats['dihedral',style] = "%16.8f %16.8f %16.8f %16.8f %16.8f"
                style = 'charmm'
                obj2 = dihedrals.ryckClass(root, tp)
                (ns, ks, deltas) = obj2.cos()
                kryck = obj2.k
                shift = kryck[0]-kryck[1]+kryck[2]-kryck[3]+kryck[4]-kryck[5]
                ns.append(0)
                ks.append(0.5*shift)
                deltas.append(0)
                for i in xrange(len(ns)):
                    k = ks[i] * units[1]
                    if k == 0. and (i > 0 or any([ki != 0. for ki in ks])): continue
                    n = ns[i]
                    delta = deltas[i]
                    factor14 = 0.
                    dihedral_types.add((style, tp, k, int(n), int(delta), factor14))
                if ('dihedral', style) not in out_formats:
                    out_formats['dihedral',style] = "%16.8f %6d %6d %16.8f"
            dihedral_style.add(style)
    dihedral_types = list(dihedral_types)
    dihedral_style = list(dihedral_style)
    ndihedral_types = len(dihedral_types)

    improper_types = set()
    improper_style = set()
    for item in root.findall('./molecule'):
        mname = item.get('name')
        mol = molecule.moleculeElement(item, mname)
        names = mol.list('improper')
        for name in names:
            obj = molecule.improperClass(item, name)
            tp = obj.type
            item2 = root.find('./impropers/*[@type="'+tp+'"]')
            tag = item2.tag.lower()
            style = tag.upper()
            if tag == 'cos' or tag == 'ryck':
                style = 'cvff'
                if tag == 'cos':
                    obj2 = impropers.cosClass(root, tp)
                    ks = obj2.k
                    ns = obj2.n
                    deltas = obj2.delta
                else:
                    obj2 = impropers.ryckClass(root, tp)
                    (ns, ks, deltas) = obj2.cos()
                    kryck = obj2.k
                    shift = kryck[0]-kryck[1]+kryck[2]-kryck[3]+kryck[4]-kryck[5]
                    ns.append(0)
                    ks.append(0.5*shift)
                    deltas.append(0)
                for i in xrange(len(ns)):
                    k = ks[i] * units[1]
                    if k == 0: continue
                    n = ns[i]
                    if deltas[i] == 0: delta = 1
                    elif deltas[i] == 180: delta = -1
                    else:
                        helpers.warning("LAMMPS cannot handle this improper. I skip.")
                        continue
                    improper_types.add((style, tp, k, int(n), int(delta)))
                if ('improper', style) not in out_formats:
                    out_formats['improper',style] = "%16.8f %6d %6d"
            elif tag == 'harm':
                style = 'harmonic'
                obj2 = impropers.harmClass(root, tp)
                k = obj2.k * units[1]
                theta0 = obj2.theta0
                improper_types.add((style, tp, k, theta0))
                if ('improper', style) not in out_formats:
                    out_formats['improper',style] = "%16.8f %16.8f"
            improper_style.add(style)
    improper_types = list(improper_types)
    improper_style = list(improper_style)
    nimproper_types = len(improper_types)

    pairtxt = "# Pair coeffs\n"
    for itp, tp1 in enumerate(atom_types):
        tp1 = tp1[0]
        for jtp, tp2 in enumerate(atom_types[itp:]):
            tp2 = tp2[0]
            #print tp1, tp2, type(tp1), ' '.join([tp1,tp2])
            ptp = bonds.sequence(' '.join([tp1,tp2]))
            for p in pair_types:
                if p[1] != ptp: continue
                pairtxt += "pair_coeff %4d %4d " % (itp+1, itp+jtp+1)
                if len(pair_style) > 1: pairtxt += p[0] + ' '
                pairtxt += out_formats['pair', p[0]] % p[2:]
                pairtxt += " # %s\n" % p[1]

    with open(datalmp, "w") as f:
        print >>f, "# LAMMPS data file generated by ForConX"
        print >>f
        if natom_types > 0: print >>f, "%6d atom types" % natom_types
        if nbond_types > 0: print >>f, "%6d bond types" % nbond_types
        if nangle_types > 0: print >>f, "%6d angle types" % nangle_types
        if ndihedral_types > 0: print >>f, "%6d dihedral types" % ndihedral_types
        if nimproper_types > 0: print >>f, "%6d improper types" % nimproper_types
        if natoms > 0:
            atoms, box = read_pdb("forconx.pdb")
            print >>f, "%6d atoms" % natoms
            print >>f
            print >>f, "0. %10.6f xlo xhi" % box[0]
            print >>f, "0. %10.6f ylo yhi" % box[1]
            print >>f, "0. %10.6f zlo zhi" % box[2]
            if len(box) > 3:
                print >>f, "%10.6f %10.6f %10.6f xy xz yz”" % tuple(box[3:6])

        #if nbonds > 0: print >>f, "%6d bonds" % nbonds
        #if nangles > 0: print >>f, "%6d angles" % nangles
        #if ndihedrals > 0: print >>f, "%6d dihedrals" % ndihedrals
        #if nimpropers > 0: print >>f, "%6d impropers" % nimpropers
        if natom_types > 0:
            print >>f
            print >>f, "Masses"
            print >>f
            for itp, tp in enumerate(atom_types):
                print >>f, "%6d %10.6f # %s" % (itp+1, tp[1], tp[0])
        if nbond_types > 0:
            print >>f
            print >>f, "Bond Coeffs #",
            if len(bond_style) > 1: print >>f, "hybrid",
            for style in bond_style: print >>f, style,
            print >>f
            print >>f
            for itp, tp in enumerate(bond_types):
                print >>f, "%6d" % (itp+1),
                if len(bond_style) > 1: print >>f, tp[0],
                print >>f, out_formats['bond',tp[0]] % tp[2:], "# %s" % tp[1]
        if nangle_types > 0:
            print >>f
            print >>f, "Angle Coeffs #",
            if len(angle_style) > 1: print >>f, "hybrid",
            for style in angle_style: print >>f, style,
            print >>f
            print >>f
            for itp, tp in enumerate(angle_types):
                print >>f, "%6d" % (itp+1),
                if len(angle_style) > 1: print >>f, tp[0],
                print >>f, out_formats['angle',tp[0]] % tp[2:], "# %s" % tp[1]
        if ndihedral_types > 0:
            print >>f
            print >>f, "Dihedral Coeffs #",
            if len(dihedral_style) > 1: print >>f, "hybrid",
            for style in dihedral_style: print >>f, style,
            print >>f
            print >>f
            for itp, tp in enumerate(dihedral_types):
                print >>f, "%6d" % (itp+1),
                if len(dihedral_style) > 1: print >>f, tp[0],
                print >>f, out_formats['dihedral',tp[0]] % tp[2:], "# %s" % tp[1]
        if nimproper_types > 0:
            print >>f
            print >>f, "Improper Coeffs #",
            if len(improper_style) > 1: print >>f, "hybrid",
            for style in improper_style: print >>f, style,
            print >>f
            print >>f
            for itp, tp in enumerate(improper_types):
                print >>f, "%6d" % (itp+1),
                if len(improper_style) > 1: print >>f, tp[0],
                print >>f, out_formats['improper',tp[0]] % tp[2:], "# %s" % tp[1]
        if natoms > 0:
            print >>f
            print >>f, "Atoms"
            print >>f
            imols = {}
            for atom in atoms:
                #item = root.find('./molecule[@name="'+atom[3]+'"]')
                items = root.findall('./molecule')
                for item in items:
                    if item.get('name').upper().startswith(atom[3].upper()): break
                obj = molecule.atomClass(item, atom[1])
                tp = obj.type
                m = obj.mass
                q = obj.charge
                itp = atom_types.index((tp, m, q))
                print >>f, "%8d %6d %16.8f %16.8f %16.8f %16.8f " % (atom[0], itp+1, atom[4], atom[5], atom[6], q),
                imtp = root.findall('./molecule').index(item)
                item2 = item.find('./atom[@name="'+atom[1]+'"]')
                iinm = item.findall('./atom').index(item2)
                if (atom[2], imtp) in imols:
                    imol = imols[atom[2], imtp]
                else:
                    imol = len(imols)
                    imols[atom[2], imtp] = imol
                print >>f, "%8d %4d %6d # %s (%s, %s)" % (imol+1, imtp+1, iinm+1, atom[1], atom[3], tp)
            del imols

    for item in root.findall('./molecule'):
        mname = item.get('name')
        mol = molecule.moleculeElement(item, mname)
        atomnames = mol.list('atom')
        bondnames = mol.list('bond')
        anglenames = mol.list('angle')
        dihedralnames = mol.list('dihedral')
        ndihedralnames = 0
        for i, name in enumerate(dihedralnames):
            obj = molecule.dihedralClass(item, name)
            tp = obj.type
            for itp, dihedraltype in enumerate(dihedral_types):
                if dihedraltype[1] == tp: ndihedralnames += 1
        impropernames = mol.list('improper')
        nimpropernames = 0
        for i, name in enumerate(impropernames):
            obj = molecule.improperClass(item, name)
            tp = obj.type
            for itp, impropertype in enumerate(improper_types):
                if impropertype[1] == tp: nimpropernames += 1
        with open(mname.lower()+mol_suffix, "w") as f:
            print >>f, "# LAMMPS molecule template file for %s generated by ForConX" % mname
            print >>f
            if len(atomnames) > 0: print >>f, "%8d atoms" % len(atomnames)
            if len(bondnames) > 0: print >>f, "%8d bonds" % len(bondnames)
            if len(anglenames) > 0: print >>f, "%8d angles" % len(anglenames)
            if len(dihedralnames) > 0: print >>f, "%8d dihedrals" % ndihedralnames
            if len(impropernames) > 0: print >>f, "%8d impropers" % nimpropernames
            if len(atomnames) > 0:
                print >>f
                print >>f, "Types"
                print >>f
                for i, name in enumerate(atomnames):
                    obj = molecule.atomClass(item, name)
                    tp = obj.type
                    m = obj.mass
                    q = obj.charge
                    itp = atom_types.index((tp, m, q))
                    #print >>f, "%8d %6d # %s" % (i+1, itp+1, tp)
                    print >>f, "%8d %6d" % (i+1, itp+1)
                #print >>f
                #print >>f, "Masses"
                #print >>f
                #for i, name in enumerate(atomnames):
                #    obj = molecule.atomClass(item, name)
                #    tp = obj.type
                #    m = obj.mass
                #    #print >>f, "%8d %10.6f # %s" % (i+1, m, tp)
                #    print >>f, "%8d %10.6f" % (i+1, m)
                print >>f
                print >>f, "Charges"
                print >>f
                for i, name in enumerate(atomnames):
                    obj = molecule.atomClass(item, name)
                    tp = obj.type
                    q = obj.charge
                    #print >>f, "%8d %10.6f # %s" % (i+1, q, tp)
                    print >>f, "%8d %10.6f" % (i+1, q)
            if len(bondnames) > 0:
                print >>f
                print >>f, "Bonds"
                print >>f
                for i, name in enumerate(bondnames):
                    obj = molecule.bondClass(item, name)
                    tp = obj.type
                    for itp, bondtype in enumerate(bond_types):
                        if bondtype[1] == tp: break
                    i1, i2 = map(atomnames.index, name.split())
                    #print >>f, "%8d %6d %8d %8d # %s" % (i+1, itp+1, i1+1, i2+1, tp)
                    print >>f, "%8d %6d %8d %8d" % (i+1, itp+1, i1+1, i2+1)
            if len(anglenames) > 0:
                print >>f
                print >>f, "Angles"
                print >>f
                for i, name in enumerate(anglenames):
                    obj = molecule.angleClass(item, name)
                    tp = obj.type
                    for itp, angletype in enumerate(angle_types):
                        if angletype[1] == tp: break
                    i1, i2, i3 = map(atomnames.index, name.split())
                    #print >>f, "%8d %6d %8d %8d %8d # %s" % (i+1, itp+1, i1+1, i2+1, i3+1, tp)
                    print >>f, "%8d %6d %8d %8d %8d" % (i+1, itp+1, i1+1, i2+1, i3+1)
            if len(dihedralnames) > 0:
                print >>f
                print >>f, "Dihedrals"
                print >>f
                i = 0
                for name in dihedralnames:
                    obj = molecule.dihedralClass(item, name)
                    tp = obj.type
                    for itp, dihedraltype in enumerate(dihedral_types):
                        if dihedraltype[1] == tp:
                            i1, i2, i3, i4 = map(atomnames.index, name.split())
                            #print >>f, "%8d %6d %8d %8d %8d %8d # %s" % (i+1, itp+1, i1+1, i2+1, i3+1, i4+1, tp)
                            print >>f, "%8d %6d %8d %8d %8d %8d" % (i+1, itp+1, i1+1, i2+1, i3+1, i4+1)
                            i += 1
            if len(impropernames) > 0:
                print >>f
                print >>f, "Impropers"
                print >>f
                i = 0
                for name in impropernames:
                    obj = molecule.improperClass(item, name)
                    tp = obj.type
                    for itp, impropertype in enumerate(improper_types):
                        if impropertype[1] == tp:
                            i1, i2, i3, i4 = map(atomnames.index, name.split())
                            #print >>f, "%8d %6d %8d %8d %8d %8d # %s" % (i+1, itp+1, i1+1, i2+1, i3+1, i4+1, tp)
                            print >>f, "%8d %6d %8d %8d %8d %8d" % (i+1, itp+1, i1+1, i2+1, i3+1, i4+1)
                            i += 1
            if len(atomnames) > 0:
                print >>f
                print >>f, "Special Bond Counts"
                print >>f
                neigh = {}
                for i, name in enumerate(atomnames):
                    neigh[name] = [set([name]),set(),set(),set()]
                    for base in neigh[name][0]:
                        for partners in bondnames:
                            partner_list = partners.split()
                            if base not in partner_list: continue
                            partner_list.remove(base)
                            neigh[name][1].add(partner_list[0])
                    neigh[name][1] = neigh[name][1].difference(neigh[name][0])
                    for base in neigh[name][1]:
                        for partners in bondnames:
                            partner_list = partners.split()
                            if base not in partner_list: continue
                            partner_list.remove(base)
                            neigh[name][2].add(partner_list[0])
                    neigh[name][2] = neigh[name][2].difference(neigh[name][0])
                    neigh[name][2] = neigh[name][2].difference(neigh[name][1])
                    for base in neigh[name][2]:
                        for partners in bondnames:
                            partner_list = partners.split()
                            if base not in partner_list: continue
                            partner_list.remove(base)
                            neigh[name][3].add(partner_list[0])
                    neigh[name][3] = neigh[name][3].difference(neigh[name][0])
                    neigh[name][3] = neigh[name][3].difference(neigh[name][1])
                    neigh[name][3] = neigh[name][3].difference(neigh[name][2])
                    #print >>f, "%8d %6d %6d %6d # %s" % (i+1, len(neigh[name][1]), len(neigh[name][2]), len(neigh[name][3]), name)
                    print >>f, "%8d %6d %6d %6d" % (i+1, len(neigh[name][1]), len(neigh[name][2]), len(neigh[name][3]))
                print >>f
                print >>f, "Special Bonds"
                print >>f
                for i, name in enumerate(atomnames):
                    print >>f, "%8d" % (i+1),
                    for partner in neigh[name][1]:
                        ipartner = atomnames.index(partner)
                        print >>f, "%8d" % (ipartner+1),
                    for partner in neigh[name][2]:
                        ipartner = atomnames.index(partner)
                        print >>f, "%8d" % (ipartner+1),
                    for partner in neigh[name][3]:
                        ipartner = atomnames.index(partner)
                        print >>f, "%8d" % (ipartner+1),
                    print >>f

    with open(inlmp, "w") as f:
        print >>f, "# LAMMPS command file generated by ForConX"
        print >>f
        print >>f, "units real"
        print >>f, "boundary p p p"
        print >>f
        print >>f, "molecule mymols",
        for item in root.findall('./molecule'):
            print >>f, "%s%s" % (item.get('name').lower(), mol_suffix),
        print >>f
        print >>f
        print >>f, "atom_style     ",
        if len(atom_style) > 1: print >>f, "hybrid",
        for style in atom_style: print >>f, style,
        print >>f
        if len(bond_types) > 0:
            print >>f, "bond_style     ",
            if len(bond_style) > 1: print >>f, "hybrid",
            for style in bond_style: print >>f, style,
            print >>f
        if len(angle_types) > 0:
            print >>f, "angle_style    ",
            if len(angle_style) > 1: print >>f, "hybrid",
            for style in angle_style: print >>f, style,
            print >>f
        if len(dihedral_types) > 0:
            print >>f, "dihedral_style ",
            if len(dihedral_style) > 1: print >>f, "hybrid",
            for style in dihedral_style: print >>f, style,
            print >>f
        if len(improper_types) > 0:
            print >>f, "improper_style ",
            if len(improper_style) > 1: print >>f, "hybrid",
            for style in improper_style: print >>f, style,
            print >>f
        if len(pair_types) > 0:
            print >>f, "pair_style     ",
            if len(pair_style) > 1: print >>f, "hybrid/overlay",
            for style, options in pair_style: print >>f, style, options,
            print >>f
        print >>f, "pair_modify mix %s" % mixing

        print >>f, "kspace_style   ", "pppm 1.e-4"
        print >>f, "kspace_modify  ", "gewald 0.41", "mesh 54 54 54"
        print >>f
        print >>f, "special_bonds lj 0. 0. %f coul 0. 0. %f" % (factor_lj, factor_coul)
        print >>f, "read_data data.lmp"

        if pairlmp is not None:
            g = open(pairlmp, "w")
            print >>g, pairtxt
            g.close()
            print >>f, "include", pairlmp
        else:
            print >>f
            print >>f, pairtxt

        print >>f
        if len(shaked_bonds) > 0:
            print >>f, "fix SHAKE all shake 1.e-4 100 0 b",
            for tp in shaked_bonds:
                itp = bond_types.index(tp)
                print >>f, itp+1,


