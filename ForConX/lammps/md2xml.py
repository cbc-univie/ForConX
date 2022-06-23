import numpy as np
import xml.etree.ElementTree as ET
from ..md_xml import input_output
from ..md_xml import molecule
from ..md_xml import bonds, angles, dihedrals, impropers, nonbonded

def readline(f):
    line = ""
    while line == "":
        try: line = f.next().strip()
        except StopIteration: return None
        if line == "": continue
        while line[-1] == '&':
            line += f.next().strip()
        end = 0
        mode = 'normal'
        for c in line:
            if c == '"':
                if mode == 'normal': mode = 'string'
                elif mode == 'string': mode = 'normal'
            elif c == '#':
                if mode == 'normal': break
            end += 1
        line = line[:end]
    return line

def read_mol(f):
    f.next() # skip first line
    line = ""
    word = ""
    natoms = nbonds = nangles = ndihedrals = nimpropers = 0
    while word not in ["Types", "Charges", "Bonds", "Angles", "Dihedrals", "Impropers"]:
        try:
            items = line.split()
            if items[1] == "atoms": natoms = int(items[0])
            elif items[1] == "bonds": nbonds = int(items[0])
            elif items[1] == "angles": nangles = int(items[0])
            elif items[1] == "dihedrals": ndihedrals = int(items[0])
            elif items[1] == "impropers": nimpropers = int(items[0])
        except (ValueError, IndexError):
            pass
        line = f.next().strip()
        if line != "": 
            word = line.split()[0]
        else: 
            word = ""
    typelist = []
    chargelist = []
    bondlist = []
    anglelist = []
    dihedrallist = []
    improperlist = []
    while True:
        if word == "Types":
            while len(typelist) < natoms:
                line = f.next().strip()
                if line is None: break
                elif line == "": continue
                items = line.split()
                typelist.append(int(items[1]))
        elif word == "Charges":
            while len(chargelist) < natoms:
                line = f.next().strip()
                if line is None: break
                elif line == "": continue
                items = line.split()
                chargelist.append(float(items[1]))
        elif word == "Bonds":
            while len(bondlist) < nbonds:
                line = f.next().strip()
                if line is None: break
                elif line == "": continue
                items = line.split()
                bondlist.append(map(int, tuple(items[1:4])))
        elif word == "Angles":
            while len(anglelist) < nangles:
                line = f.next().strip()
                if line is None: break
                elif line == "": continue
                items = line.split()
                anglelist.append(map(int, tuple(items[1:5])))
        elif word == "Dihedrals":
            while len(dihedrallist) < ndihedrals:
                line = f.next().strip()
                if line is None: break
                elif line == "": continue
                items = line.split()
                dihedrallist.append(map(int, tuple(items[1:6])))
        elif word == "Impropers":
            while len(improperlist) < nimpropers:
                line = f.next().strip()
                if line is None: break
                elif line == "": continue
                items = line.split()
                improperlist.append(map(int, tuple(items[1:6])))
        try: line = f.next().strip()
        except StopIteration: break
        if line != "": word = line.split()[0]
        else: word = ""
    # check if the bond/angle types always match with same atom types
    # don't check dihedrals/impropers because dihedral_style charmm and improper_style cvff are represented by multiple dihedrals/impropers
    bondmap = {}
    for bond in bondlist:
        for key in [(bond[1], bond[2]), (bond[2], bond[1])]:
            if key in bondmap:
                if bondmap[key] != bond[0]:
                    print "WARNING: inconsistent bond types"
            else: bondmap[key] = bond[0]
    anglemap = {}
    for angle in anglelist:
        for key in [(angle[1], angle[2], angle[3]), (angle[3], angle[2], angle[1])]:
            if key in anglemap:
                if anglemap[key] != angle[0]:
                    print "WARNING: inconsistent angle types"
            else: anglemap[key] = angle[0]
    return typelist, chargelist, bondlist, anglelist, dihedrallist, improperlist

def read_data(f, units, template):
    atoms = []
    nmols = {}
    masses = []
    bondtypes = []
    angletypes = []
    dihedraltypes = []
    impropertypes = []
    section = ""
    for line in f:
        line = line.strip()
        if line == "": continue
        word = line.split()[0]
        if word in ["Atoms", "Masses", "Bond", "Angle", "Dihedral", "Improper"]:
            section = word
        elif section == 'Atoms':
            items = line.split()
            if template != "":
                # iatom itype x y z q imol imoltype iinmol
                iatom = int(items[0])
                itype = int(items[1])
                x = float(items[2]) * units[0]
                y = float(items[3]) * units[0]
                z = float(items[4]) * units[0]
                q = float(items[5])
                imol = int(items[6])
                imoltype = int(items[7])
                iamol = int(items[8])
            else:
                # iatom imol itype q x y z
                iatom = int(items[0])
                imol = int(items[1])
                itype = int(items[2])
                q = float(items[3])
                x = float(items[4]) * units[0]
                y = float(items[5]) * units[0]
                z = float(items[6]) * units[0]
                imoltype = -1
                iamol = -1
            atoms.append((iatom, itype, x, y, z, q, imol, imoltype, iamol))
            if iamol == 1:
                if imoltype in nmols: nmols[imoltype] += 1
                else: nmols[imoltype] = 1
        elif section == 'Masses':
            masses.append(float(line.split()[1]))
        elif section == 'Bond':
            bondtypes.append(line.split()[1:])
        elif section == 'Angle':
            angletypes.append(line.split()[1:])
        elif section == 'Dihedral':
            dihedraltypes.append(line.split()[1:])
        elif section == 'Improper':
            impropertypes.append(line.split()[1:])
    # check
    chargemap = {}
    for a in atoms:
        itype = a[1]
        q = a[5]
        if itype in chargemap:
            if chargemap[itype] != q:
                print "WARNING inconsistent charges"
        else:
            chargemap[itype] = q
    return atoms, nmols, masses, bondtypes, angletypes, dihedraltypes, impropertypes

def read_input(f, root):
    moleculelist = {}
    bondstyle = set()
    anglestyle = set()
    dihedralstyle = set()
    improperstyle = set()
    pairstyle = []
    pairtypes = []
    datafilename = ""
    template = ""
    elec14 = vdw14 = None
    mixing = None
    line = readline(f)
    units = None # (distance, energy) conversion factors: multiply by them what you read in lammps to get a result in (Angstrom, kcal/mol)
    while line is not None:
        groups = line.split('"')
        items = []
        for i in range(len(groups)):
            if i % 2 == 0:
                items += groups[i].split()
            else:
                items.append(groups[i])
        command = items[0]
        if command == "units":
            co = input_output.mdElement(root,'output')
            if items[1] == "real": # angstrom, kcal/mol
                units = (co.convert_distance("ANGSTROEM"), co.convert_energy("KCAL"))
            elif items[1] == "metal": # angstrom, eV
                units = (co.convert_distance("ANGSTROEM"), co.convert_energy("EV"))
            else:
                print "WARNING, units are not known. I will do no conversion."
        elif command == "molecule":
            moleculelist[items[1]] = tuple(items[2:])
        elif command == "atom_style":
            if items[1:] == ["hybrid", "charge", "template"]:
                template = items[4]
            elif items[1:] != ["full"]:
                print " WARNING atom_style should be full or hybrid charge template"
        elif command == "bond_style":
            bondstyle = set(items[1:])
            if len(bondstyle) > 1: bondstyle.remove("hybrid")
        elif command == "angle_style":
            anglestyle = set(items[1:])
            if len(anglestyle) > 1: anglestyle.remove("hybrid")
        elif command == "dihedral_style":
            dihedralstyle = set(items[1:])
            if len(dihedralstyle) > 1: dihedralstyle.remove("hybrid")
        elif command == "improper_style":
            improperstyle = set(items[1:])
            if len(improperstyle) > 1: improperstyle.remove("hybrid")
        elif command == "pair_style":
            pairstyle = list(items[1:])
            if len(pairstyle) > 1:
                if pairstyle[0].startswith("hybrid"): # hybrid or overlay
                    pst = []
                    for st in pairstyle:
                        try: float(st)
                        except TypeError: pst.append(st)
                    pairstyle = pst
                else:
                    pairstyle = [pairstyle[0]]
        elif command == "special_bonds":
            if items[1] in ["coul", "lj/coul"]: elec14 = float(items[4])
            if items[1] in ["lj", "lj/coul"]:  vdw14 = float(items[4])
            if len(items) > 5:
                if items[5] == "coul": elec14 = float(items[8])
                elif items[5] == "lj": vdw14 = float(items[8])
        elif command == "read_data":
            datafilename = items[1]
        elif command == "pair_coeff":
            pairtypes.append(tuple(items[1:]))
        elif command == "pair_modify":
            if "mix" in items[1:]:
                mixing = items[items.index("mix")+1]
        elif command == "include":
            f1 = open(items[1], "r")
            moleculelist1, template1, bondstyle1, anglestyle1, dihedralstyle1, improperstyle1, pairtypes1, pairstyle1, elec141, vdw141, datafilename1, mixing1, units1 = read_input(f1, root)
            f1.close()
            moleculelist.update(moleculelist1)
            if len(bondstyle1) != 0: bondstyle = bondstyle1
            if len(anglestyle1) != 0: anglestyle = anglestyle1
            if len(dihedralstyle1) != 0: dihedralstyle = dihedralstyle1
            if len(improperstyle1) != 0: improperstyle = improperstyle1
            if len(pairstyle1) != 0: pairstyle = pairstyle1
            pairtypes += pairtypes1
            if elec141 is not None: elec14 = elec141
            if vdw141 is not None: vdw14 = vdw141
            if datafilename1 != "": datafilename = datafilename1
            if template1 != "": template = template1
            if mixing1 is not None: mixing = mixing1
            if units1 is not  None: units = units1
        line = readline(f)
    bondstyle = list(bondstyle)
    anglestyle = list(anglestyle)
    dihedralstyle = list(dihedralstyle)
    improperstyle = list(improperstyle)
    return moleculelist, template, bondstyle, anglestyle, dihedralstyle, improperstyle, pairtypes, pairstyle, elec14, vdw14, datafilename, mixing, units

def write_pdb(f, atoms, molnames):
    print >>f, "REMARK Created by ForConX"
    for iatom, itype, x, y, z, q, imol, imoltype, iamol in atoms:
        mname = molnames[imoltype-1][:4].upper()
        if '.' in mname: mname = mname[:mname.index('.')]
        print >>f, "ATOM %6d %4s %-4s %4d %8.3f %8.3f %8.3f 1.0 0.0 %-4s" % (iatom, "A%d" % iamol, mname, imol, x, y, z, mname)

def read(root):
    inputitem = root.find("input/command")
    inputfilename = inputitem.get("file")
    inputfile = open(inputfilename, "r")
    moleculelist, template, bondstyle, anglestyle, dihedralstyle, improperstyle, pairtypes, pairstyle, elec14, vdw14, datafilename, mixing, units = read_input(inputfile, root)
    if elec14 is None: elec14 = 0.
    if vdw14 is None:  vdw14 = 0.
    if mixing is None: mixing = 'geometric'
    if units is None: units = (1., 1.)
    inputfile.close()
    datafile = open(datafilename, "r")
    atoms, nmols, masses, bondtypes, angletypes, dihedraltypes, impropertypes = read_data(datafile, units, template)
    datafile.close()
    bondmap = {}
    anglemap = {}
    dihedralmap = {}
    impropermap = {}
    for imol, molfile in enumerate(moleculelist[template]):
        mname = molfile[:4].upper()
        if '.' in mname: mname = mname[:mname.index('.')]
        mol = root.find("molecule[@name=\"%s\"]" % mname)
        if mol is None:
            mol = ET.Element("molecule", name=mname)
            root.append(mol)
        else:
            mol.clear()
            mol.set("name", mname)
        mol = molecule.moleculeElement(mol, mname)
        mol.name = mname
        mol.nmol = nmols[imol+1]
        mol.mol.set("nmol", str(mol.nmol))
        mf = open(molfile, "r")
        typelist, chargelist, bondlist, anglelist, dihedrallist, improperlist = read_mol(mf)
        mf.close()
        diheddone = set()
        improdone = set()
        for i in range(len(typelist)):
            atom = molecule.atomClass(mol.mol, "A%d" % (i+1))
            atom.type = type="T%d" % typelist[i]
            atom.charge = chargelist[i]
            atom.mass = masses[typelist[i]-1]
        for bond in bondlist:
            if bond[-1] < bond[1]: bond[1:] = bond[:0:-1]
            bondname = " ".join(["A%d" % i for i in bond[1:]])
            bondelt = molecule.bondClass(mol.mol, bondname)
            bmap = bonds.sequence(" ".join(["T%d" % typelist[i-1] for i in bond[1:]]))
            bondelt.type = bmap
            bondmap[bond[0]] = bmap
        for angle in anglelist:
            if angle[-1] < angle[1]: angle[1:] = angle[:0:-1]
            anglename = " ".join(["A%d" % i for i in angle[1:]])
            angleelt = molecule.angleClass(mol.mol, anglename)
            amap = angles.sequence(" ".join(["T%d" % typelist[i-1] for i in angle[1:]]))
            angleelt.type = amap
            anglemap[angle[0]] = amap
        for dihedral in dihedrallist:
            if dihedral[-1] < dihedral[1]: dihedral[1:] = dihedral[:0:-1]
            dihedralname = " ".join(["A%d" % i for i in dihedral[1:]])
            dihedralelt = molecule.dihedralClass(mol.mol, dihedralname)
            dmap = dihedrals.sequence(" ".join(["T%d" % typelist[i-1] for i in dihedral[1:]]))
            dihedralelt.type = dmap
            dihedralmap[dihedral[0]] = dmap
            if dihedralname in diheddone: continue
            diheddone.add(dihedralname)
        for improper in improperlist:
            # look for the central atom (assuming it is unique)
            for central in [1,2,3,4]:
                ok = True
                for neigh in improper[2:]:
                    found = False
                    for bond in bondlist:
                        if bond[1] == neigh and bond[2] == improper[central] or bond[1] == improper[central] and bond[2] == neigh:
                            found = True
                            break
                    if not found:
                        ok = False
                        break
                if ok: break
            if central != 1:
                if central != 4: print "WARNING: central atom of LAMMPS impropers should be first or last"
                improper[1:] = improper[:0:-1]
            central = 1
            impropername = " ".join(["A%d" % i for i in improper[1:]])
            improperelt = molecule.improperClass(mol.mol, impropername)
            dmap = " ".join(["T%d" % typelist[i-1] for i in improper[1:]])
            improperelt.type = dmap
            impropermap[improper[0]] = dmap
            if impropername in improdone: continue
            improdone.add(impropername)
            improperelt.central = "A%d" % improper[central]
    bondlist = root.find("bonds")
    if bondlist is not None: bondlist.clear()
    for ibondtype, bondtype in enumerate(bondtypes):
        if len(bondstyle) > 1:
            bondst = bondtype[0]
            bondtype = bondtype[1:]
        else:
            bondst = bondstyle[0]
        if bondst == 'harmonic':
            bond = bonds.harmClass(root, bondmap[ibondtype+1])
            bond.k = float(bondtype[0]) * units[1]/units[0]**2
            bond.r0 = float(bondtype[1]) * units[0]
    anglelist = root.find("angles")
    if anglelist is not None: anglelist.clear()
    for iangletype, angletype in enumerate(angletypes):
        if len(anglestyle) > 1:
            anglest = angletype[0]
            angletype = angletype[1:]
        else:
            anglest = anglestyle[0]
        if anglest == 'harmonic':
            angle = angles.harmClass(root, anglemap[iangletype+1])
            angle.k = float(angletype[0]) * units[1]
            angle.theta0 = float(angletype[1])
    dihedrallist = root.find("dihedrals")
    if dihedrallist is not None: dihedrallist.clear()
    diheddone = set()
    for idihedraltype, dihedraltype in enumerate(dihedraltypes):
        if len(dihedralstyle) > 1:
            dihedralst = dihedraltype[0]
            dihedraltype = dihedraltype[1:]
        else:
            dihedralst = dihedralstyle[0]
        dmap = dihedralmap[idihedraltype+1]
        if dmap in diheddone: continue
        diheddone.add(dmap)
        if dihedralst == 'charmm':
            ks = []
            ns = []
            ds = []
            for id2, d2 in enumerate(dihedraltypes):
                if dihedralmap[id2+1] == dmap:
                    if len(dihedralstyle) > 1:
                        d2st = d2[0]
                        d2 = d2[1:]
                    else:
                        d2st = dihedralstyle[0]
                    if d2st != 'charmm': print "WARNING: multiple dihedrals with different styles"
                    ks.append(float(d2[0]) * units[1])
                    ns.append(float(d2[1]))
                    ds.append(float(d2[2]))
                    if float(d2[3]) != 0.: print "WARNING: factor 14 should be defined in special_bonds, not in pair style charmm"
            dihedral = dihedrals.cosClass(root, dmap)
            dihedral.k = ks
            dihedral.n = ns
            dihedral.delta = ds
        elif dihedralst == 'opls': # convert to Ryckaert-Bellemans
            ks = [0.5 * float(dihedraltype[0]) * units[1],
                 -0.5 * float(dihedraltype[0]) * units[1],
                  0.5 * float(dihedraltype[0]) * units[1], 
                 -0.5 * float(dihedraltype[0]) * units[1]]
            shift = ks[0] - ks[1] + ks[2] - ks[3]
            dihedral = dihedrals.ryckClass(root, dmap)
            dihedral.k = [shift] + ks
    improperlist = root.find("impropers")
    if improperlist is not None: improperlist.clear()
    improdone = set()
    for iimpropertype, impropertype in enumerate(impropertypes):
        if len(improperstyle) > 1:
            improperst = impropertype[0]
            impropertype = impropertype[1:]
        else:
            improperst = improperstyle[0]
        if improperst == 'cvff':
            dmap = impropermap[iimpropertype+1]
            if dmap in improdone: continue
            improdone.add(dmap)
            ks = []
            ns = []
            ds = []
            for id2, d2 in enumerate(impropertypes):
                if impropermap[id2+1] == dmap:
                    if len(improperstyle) > 1:
                        d2st = d2[0]
                        d2 = d2[1:]
                    else:
                        d2st = improperstyle[0]
                    if d2st != 'cvff': print "WARNING: multiple impropers with different styles"
                    ks.append(float(d2[0]) * units[1])
                    ns.append(float(d2[1]))
                    ds.append(float(d2[2]))
            improper = impropers.cosClass(root, dmap)
            improper.k = ks
            improper.n = ns
            improper.delta = ds
        elif improperst == 'harmonic':
            improper = impropers.harmClass(root, impropermap[iimpropertype+1])
            improper.k = float(impropertype[0]) * units[1]
            improper.theta0 = float(impropertype[1])
    nonbondeditem = root.find("nonbonded")
    if nonbondeditem is not None: nonbondeditem.clear()
    vdwmap = {}
    vdws = []
    for p in pairtypes:
        iborn = map(int, p[0].split("*"))
        jborn = map(int, p[1].split("*"))
        ivdw = len(vdws)
        iarg = 2
        while iarg < len(p):
            if len(pairstyle) > 1:
                pstyle = p[iarg]
                iarg += 1
            else:
                pstyle = pairstyle[0]
            if pstyle.startswith("lj/charmm/coul") or pstyle.startswith("lj/cut") or pstyle.startswith("lj/long"):
                coeffs = []
                coeffs.append(float(p[iarg]) * units[1]) # epsilon
                iarg += 1
                coeffs.append(float(p[iarg]) * units[0]) # sigma
                iarg += 1
                try: # todo: use 14 values if charmm?
                    float(p[iarg])
                    iarg += 1
                    float(p[iarg])
                    iarg += 1
                except (ValueError, IndexError):
                    pass
                vdws.append(coeffs)

        for i in range(iborn[0], iborn[-1]+1):
            for j in range(jborn[0], jborn[-1]+1):
                if j < i: continue
                vdwmap[i,j] = ivdw
    for key, val in vdwmap.items():
        tp = "T%d T%d" % key
        vdw = nonbonded.vdwClass(root, tp)
        vdw.epsilon = vdws[val][0]
        vdw.sigma = vdws[val][1]
        vdw.elec14 = elec14
        vdw.vdw14 = vdw14
    nonbondedelt = nonbonded.nonbondedElement(root)
    nonbondedelt.mixing_epsilon = "geometric"
    nonbondedelt.mixing_sigma = mixing

    pdbfile = open("forconx.pdb", "w")
    write_pdb(pdbfile, atoms, moleculelist[template])
    pdbfile.close()
    #ET.dump(root)


