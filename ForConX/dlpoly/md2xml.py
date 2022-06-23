import xml.etree.ElementTree as ET
#======================================================================================
# md2xml.read
#======================================================================================
def read(root):
    """
    """
    from ..md_xml import helpers
    
    current_input = root.find('./input/field')
    try:
        field_file = current_input.get('file')
        f = open(field_file)
        print '> Field file = %s ...'%field_file
    except IOError:
        helpers.error('Field file %s cannot be opened!'%field_file)
    except:
        helpers.error('<field file="..."> is ill-defined!')
    
    # Reading header    
    while 1:
        buffer = f.readline()
        buffer = buffer.upper()
        if buffer[0:1] == "#":
            continue
        if buffer[0:5] == "UNITS":
            line_element = buffer.split()
            energy_unit = line_element[1].upper()
            pointer = root.find('input/energy')
            if pointer is None:
                pointer = ET.SubElement(root.find('input'),'energy')
            pointer.set('unit',energy_unit)
            print "\t\tEnergy unit = ",energy_unit
        if buffer[0:5] == "CLOSE":
            break
        if buffer[0:7] == "MOLECUL":
            break

    # Reading molecule
    print"\n\t----------------------------------"
    print"\t1.1 Reading molecules"
    print"\t----------------------------------"
    moleculelist = []
    residue   = {}    
    line_element = buffer.split()
    while 1:
        rewind_position = f.tell()
        buffer = f.readline()
        buffer = buffer.upper()
        if buffer[0:3]=="VDW":
            break
        else:
            f.seek(rewind_position)
            buffer = f.readline()
            line_element = buffer.split()
            molname = line_element[0].strip()
            print '\t\t<molecule name="%s">'%molname
            topology = []
            while 1:
                buffer = f.readline()
                buffer = buffer.upper()
                if buffer[0:6]=="FINISH":
                    break
                topology.append(buffer.rstrip())
            residue[molname] = topology
            moleculelist.append(molname)

    # Reading vdw pairs
    print"\n\t----------------------------------"
    print"\t1.2 Reading van-der-Waals"
    print"\t----------------------------------"
    vdw = []    
    while 1:
        buffer = f.readline()
        buffer = buffer.upper()
        if buffer[0:5]=="CLOSE":
            break
        vdw.append(buffer.rstrip())
    print"\t\tNumber of vdW - pair interactions = ",len(vdw)    
    f.close()
    
    # Converting corresponding dictionaries to XML
    print"\n\t----------------------------------"
    print"\t1.3 Converting molecules"
    print"\t----------------------------------"
    from ..md_xml import helpers
    for mol in root.findall('molecule'):
        molname = (mol.get("name")).strip()
        print'\t\t<molecule name=\"%s\">' %molname
        if molname not in residue.keys():
            helpers.error('Molecule is not found in FIELD file!')
        for i in mol.findall('*'):
            mol.remove(i)
            
        parameter = decompose_residue(root,mol,residue[molname])    
        atomClass = read_atom(mol,molname,parameter["atoms"])
        print"\t\t\tNumber of atoms     = %5s" % len(atomClass)

        if parameter.has_key("shell"):
            print "\t\t\tMaking atoms polarizable ..."
            pointer = root.find('input/polarizability')
            if pointer is None:
                pointer = ET.SubElement(root.find('input'),'polarizability')
            pointer.set('unit','ANGSTROEM')
            read_shell(root,mol,parameter["shell"],atomClass)
        
        nbond = read_bond(root,mol,parameter["bonds"],atomClass)
        print"\t\t\tNumber of bonds     = %5s" % nbond

        nangle = read_angle(root,mol,parameter["angles"],atomClass)
        print"\t\t\tNumber of angles    = %5s" % nangle

        (ndihedral,nimproper) = read_dihedral(root,mol,parameter["dihedrals"],atomClass)
        print"\t\t\tNumber of dihedrals = %5s  " %ndihedral
        print"\t\t\tNumber of impropers = %5s\n" %nimproper

    print"\n\t----------------------------------"
    print"\t1.4 Converting nonbonded interactions"
    print"\t----------------------------------"
    read_vdw(root,vdw)
    
#===============================================================================================
# md2xml.decompose_residue
#===============================================================================================
def decompose_residue(root,mol,residue):
    """
    """
    parameter = {}
    name = "dump"
    parameter["dump"] = []
    for line in residue:
        buffer = line.upper()
        if buffer[0:7]=="NUMMOLS":
            mol.set('nmol',str(int(buffer[8:].strip())))
        if buffer[0:5]=="ATOMS":
            name = "atoms"
            parameter["atoms"] = []
            continue
        if buffer[0:5]=="BONDS":
            name = "bonds"
            parameter["bonds"] = []
            continue
        if buffer[0:11]=="CONSTRAINTS":    
            name = "bonds"
            continue
        if buffer[0:5]=="SHELL":
            name = "shell"
            parameter["shell"] = []
            continue
        if buffer[0:6]=="ANGLES":    
            name = "angles"
            parameter["angles"] = []
            continue
        if buffer[0:9]=="DIHEDRALS":
            name = "dihedrals"
            parameter["dihedrals"] = []
            continue
        parameter[name].append(buffer.upper())    
    return parameter

#===============================================================================================
# md2xml.read_atom
#===============================================================================================
def read_atom(mol,molname,parameter_atoms):
    """
    """
    from ..md_xml import molecule
    name_list = {}
    atomClass = []
    current_molecule = molecule.moleculeElement(mol,molname)
    for line in parameter_atoms:
        line_element = line.split()
        type = line_element[0]            
        mass = float(line_element[1])
        charge = float(line_element[2])
        repeat = int(line_element[3])
        for i in range(repeat):
            atomname = current_molecule.unique_atomname(type)
            # Building atom objects
            current_atom = molecule.atomClass(mol,atomname)
            current_atom.type = type
            current_atom.mass = mass
            current_atom.charge = charge
            atomClass.append(current_atom)
    return atomClass

#===============================================================================================
# md2xml.read_shell
#===============================================================================================
def read_shell(root,mol,parameter_shell,atomClass):
    """
    """
    from ..md_xml import molecule

    for line in parameter_shell:
        line_element = line.split()
        index_i = int(line_element[0])
        index_j = int(line_element[1])
        current_i = atomClass[index_i-1]
        current_j = atomClass[index_j-1]
        mass_i = current_i.mass
        mass_j = current_j.mass
        # First atom i should be a real atom        
        if mass_i < mass_j:
            tmp = current_i
            current_i = current_j
            current_j = tmp
        charge_i = current_i.charge
        charge_j = current_j.charge
        current_i.charge = charge_i + charge_j
        current_i.mass   =  mass_i  + mass_j

        # alpha = q^2 / 4 pi epsilon0 k * eV factor
        k = float(line_element[2])
        current_i.alpha = charge_j*charge_j*1602.2/(4*3.141592654*8.85*k)
        current_j.remove()
    return

#===============================================================================================
# md2xml.read_bond
#===============================================================================================
def read_bond(root,mol,parameter_bonds,atomClass):
    """
    """
    from ..md_xml import bonds
    from ..md_xml import molecule
    from ..md_xml import helpers 
    
    nbond = 0
    for line in parameter_bonds:
        line_element = line.split()
        if len(line_element)==3:
            # constraints            
            index_i = int(line_element[0])
            index_j = int(line_element[1])
        else:
            # bonds            
            index_i = int(line_element[1])
            index_j = int(line_element[2])
        current_i = atomClass[index_i-1]
        current_j = atomClass[index_j-1]
        try:
            name = ' '.join([current_i.name,current_j.name])
            type = ' '.join([current_i.type,current_j.type])
        except:
#           Bond between shell atoms            
            continue
            
        if line_element[0]=="HARM":
            nbond += 1
            current_bond = molecule.bondClass(mol,name)
            current_bond.type = type
            bond    = bonds.harmClass(root,type)
            bond.k  = float(line_element[3]) / 2.
            bond.r0 = float(line_element[4])
        elif len(line_element)==3:
            nbond += 1
            current_bond = molecule.bondClass(mol,name)
            current_bond.type = type
            bond = bonds.harmClass(root,type)
            bond.k = None
            bond.r0 = float(line_element[2])
        elif line_element[0]=="MORS":
            nbond += 1
            current_bond = molecule.bondClass(mol,name)
            current_bond.type = type
            bond = bonds.morsClass(root,type)
            bond.D0   = float(line_element[3])
            bond.r0   = float(line_element[4])
            bond.beta = float(line_element[5])  
        else:
            helpers.warning('The following potential is not supported by ForConX!\n%s'%line)
    return nbond

#===============================================================================================
# md2xml.read_angle
#===============================================================================================
def read_angle(root,mol,parameter_angles,atomClass):
    """
    """
    from ..md_xml import angles
    from ..md_xml import molecule
    from ..md_xml import helpers
    
    nangle = 0
    for line in parameter_angles:
        line_element = line.split()
        current_i = atomClass[int(line_element[1])-1]
        current_j = atomClass[int(line_element[2])-1]
        current_k = atomClass[int(line_element[3])-1]
        name = ' '.join([current_i.name,current_j.name,current_k.name])
        type = ' '.join([current_i.type,current_j.type,current_k.type])
        
        if line_element[0]=="HARM":
            nangle += 1
            current_angle = molecule.angleClass(mol,name)
            current_angle.type = type
            angle        = angles.harmClass(root,type)
            angle.k      = float(line_element[4])/2.
            angle.theta0 = float(line_element[5]) 
        elif line_element[0]=="HCOS":
            import numpy
            helpers.warning('Harmonic cosine %s converted to harmonic potential.'%type) 
            nangle += 1
            current_angle = molecule.angleClass(mol,atomnames)
            current_angle.type = type
            angle = angles.harmClass(root,atomtypes)
            angle.theta0 = float(line_element[5])                        
            angle.k = float(line_element[4])/2. * numpy.sin(angle.theta0*3.14159265/180.)
        else:
            helpers.warning('The following potential is not supported by ForConX!\n%s\n'%line)
    return nangle

#===============================================================================================
# md2xml.read_dihedrals
#===============================================================================================
def read_dihedral(root,mol,parameter_dihedrals,atomClass):
    """
    """
    from ..md_xml import molecule
    from ..md_xml import bonds
    from ..md_xml import dihedrals
    from ..md_xml import impropers
    from ..md_xml import nonbonded
    from ..md_xml import helpers 
    
    ndihedral = 0
    nimproper = 0
    for line in parameter_dihedrals:
        line_element = line.split()
        atom_i = atomClass[int(line_element[1])-1]
        atom_j = atomClass[int(line_element[2])-1]
        atom_k = atomClass[int(line_element[3])-1]
        atom_l = atomClass[int(line_element[4])-1]
        elec14 = float(line_element[8])

        # dihedral
        if elec14 > 0.:
            ndihedral += 1
            name = dihedrals.sequence(' '.join([atom_i.name,atom_j.name,atom_k.name,atom_l.name]))
            type = dihedrals.sequence(' '.join([atom_i.type,atom_j.type,atom_k.type,atom_l.type]))
            current_dihedral = molecule.dihedralClass(mol,name)
            current_dihedral.type = type

            type_14 = bonds.sequence(' '.join([atom_i.type,atom_l.type]))
            current_vdw = nonbonded.vdwClass(root,type_14)
            current_vdw.elec14 = float(line_element[8])
            current_vdw.vdw14  = float(line_element[9])

            # Ryckaert-Bellemans, Chem. Phys. Lett. 30 (1975), 123  (untestest)  
            if line_element[0]=="RYCK":
                a = float(line_element[5]) 
                current_ryck   = dihedrals.ryckClass(root,type)
                current_ryck.k = [1116*a, 1462*a, -1578*a, -368*a, 3156*a, -3788*a]
            elif line_element[0]=="COS":
                current_cos = dihedrals.cosClass(root,type)
                current_cos.n = [ float(line_element[7]) ]
                current_cos.k = [ 0.5*float(line_element[5]) ]
                current_cos.delta = [ float(line_element[6]) ]
            elif line_element[0]=="COS3":
                current_cos = dihedrals.cosClass(root,type)
                current_cos.n = [ 1  , 2  , 3  ]
                current_cos.k = [ 0.5*float(line_element[5]),
                                  0.5*float(line_element[6]),
                                  0.5*float(line_element[7]) ]
                current_cos.delta = [ 0., 180., 0. ] 
            elif line_element[0]=="OPLS":
                current_cos = dihedrals.cosClass(root,type)
                current_cos.n = [ 1 , 2 , 3 ]
                current_cos.k = [ float(line_element[6]),
                                  float(line_element[7]),
                                  float(line_element[10]) ]
                current_cos.delta = [ float(line_element[11]),
                                      2.*float(line_element[11]),
                                      3.*float(line_element[11]) ]
            else:
                helpers.warning('The following dihedral potential is not supported by ForConX!\n%s'%line)

        # improper
        else:
            nimproper += 1
            name = impropers.sequence(' '.join([atom_j.name,atom_i.name,atom_k.name,atom_l.name]))
            type = impropers.sequence(' '.join([atom_j.type,atom_i.type,atom_k.type,atom_l.type]))
            
            current_improper = molecule.improperClass(mol,name)
            current_improper.central = atom_j.name
            current_improper.type = type

            if line_element[0]=="HCOS":
                import numpy
                helpers.warning("Harmonic cosine of %s converted to harmonic angle"%type)
                current_harm  = impropers.harmClass(root,type)
                
                # theta0 degree -> radian
                theta0 = float(line_element[5])*0.017453293
                current_harm.theta0 = float(line_element[6])
                current_harm.k      = float(line_element[5])*numpy.sin(theta0)*numpy.sin(theta0)*0.5
            elif line_element[0]=="HARM":
                current_harm  = impropers.harmClass(root,type)
                current_harm.theta0 = float(line_element[6])
                current_harm.k      = float(line_element[5])*0.5
            elif line_element[0]=="COS3":
                current_cos = impropers.cosClass(root,type)
                current_cos.n = [ 1  , 2  , 3  ]
                current_cos.k = [ 0.5*float(line_element[5]),
                                  0.5*float(line_element[6]),
                                  0.5*float(line_element[7]) ]
                current_cos.delta = [ 0., 180., 0. ] 
            else:
                helpers.warning('The following potential is not supported by ForConX!\n%s\n'%line)
    return ndihedral,nimproper

#===============================================================================================
# md2xml.read_vdw
#===============================================================================================
def read_vdw(root,vdw):
    """
    """
    from ..md_xml import nonbonded
    current_nonbonded = nonbonded.nonbondedElement(root)
    current_nonbonded.mixing_epsilon = "geometric"
    current_nonbonded.mixing_sigma   = "arithmetic"

#   Iterating over all vdW-interactions
    print"\t\tLennard-Jones pairs ...",
    for line in vdw:
        line_element = (line.upper()).split()
        if line_element[2]=="LJ" or line_element[2]=="-LJ":
            current_epsilon = float(line_element[3])
            current_sigma   = float(line_element[4])
        if line_element[2]=="12-6" or line_element[2]=="-126":
            current_A = float(line_element[3])
            current_B = float(line_element[4])
            sigma6 = current_A/current_B
            current_sigma = sigma6**0.16666666666   
            current_epsilon = current_B/(4.*sigma6)
            
        (type_i, type_j) = sorted([line_element[0],line_element[1]])
        vdw = nonbonded.vdwClass(root,' '.join([type_i,type_j]))
        vdw.epsilon = current_epsilon
        vdw.sigma   = current_sigma
    print len(root.findall("nonbonded/vdw"))


