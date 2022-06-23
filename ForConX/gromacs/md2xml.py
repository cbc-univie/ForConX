#!/usr/bin/python
import xml.etree.ElementTree as ET
import re,sys
import os.path
from operator import itemgetter, attrgetter

#==================================================================================================
# md2xml.preprocess
#==================================================================================================
def preprocess(file):
    """
    preprocess handles include statements in the topology file.
    All input files will be merged and comment and blank lines deleted.
    """
    lines = []
    f = open(file)
    for line in f:
        line = (line.split(";"))[0].strip()
        if len(line)<2:
            continue
        if line[0:8]=="#include":
            line = line.split()
            filename = re.sub('\"', '', line[1])
            print ">    Open %s ..."%filename
            lines += (preprocess(filename))
        else:
            lines.append(line) 
    f.close()
    return lines
    
#==================================================================================================
# md2xml.read_top
#==================================================================================================
def read_top(root,f):
    """
    This routine reads the complete topology file.
    Due to a slightly different file format, OPLS force field are handled with special care.

    """
    from ..md_xml import helpers
    from ..md_xml import molecule
    
    print"\n\t-----------------------------------------"
    print"\t 1.1 Reading all itp-force field files"
    print"\t-----------------------------------------"
    print">\t Reading keywords"
    (toplines,fudge_lj,fudge_qq,combination_rule) = decompose_top(f)
    print"\n\t\t fudge_lj = %10.5f"%fudge_lj
    print"\t\t fudge_qq = %10.5f\n"%fudge_qq

    # OPLS keyword changes the force field structure
    oplstype = {}
    if "[ FF_OPLS ]" in toplines["ATOMTYPES"]:
        for line in toplines["ATOMTYPES"]:
            line_element = line.split()
            if "opls" in line_element[0]:
                oplstype[line_element[0]] = line_element[1]
    
    print"\n\t---------------------------------------------------"
    print"\t1.2 Converting molecules to XML"
    print"\t---------------------------------------------------"
    for mol in root.findall('molecule'):
        molname = mol.get('name')
        moltype = "MOLECULETYPE "+ molname
        if moltype not in toplines.keys():
            helpers.error('<molecule name="%s"> cannot be found in topology file'%molname)
        else:
            # Clear existing information        
            for i in mol.findall('*'):
                mol.remove(i)
        print'\n\t<molecule name="%s">'%molname
        residue = decompose_residue(toplines[moltype])

        # Read atoms and create <molecule name=molname> 
        (name_list,natom) = read_atom(mol,residue["ATOMS"],oplstype)
        print "\t\tNumber of atoms     = %5s" % natom

        # Bonds
        (nbond,toplines["BONDTYPES"]) = read_molecule_bond(mol,name_list,residue["BONDS"],toplines["BONDTYPES"],residue["ATOMS"])
        print"\t\tNumber of bonds     = %5s" % nbond
        
        # Angles
        (nangle,toplines["ANGLETYPES"]) = read_molecule_angle(mol,name_list,residue["ANGLES"],toplines["ANGLETYPES"])
        print"\t\tNumber of angles    = %5s" % nangle

        # Dihedrals and impropers are handled equally in GROMACS
        if "DIHEDRALS" in residue.keys():
            (ndihedral, nimproper,toplines["DIHEDRALTYPES"]) = read_molecule_dihedral_improper(mol,name_list,residue["DIHEDRALS"],toplines["DIHEDRALTYPES"])
            print"\t\tNumber of dihedrals = %5s" % ndihedral
            print"\t\tNumber of impropers = %5s" % nimproper

        # Pairs information are added to the nonbonded dictionary
        if "PAIRS" in residue.keys():
            toplines["ATOMTYPES"].append("[ moleculetype %s ]"%molname)
            toplines["ATOMTYPES"].append(residue["PAIRS"])
        
    print"\n\t--------------------------------------------------------------------------------------------------"
    print"\t1.3 Converting force field parameters to XML"
    print"\t--------------------------------------------------------------------------------------------------"
    nbond = read_bonds(root,toplines["BONDTYPES"])
    print"\t\tNumber of bonds     = %5s" % nbond

    nangle = read_angles(root,toplines["ANGLETYPES"])
    print"\t\tNumber of angles    = %5s" % nangle

    ndihedral = read_dihedrals(root,toplines["DIHEDRALTYPES"])
    print"\t\tNumber of dihedrals = %5s" % ndihedral

    nimpropers = read_impropers(root,toplines["DIHEDRALTYPES"])
    print"\t\tNumber of impropers = %5s" % nimpropers
        
    print"\n\t--------------------------------------------------------------------------------------------------"
    print"\t1.4 Assigning non-bonded parameter"
    print"\t--------------------------------------------------------------------------------------------------"
    iatoms = van_der_waals(root,toplines["ATOMTYPES"],fudge_lj,fudge_qq,combination_rule)
    print "\t\tNumber of nonbonded parameter   = %6s" % iatoms
    return
        
#==================================================================================================
# md2xml.decompose_topology
#==================================================================================================
def decompose_top(f):
    """
    The overall topology file is decomposed into several sections which are used in the 
    subsequent subroutines.
    """
    from ..md_xml import helpers
    fudge_lj = 0.0
    fudge_qq = 0.0
    combination_rule = 0.0
    
    toplines = {}
    flag = "DUMP"
    toplines["ATOMTYPES"] = []
    toplines["BONDTYPES"] = []
    toplines["ANGLETYPES"] = []
    toplines["DIHEDRALTYPES"] = []
    for i in range(len(f)):
        line = f[i]
        if "system" in line:
            print "\t\t[ system ] ..."
            flag = "DUMP"
            continue
        if "defaults" in line:
            print "\t\t[ defaults ] ..."
            flag = "DUMP"
            i += 1
            line_element = f[i].split()
            fudge_lj = float(line_element[3])
            fudge_qq = float(line_element[4])
            combination_rule = int(line_element[1])
            continue
        if "moleculetype" in line:
            i += 1
            molname = f[i].split()[0]
            print "\t\t[ moleculetype ] ... ",molname
            flag = "MOLECULETYPE "+molname
            continue
        if "#define _FF_OPLS" in line:
            print "\t\t OPLS force field"
            toplines["ATOMTYPES"].append("[ FF_OPLS ]")
        if "atomtypes" in line:
            print "\t\t[ atomtypes ] ..."
            flag = "ATOMTYPES"
        if "bondtypes" in line:
            print "\t\t[ bondtypes ] ..."
            flag = "BONDTYPES"
            continue
        if "angletypes" in line:
            print "\t\t[ angletypes ] ..."
            flag = "ANGLETYPES"
            continue
        if "dihedraltypes" in line:
            print "\t\t[ dihedraltypes ] ..."
            flag = "DIHEDRALTYPES"
            continue
        if "constrainttypes" in line:
            print "\t\t[ constrainttypes ] ..."
            flag = "CONSTRAINTSTYPES"
            continue

        if flag not in toplines.keys():
            toplines[flag] = []
        toplines[flag].append(line)

    if fudge_lj == 0.0 or fudge_qq == 0.0 or combination_rule == 0:
        helpers.error("[ defaults ] is ill-defined!")
            
    return toplines,fudge_lj, fudge_qq, combination_rule

#==================================================================================================
# md2xml.decompose_residue
#==================================================================================================
def decompose_residue(toplines):
    """
    decompose molecular information in bond, angle, dihedral
    """
    residue = {}
    flag = "DUMP"
    for line in toplines:
        if "atoms" in line:
            flag = "ATOMS"
            continue
        if "bonds" in line:
            flag = "BONDS"
            continue
        if "angles" in line:
            flag = "ANGLES"
            continue
        if "dihedrals" in line:
            flag = "DIHEDRALS"
            continue
        if "pairs" in line:
            flag = "PAIRS"
            continue
#### added thole and exclusions
        if "thole_polarization" in line:
            flag = "THOLE"
	    continue
        if "exclusions" in line:
            flag = "EXCLUSIONS"
	    continue
        if flag not in residue.keys():
            residue[flag] = []
        residue[flag].append(line)

    if len(residue["ANGLES"])==0:
        residue["ANGLES"] = "AUTOGENERATE"
    return residue

#==================================================================================================
# md2xml.read_atom
#==================================================================================================
def read_atom(mol,atom,oplstype):
    """
    Virtual atoms are detected by the mass.
    Drude atoms are detected by the type.
    """
    from ..md_xml import molecule

    natom = 0
    name_list = {}
    type_list = {}

    current_molecule = molecule.moleculeElement(mol,mol.get('name'))
    for line in atom:
    	dcheck = False
        line_element = line.split()
        index        = line_element[0]
        index2       = int(line_element[0])
        atomname     = current_molecule.unique_atomname(line_element[4])
        dtype        = line_element[1]
        mass         = float(line_element[7])
        charge       = float(line_element[6])

        if mass>0.1 and index2<len(atom):
            line2 = atom[index2]
            line2_element = line2.split()
            atomname2     = current_molecule.unique_atomname(line2_element[4])
            if atomname2 == "D"+atomname:
                current_atom = molecule.atomClass(mol,atomname)
                dmass = float(line2_element[7])
                dcharge = float(line2_element[6])
                alpha  = (dcharge*dcharge) / 3013.49
                current_atom.alpha = alpha
                current_atom.mass = mass + dmass
                current_atom.charge = charge + dcharge
            elif dtype == "D":
                dcheck = True
            else:
                current_atom = molecule.atomClass(mol,atomname)
                current_atom.mass = mass
                current_atom.charge = charge
        elif mass>0.1 and index2==len(atom):
            if dtype == "D":
                dcheck = True
                continue
            else:
                current_atom = molecule.atomClass(mol,atomname)
                current_atom.mass = mass
                current_atom.charge = charge

        else:
            # Virtual particle        
            current_atom = molecule.virtualClass(mol,atomname)
            current_atom.charge = charge
        # conversion to opls types?    
        if dcheck:
            continue
        else:
            if line_element[1] in oplstype.keys():
                current_atom.type = oplstype[line_element[1]]
            else:
                current_atom.type = line_element[1]
            name_list[index] = atomname
            type_list[index] = current_atom.type
            natom += 1

    return name_list,natom


#==================================================================================================
# md2xml.read_molecule_bond
#==================================================================================================
def read_molecule_bond(mol,name_list,bond,bondtypes,atom):
    """
    This routine adds all <molecule/bond> entries to the XML structure.

    Since GROMACS offers two possibilities to store bond potential information:
    (a) [ bonds ]
    (b) [ bondtypes ]

    and [ bonds ] is processed here, the force field parameters of this section is added to
    bondtypes and processes in the corresponding function
    
    Bonds between atoms and drudes are not extra listed in XML.
    """
    from ..md_xml import molecule
    from ..md_xml import helpers
    from ..md_xml import bonds

    nbond = 0
    for line in bond:
        line_element = line.split()
        name_i = name_list[line_element[0]]
        
        index = int(line_element[1])-1
        line2  = atom[index]
        line2_element = line2.split()
        
        if line2_element[4] == "D"+name_i:
            continue
        else:
            name_j = name_list[line_element[1]]
            atomnames = bonds.sequence(' '.join([name_i,name_j]))
            
            atom_i = molecule.atomClass(mol,name_i)
            atom_j = molecule.atomClass(mol,name_j)
            current_bond = molecule.bondClass(mol,atomnames)
            current_bond.type = bonds.sequence(' '.join([atom_i.type,atom_j.type]))
            nbond += 1

        # if bond potential are defined here, they will be added to bondtypes
        # and handled later
        if len(line_element)>3:
            if line_element[2] == "1":
                tmp =  "%10s 1 " %current_bond.type
                tmp += "%10.5f " %float(line_element[3])
                tmp += "%10.5f\n" %float(line_element[4])
                bondtypes.append(tmp)
            elif line_element[2] == "3":
                tmp =  "%10s 3 " %current_bond.type
                tmp += "%10.5f " %float(line_element[3])
                tmp += "%10.5f\n" %float(line_element[4])
                tmp += "%10.5f\n" %float(line_element[5])
                bondtypes.append(tmp)
            else:
                helpers.warning("Bond type in %s is currently not supported!"%line)
    return nbond, bondtypes

#==================================================================================================
# md2xml.read_molecule_angle
#==================================================================================================
def read_molecule_angle(mol,name_list,angle,angletypes):
    """
    This routine acts completely similar as read_bonds
    Currently supported potential is only harmonic.

    """
    from ..md_xml import molecule
    from ..md_xml import helpers
    from ..md_xml import angles
    
    nangle = 0   
    for line in angle:
        line_element = line.split()

        name_i = name_list[line_element[0]]
        name_j = name_list[line_element[1]]
        name_k = name_list[line_element[2]]
        atomnames = angles.sequence(' '.join([name_i,name_j,name_k]))

        atom_i = molecule.atomClass(mol,name_i)
        atom_j = molecule.atomClass(mol,name_j)
        atom_k = molecule.atomClass(mol,name_k)
        current_angle = molecule.angleClass(mol,atomnames)
        current_angle.type = angles.sequence(' '.join([atom_i.type,atom_j.type,atom_k.type]))
        nangle += 1
        
        # if angle potential are defined here, they will be added to angletypes
        # and handled later
        if len(line_element)>4:
            if line_element[3]=="1" or line_element[3]=="2":
                tmp =  "%15s %4s "%(current_angle.type,line_element[3])
                tmp += "%10.5f "  %float(line_element[4])
                tmp += "%10.5f\n" %float(line_element[5])
                angletypes.append(tmp)
            elif line_element[3] == "5":
                tmp =  "%15s 5 " %current_angle.type
                tmp += "%10.5f " %float(line_element[4])
                tmp += "%10.5f " %float(line_element[5])
                tmp += "%10.5f " %float(line_element[6])
                tmp += "%10.5f\n" %float(line_element[7])
                angletypes.append(tmp)
            else:
                helpers.warning("Angle type in %s is currently not supported!"%line)
    return nangle, angletypes

#==================================================================================================
# md2xml.read_dihedral_improper
#==================================================================================================
def read_molecule_dihedral_improper(mol,name_list,dihedral,dihedraltypes):
    """
    Dihedrals and impropers are defined with the same keyword. The difference is the functional form
    Therefore, both are read together. Again the same proceeding as for bonds and angles.
    Currently supported potentials:
    - multiple cosine
    - Ryckaert-Bellmanns
    - Fourier dihedral -> is converted to Ryckaert-Bellmanns
    - Harmonic (only for impropers)

    """
    from ..md_xml import molecule
    from ..md_xml import dihedrals
    from ..md_xml import impropers
    from ..md_xml import helpers
                          
    ndihedral = 0
    nimproper = 0
    for line in dihedral:
        line_element = line.split()
        
        name_i = name_list[line_element[0]]
        name_j = name_list[line_element[1]]
        name_k = name_list[line_element[2]]
        name_l = name_list[line_element[3]]
 
        atom_i = molecule.atomClass(mol,name_i)
        atom_j = molecule.atomClass(mol,name_j)
        atom_k = molecule.atomClass(mol,name_k)
        atom_l = molecule.atomClass(mol,name_l)

        if line_element[4]=='1' or line_element[4]=='3' or line_element[4]=='5':
            ndihedral += 1
            atomnames = dihedrals.sequence(' '.join([name_i,name_j,name_k,name_l]))
            current_dihedral = molecule.dihedralClass(mol,atomnames)
            current_dihedral.type = dihedrals.sequence(' '.join([atom_i.type,atom_j.type,atom_k.type,atom_l.type]))
        elif line_element[4]=='2' or line_element[4]=='4':
            nimproper += 1
            atomnames = impropers.sequence(' '.join([name_i,name_j,name_k,name_l]))
            current_improper = molecule.improperClass(mol,atomnames)
            current_improper.type = impropers.sequence(' '.join([atom_i.type,atom_j.type,atom_k.type,atom_l.type]))
            current_improper.central = name_i
        else:
            helpers.warning("Dihedral type in %s is currently not supported!"%line)
  
        if len(line_element)>5:
            if line_element[4]=='1':
                tmp  = "%20s 1" %current_dihedral.type
                tmp += "%10.5f "%float(line_element[5])
                tmp += "%10.5f "%float(line_element[6])
                tmp += "%10.5f\n"%int(line_element[7])
            elif line_element[4]=='2':
                tmp  = "%20s 2" %current_improper.type
                tmp += "%10.5f "%float(line_element[5])
                tmp += "%10.5f\n"%float(line_element[6])
            elif line_element[4]=='3':
                tmp  = "%20s 3" %current_dihedral.type
                tmp += "%10.5f "%float(line_element[5])
                tmp += "%10.5f "%float(line_element[6])
                tmp += "%10.5f "%float(line_element[7])
                tmp += "%10.5f "%float(line_element[8])
                tmp += "%10.5f "%float(line_element[9])
                tmp += "%10.5f\n"%float(line_element[10])
            elif line_element[4]=='4':
                tmp  = "%20s 4" %current_improper.type
                tmp += "%10.5f "%float(line_element[5])
                tmp += "%10.5f "%float(line_element[6])
                tmp += "%10.5f\n"%int(line_element[7])
            elif line_element[5]=='5':
                tmp  = "%20s 5" %current_dihedral.type
                tmp += "%10.5f "%float(line_element[5])
                tmp += "%10.5f "%float(line_element[6])
                tmp += "%10.5f "%float(line_element[7])
                tmp += "%10.5f\n"%float(line_element[8])
            dihedraltypes.append(tmp)
    return ndihedral, nimproper,dihedraltypes

#==================================================================================================
# md2xml.read_bonds
#==================================================================================================
def read_bonds(root,bondtypes):
    """
    This routine uses bondtypes to create the <bonds> section.
    """
    from ..md_xml import bonds
    from ..md_xml import helpers
    import parameter

    for bond_xml in root.findall('molecule/bond'):
        type_xml = bonds.sequence(bond_xml.get('type'))
        for line in bondtypes:
            if line[0:1]=="#":
                continue
            line_element = line.split()
            type_i = line_element[0].upper()
            type_j = line_element[1].upper()
            type_ff = bonds.sequence(' '.join([type_i,type_j]))

            if parameter.match_type(type_ff.split(),type_xml.split()):
                if line_element[2]=='1':
                    current_harm = bonds.harmClass(root,type_xml)
                    new_r0 = float(line_element[3])
                    new_k  = float(line_element[4])*0.5
                    duplicate = current_harm.equal(new_k,new_r0)
                    if not duplicate and current_harm.r0 != None:
                        helpers.warning("Bond potential %s has already force field parameters"%type_xml)
                    else:
                        current_harm.k  = new_k
                        current_harm.r0 = new_r0
                elif line_element[2]=='3':
                    current_morse = bonds.morsClass(root,atomtype)
                    new_r0   = float(line_element[3])
                    new_D0   = float(line_element[4])
                    new_beta = float(line_element[5])
                    duplicate = current_harm.equal(new_k,new_r0)
                    if not duplicate and current_harm.r0 != None:
                        helpers.warning("Bond potential %s has already force field parameters"%type_xml)
                    else:
                        current_morse.r0   = new_r0
                        current_morse.D0   = new_D0
                        current_morse.beta = new_beta
                else:
                    helpers.warning("Bond potential for %s is not supported by ForConX!"%type_xml)
    current_bonds = bonds.bondsElement(root)
    return len(current_bonds.list("HARM MORS"))

#==================================================================================================
# md2xml.read_angles
#==================================================================================================
def read_angles(root,angletypes):
    """
    This routine uses angletypes to create the <angles> section.
    """
    from ..md_xml import angles
    from ..md_xml import helpers
    import parameter

    for angle_xml in root.findall('molecule/angle'):
        type_xml = angles.sequence(angle_xml.get('type'))
        for line in angletypes:
            if line[0:1]=="#":
                continue
            line_element = line.split()
            type_i = line_element[0].upper()
            type_j = line_element[1].upper()
            type_k = line_element[2].upper()
            type_ff = angles.sequence(' '.join([type_i,type_j,type_k]))

            if parameter.match_type(type_ff.split(),type_xml.split()):
                if line_element[3]=='1':
                    current_harm = angles.harmClass(root,type_xml)
                    new_theta0   = float(line_element[4])
                    new_k        = float(line_element[5])*0.500000000
                    duplicate = current_harm.equal(new_k,new_theta0)
                    if not duplicate and current_harm.theta0 != None:
                        helpers.warning("Angle potential %s has already force field parameters"%type_xml)
                        sys.exit()
                    else:
                        current_harm.theta0 = new_theta0
                        current_harm.k      = new_k
                elif line_element[3]=='2':
                    import numpy
                    helpers.warning("Harmonic cosine of %s converted to harmonic angle"%type_xml)
                    current_harm = angles.harmClass(root,type_xml)
                    # theta0 degree -> radian
                    new_theta0 = float(line_element[4])*0.017453293
                    new_k      =  float(line_element[5])*numpy.sin(new_theta0)*numpy.sin(new_theta0)
                    duplicate = current_harm.equal(new_k,new_theta0)
                    if not duplicate and current_harm.theta0 != None:
                        helpers.warning("Angle potential %s has already force field parameters"%type_xml)
                    else:
                        current_harm.theta0 = new_theta0
                        current_harm.k      = new_k
                elif line_element[3]=='5':
                    current_harm        = angles.harmClass(root,type_xml)
                    current_harm.theta0 = float(line_element[4])
                    current_harm.k      = float(line_element[5])
                    
                    current_urey = angles.ureyClass(root,type_xml)
                    current_urey.r0 = float(line_element[6])
                    current_urey.k  = float(line_element[7])
                else:
                    helpers.warning("Angle potential for %s is not supported by ForConX!"%type_xml)
    current_angles = angles.anglesElement(root)
    return len(current_angles.list("HARM UREY"))

#==================================================================================================
# md2xml.read_dihedrals
#==================================================================================================
def read_dihedrals(root,dihedraltypes):
    """
    """
    from ..md_xml import dihedrals
    from ..md_xml import helpers
    import parameter

    for dihedral_xml in root.findall('molecule/dihedral'):
        type_xml = dihedrals.sequence(dihedral_xml.get('type'))
        for line in dihedraltypes:
            if line[0:1]=="#":
                continue
            line_element = line.split()
            type_i = line_element[0].upper()
            type_j = line_element[1].upper()
            type_k = line_element[2].upper()
            type_l = line_element[3].upper()
            type_ff = dihedrals.sequence(' '.join([type_i,type_j,type_k,type_l]))

            if parameter.match_type(type_ff.split(),type_xml.split()):
                if line_element[4]=='1':
                    cos = dihedrals.cosClass(root,type_xml)
                    old_n = cos.n
                    old_k = cos.k
                    old_delta = cos.delta

                    current_n = float(line_element[7])
                    current_k = float(line_element[6])
                    current_delta = float(line_element[5])

                    print old_n, current_n, old_n==current_n
                    
                    if old_n != current_n or old_k != current_k:
                        
                        cos.n = old_n + [current_n]
                        cos.k = old_k + [current_k]
                        cos.delta = old_delta + [current_delta]
                elif line_element[4]=='3':
                    ryck = dihedrals.ryckClass(root,type_xml)
                    k = []
                    for i in range(5,11):
                        k.append(float(line_element[i]))
                    ryck.k = k
                else:
                    helpers.warning("Dihedral potential for %s is not supported by ForConX!"%type_xml)
    current_dihedrals = dihedrals.dihedralsElement(root)
    return len(current_dihedrals.list("COS RYCK"))

#==================================================================================================
# md2xml.read_impropers
#==================================================================================================
def read_impropers(root,dihedraltypes):
    """
    """
    from ..md_xml import impropers
    from ..md_xml import helpers
    import parameter

    for improper_xml in root.findall('molecule/improper'):
        type_xml = impropers.sequence(improper_xml.get('type'))
        for line in dihedraltypes:
            if line[0:1]=="#":
                continue
            line_element = line.split()
            type_i = line_element[0].upper()
            type_j = line_element[1].upper()
            type_k = line_element[2].upper()
            type_l = line_element[3].upper()
            type_ff = impropers.sequence(' '.join([type_i,type_j,type_k,type_l]))

            if parameter.match_type(type_ff.split(),type_xml.split()):
                if line_element[4]=='2':
                    harm = impropers.harmClass(root,type_xml)
                    harm.theta0   = float(line_element[4])
                    harm.k        = float(line_element[5])
                elif line_element[4]=='4':
                    cos = impropers.cosClass(root,type_xml)
                    old_n = cos.n
                    old_k = cos.k
                    old_delta = cos.delta

                    current_n = float(line_element[7])
                    current_k = float(line_element[6])
                    current_delta = float(line_element[5])
                    cos.n = old_n + [current_n]
                    cos.k = old_k + [current_k]
                    cos.delta = old_delta + [current_delta]
                else:
                    helpers.warning("Improper potential for %s is not supported by ForConX!"%type_xml)
    current_impropers = impropers.impropersElement(root)
    return len(current_impropers.list("COS HARM"))
    
#==================================================================================================
# md2xml.van_der_waals
#==================================================================================================
def van_der_waals(root,parameterlist,fudgelj,fudgeqq,combination_rule):
    """
    """
    from ..md_xml import molecule
    from ..md_xml import nonbonded
    from ..md_xml import bonds
    from ..md_xml import helpers
    import math
    
    inb = 0
    opls = 0
    flag = None
    lj = {}
    for line in parameterlist:
        if "[ FF_OPLS ]" in line:
            opls = 1
        if "atomtypes" in line:
            flag = "ATOMTYPES"
            continue
        if "moleculetype" in line:
            flag = "PAIRS"
            continue
        if "pairtypes" in line:
            flag = "PAIRS"
            continue
        if flag == "ATOMTYPES":
            line_element = line.split()
            if line[0:1]=="#":
                continue
            type = line_element[0+opls].upper()
            lj[type] = [float(line_element[6+opls]),float(line_element[5+opls])]
                
        if flag == "PAIRS":
            # @Marcello feel free to add content
            pass

    current_nonbonded = nonbonded.nonbondedElement(root)
    type_list = []
    for mol in root.findall("molecule"):
        current_molecule = molecule.moleculeElement(mol,mol.get('name'))
        for atom in current_molecule.list("ATOM"):
            current_atom = molecule.atomClass(mol,atom)
            type_list.append(current_atom.type)
    type_list = sorted(set(type_list))
        
    # combination rules
    if combination_rule==2:
        current_nonbonded.mixing_epsilon="geometric"
        current_nonbonded.mixing_sigma="arithmetic"
        for atomtype in type_list:
            current_atom = nonbonded.atomClass(root,atomtype)
            try:
                current_atom.epsilon = lj[atomtype][0]
                current_atom.sigma = lj[atomtype][1]
                current_atom.elec14 = fudgeqq
                current_atom.vdw14 = fudgelj
                inb += 1
            except:
                helpers.warning("Lennard-Jones parameters for %s not found!"%atomtype)
    elif (combination_rule==1 or combination_rule==3):
        # combination rule 1 or 3
        current_nonbonded.mixing_epsilon="geometric"
        current_nonbonded.mixing_sigma="geometric"
        
        for type_i in type_list:
            try:
                sigma_i   = lj[type_i][1]
                epsilon_i = lj[type_i][0]
            except:
                helpers.warning("Lennard-Jones parameters for %s not found!"%type_i)
                continue
                
            for type_j in type_list:
                try:
                    sigma_j   = lj[type_j][1]
                    epsilon_j = lj[type_j][0]
                except:
                    helpers.warning("Lennard-Jones parameters for %s not found!"%type_i)
                    continue
                atomtype = bonds.sequence(' '.join([type_i,type_j]))
                current_vdw = nonbonded.vdwClass(root,atomtype)
                current_vdw.sigma = math.sqrt(sigma_i*sigma_j)
                current_vdw.epsilon = math.sqrt(epsilon_i*epsilon_j)
                current_vdw.elec14 = fudgeqq
                current_vdw.vdw14 = fudgelj
                inb += 1
    else:
        helpers.error("Wrong mixing rule: %d "%combination_rule)
        sys.exit()
    return inb

