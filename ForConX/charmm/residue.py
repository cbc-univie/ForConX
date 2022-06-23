#======================================================================================
# residue.read_atom
#======================================================================================
def read_atom(mol,atom,mass):
    """
    reads information concerning the atoms of the CHARMM topology from string lines
    and converts them to XML.
    """
    from ..md_xml import molecule
    
    iatom = 0
    ivirtual = 0

#   processing information on atoms / virtual particles    
    for line in atom:
        line_element = line.split()
        name = line_element[1]
        type = line_element[2]
        charge = float(line_element[3])
        if mass[type] > 0.1:
            iatom += 1
            current_atom = molecule.atomClass(mol,name)
            current_atom.mass = mass[type]
            alpha_pos = line.find('ALPHA')
            if alpha_pos>0:
                tmp = line[alpha_pos:].split()
                current_atom.alpha = abs(float(tmp[1]))
        else:
            ivirtual += 1
            current_atom = molecule.virtualClass(mol,name)
        current_atom.type   = type
        current_atom.charge = charge
    return iatom, ivirtual

#======================================================================================
# residue.write_atom
#======================================================================================
def write_atom(current_molecule,polarizability_unit):
    """
    reads information concerning the atoms of the CHARMM topology from XML
    and converts them to string lines.
    """
    from ..md_xml import molecule

    atom = []
    for i in current_molecule.list("ATOM VIRTUAL"):
        (atomclass,pointer) = current_molecule.find_name(i)
        mol = current_molecule.mol
        if atomclass=="ATOM":
            current_atom = molecule.atomClass(mol,i)
            line  = "ATOM "
            line += "%6s "  % i
            line += "%6s "  % current_atom.type
            line += "%10s " % current_atom.charge
            alpha = current_atom.alpha
            if alpha is not None:
                alpha = -alpha * polarizability_unit
                line += " ALPHA %10.5f THOLE 0.0" % alpha
            line += ("\n")
        if atomclass=="VIRTUAL":
            current_atom = molecule.virtualClass(mol,i)
            line  = "ATOM "
            line += "%6s "  % i
            line += "%6s "  % current_atom.type
            line += "%10s " % current_atom.charge
            line += ("\n")
        atom.append(line)
    return atom

#======================================================================================
# residue.read_bond
#======================================================================================
def read_bond(root,mol,bond):
    """
    sets bond information in <molecule> and add corresponding bonds to <bonds>.
    """
    from ..md_xml import molecule
    from ..md_xml import bonds
    from ..md_xml import helpers
    ibond = 0

    current_molecule = molecule.moleculeElement(mol,mol.get('name'))
    atomlist    = current_molecule.list('atom')
    virtuallist = current_molecule.list('virtual')
    
    for line in bond:
        line_element = line.split()
#       Looking for several bonds in corresponding line
        del line_element[0]
        nbond = int(len(line_element)/2)
        for index in range(nbond):
            name_i = remove_prevnext(line_element[2*index])
            name_j = remove_prevnext(line_element[2*index+1])

            if (name_i in atomlist) and (name_j in atomlist):
                molecule.bondClass(mol,' '.join([name_i,name_j]))
            
                atom_i = molecule.atomClass(mol,name_i)
                atom_j = molecule.atomClass(mol,name_j)
                bonds.harmClass(root,' '.join([atom_i.type,atom_j.type]))
                ibond += 1
            elif (name_i in virtuallist) or (name_j in virtuallist):
                helpers.warning('Bond %s includes virtual atoms'%' '.join([name_i,name_j]))
    return ibond

#======================================================================================
# residue.write_bond
#======================================================================================
def write_bond(current_molecule):
    """
    """
    from ..md_xml import molecule

    bond = []
    mol = current_molecule.mol
    for i in current_molecule.list('bond'):
        current_bond = molecule.bondClass(mol,i)
        name = (current_bond.name).split()
        line  = "BOND "
        line += "%6s "  %name[0]
        line += "%6s \n"%name[1]
        bond.append(line)
    return bond

#======================================================================================
# residue.read_angle
#======================================================================================
def read_angle(root,mol,angle):
    """
    sets angle information in <molecule> and add corresponding angles to <angles>.
    """
    from ..md_xml import molecule
    from ..md_xml import angles
    
    iangle = 0
    current_molecule = molecule.moleculeElement(mol,mol.get('name'))
    #   AUTOGENERATION
    if "AUTOGENERATE" in angle:
        for name in current_molecule.generate_angle():
            molecule.angleClass(mol,name)
            atoms = name.split()
            atom_i = molecule.atomClass(mol,atoms[0])
            atom_j = molecule.atomClass(mol,atoms[1])
            atom_k = molecule.atomClass(mol,atoms[2])
            angles.harmClass(root,' '.join([atom_i.type,atom_j.type,atom_k.type]))
            iangle += 1
        return iangle
        
#   Using angle information from topology file    
    for line in angle:
        line_element = line.split()
#       Looking for several angles in line ANGL ...        
        del line_element[0]
        nangle = int(len(line_element)/3)
        for index in range(nangle):
            name_i = remove_prevnext(line_element[3*index])
            name_j = remove_prevnext(line_element[3*index+1])
            name_k = remove_prevnext(line_element[3*index+2])
            molecule.angleClass(mol,' '.join([name_i,name_j,name_k]))
            
            atom_i = molecule.atomClass(mol,name_i)
            atom_j = molecule.atomClass(mol,name_j)
            atom_k = molecule.atomClass(mol,name_k)
            angles.harmClass(root,' '.join([atom_i.type,atom_j.type,atom_k.type]))
            iangle += 1
    return iangle
            
#======================================================================================
# residue.write_angle
#======================================================================================
def write_angle(current_molecule):
    """
    """
    from ..md_xml import molecule

    angle = []
    mol = current_molecule.mol
    for i in current_molecule.list('angle'):
        current_angle = molecule.angleClass(mol,i)
        name = (current_angle.name).split()
        line  = "ANGL "
        line += "%6s "  %name[0]
        line += "%6s "  %name[1]
        line += "%6s\n" %name[2]
        angle.append(line)
    return angle

#======================================================================================
# residue.read_dihedral
#======================================================================================
def read_dihedral(root,mol,dihedral,ic):
    """
    sets dihedral information in <molecule> and add corresponding dihedrals to <dihedrals>.
    """
    from ..md_xml import molecule
    from ..md_xml import dihedrals
    from ..md_xml import helpers
    
    idihedral = 0
    dangerous_angles = []
    current_molecule= molecule.moleculeElement(mol,mol.get('name'))
    
    if "AUTOGENERATE" in dihedral:
#       Looking for dangerous angles 
        for line in ic:
            line_element = line.split()
            if (float(line_element[6])>179.):
                dangerous_angles.append([line_element[1], line_element[2], line_element[3]])
            if (float(line_element[8])>179.):
                dangerous_angles.append([line_element[2], line_element[3], line_element[4]])
#       Constructing all possible dihedrals
        for name in current_molecule.generate_dihedral():
            atoms = name.split()
            atom_i = molecule.atomClass(mol,atoms[0])
            atom_j = molecule.atomClass(mol,atoms[1])
            atom_k = molecule.atomClass(mol,atoms[2])
            atom_l = molecule.atomClass(mol,atoms[3])
            angle_ijk = [atoms[0],atoms[1],atoms[2]]
            angle_jkl = [atoms[1],atoms[2],atoms[3]]
            dangerous = False
            for line in dangerous_angles:
                if line==angle_ijk or line==angle_ijk[::-1]:
                    dangerous = True
                if line==angle_jkl or line==angle_jkl[::-1]:
                    dangerous = True

            idihedral += 1
            current_dihedral = molecule.dihedralClass(mol,name)
            current_cos = dihedrals.cosClass(root,' '.join([atom_i.type,atom_j.type,atom_k.type,atom_l.type]))
            if dangerous:
                helpers.warning('<dihedral name="%s">  contains an angle close to 180 deg. Force constants set to zero.'%atoms)
                current_dihedral.linear = "yes"
                current_cos.n = [3]
                current_cos.k = [0]
                current_cos.delta = [0]
        return idihedral
    
#   Using given DIHEDRAL information from topology file  
    for line in dihedral:
        line_element = line.split()
#       Looking for several dihedrals in line DIHE
        del line_element[0]
        ndihedral = int(len(line_element)/4)
        for index in range(ndihedral):
            name_i = remove_prevnext(line_element[4*index])
            name_j = remove_prevnext(line_element[4*index+1])
            name_k = remove_prevnext(line_element[4*index+2])
            name_l = remove_prevnext(line_element[4*index+3])
            molecule.dihedralClass(mol,' '.join([name_i, name_j, name_k, name_l]))

            atom_i = molecule.atomClass(mol,name_i)
            atom_j = molecule.atomClass(mol,name_j)
            atom_k = molecule.atomClass(mol,name_k)
            atom_l = molecule.atomClass(mol,name_l)
            dihedrals.cosClass(root,' '.join([atom_i.type,atom_j.type,atom_k.type,atom_l.type]))
            idihedral += 1
    return idihedral

#======================================================================================
# residue.write_dihedral
#======================================================================================
def write_dihedral(current_molecule):
    """
    """
    from ..md_xml import molecule

    mol = current_molecule.mol
    dihedral = []
    for i in current_molecule.list('dihedral'):
        current_dihedral = molecule.dihedralClass(mol,i)
        name = (current_dihedral.name).split()
        line  = "DIHE "
        line += "%6s "  %name[0]
        line += "%6s "  %name[1]
        line += "%6s "  %name[2]
        line += "%6s\n" %name[3]
        dihedral.append(line)
    return dihedral

#======================================================================================
# residue.read_improper
#======================================================================================
def read_improper(root,mol,improper):
    """
    sets improper information in <molecule> and add corresponding impropers to <impropers>.
    The potential type of the improper is temporarily set to cosine but will be changed,
    if defined differently in the CHARMM parameter file.
    """
    from ..md_xml import molecule
    from ..md_xml import impropers
    
    iimproper = 0
    current_molecule = molecule.moleculeElement(mol,mol.get('name'))
    for line in improper:
        line_element = line.split()
#       Looking for several impropers in line IMPR ...  
        del line_element[0]
        nimproper = int(len(line_element)/4)
        for index in range(nimproper):
            name_i = remove_prevnext(line_element[4*index])
            name_j = remove_prevnext(line_element[4*index+1])
            name_k = remove_prevnext(line_element[4*index+2])
            name_l = remove_prevnext(line_element[4*index+3])
            improper = molecule.improperClass(mol,' '.join([name_i, name_j, name_k, name_l]))
            improper.central = name_i
            
            atom_i = molecule.atomClass(mol,name_i)
            atom_j = molecule.atomClass(mol,name_j)
            atom_k = molecule.atomClass(mol,name_k)
            atom_l = molecule.atomClass(mol,name_l)
            cos = impropers.cosClass(root,' '.join([atom_i.type,atom_j.type,atom_k.type,atom_l.type]))
            iimproper += 1
    return iimproper

#======================================================================================
# residue.write_improper
#======================================================================================
def write_improper(current_molecule):
    """
    """
    from ..md_xml import molecule

    improper = []
    mol = current_molecule.mol
    for i in current_molecule.list('improper'):
        current_improper = molecule.improperClass(mol,i)
        name = (current_improper.name).split()
        line  = "IMPR "
        line += "%6s "  %name[0]
        line += "%6s "  %name[1]
        line += "%6s "  %name[2]
        line += "%6s\n" %name[3]
        improper.append(line)
    return improper

#======================================================================================
# residue.read_lonepair
#======================================================================================
def read_lonepair(mol,lonepair):
    """

    -> molecule.virtualClass
    """
    from ..md_xml import molecule

    for line in lonepair:
        line_element = (line.upper()).split()
        name_i = line_element[2]
        name_j = line_element[3]
        name_k = line_element[4]
        name_l = line_element[5]
        dist   = float(line_element[7])
        theta = float(line_element[9])
        phi   = float(line_element[11])

        virtual = molecule.virtualClass(mol,name_i)
        if line_element[1]=="RELATIVE":
            virtual.zmatrix = ' '.join([name_j,name_k,name_l])
            virtual.dist = dist
            virtual.theta = theta
            virtual.phi = phi
        if line_element[1]=="BISECTOR":
#
#                        T O     D O
#
            pass
        
#======================================================================================
# residue.read_lonepair
#======================================================================================
def write_lonepair(current_molecule,distance_unit):
    """

    -> molecule.virtualClass
    """
    from ..md_xml import molecule
    
    lonepair = []
    mol = current_molecule.mol
    for i in current_molecule.list("VIRTUAL"):
        current_virtual = molecule.virtualClass(mol,i)
        if  current_virtual.zmatrix is not None:
            line  = "LONEPAIR RELATIVE"
            name  = current_virtual.zmatrix
            dist  = current_virtual.dist * distance_unit
            theta = current_virtual.theta
            phi   = current_virtual.phi
            line += " %s %s "         % (i,name)
            line += "distance %6.3f " % dist
            line += "angle %6.1f "    % theta
            line += "dihe %6.1f\n"    % phi
            lonepair.append(line)
    return lonepair

#======================================================================================
# residue.remove_prevnext
#======================================================================================
def remove_prevnext(name):
    """
    removes link between this residue and next/previous residue:
    + indicates interresidue potential to next residue
    - to previous residue.

    So far, ForConX is not able to take this information into account. Consequently,
    the sequence of amino acids is splitted.
    """
    if name[0]=="+" or name[0]=="-":
        name = name[1:]
    return name
