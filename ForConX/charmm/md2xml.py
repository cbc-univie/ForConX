import xml.etree.ElementTree as ET
#======================================================================================
# md2xml.read
#======================================================================================
def read(root):
    """
    converts CHARMM topology and parameter file into XML.
    
    -> read_topology()    
    -> residue2xml()      
    -> read_parameter()   
    -> parameter2xml()    
    """
    from ..md_xml import helpers

    print"\n\t----------------------------------"
    print"\t1.1 Reading topology"
    print"\t----------------------------------"
    current_input = root.find('input/topology')
    try:
        topology_file = current_input.get('file') 
        f = open(topology_file)
        print '> Topology file = %s ...'%topology_file
    except IOError:
        helpers.error('<input> Topology file %s cannot be opened!'%topology_file)
    except:
        helpers.error('<input/topology file="..."> is ill-defined!')
        
    (mass,residue,autogenerate_angle,autogenerate_dihedral) = read_topology(root,f)
    f.close()
    
    print"\n\t----------------------------------"
    print"\t1.2 Converting topology to XML"
    print"\t----------------------------------"
    moleculelist = []
    for mol in root.findall('molecule'):
        molname = mol.get('name')
        if not molname in residue:
            helpers.error('<molecule name="%s"> cannot be found in topology file'%molname)
        else:
            # Clear existing information        
            for i in mol.findall('*'):
                mol.remove(i)
            
        print'\t<molecule name="%s">'%molname
        moleculelist.append(molname)
        residue2xml(root,mol,mass,residue[molname],
                    autogenerate_angle,autogenerate_dihedral)

    print"\n\t----------------------------------"
    print"\t1.3 Reading parameter"
    print"\t----------------------------------"
    current_input = root.find('input/parameter')
    try:
        parameter_file = current_input.get('file') 
        f = open(parameter_file)
        print '> Parameter file = %s ...'%parameter_file
    except IOError:
        helpers.error('<input> Parameter file %s cannot be opened!'%parameter_file)
    except:
        helpers.error('<input/parameter file="..."> is ill-defined!')
        
    (parameter,elec14) = read_parameter(root,f)
    f.close()

    print"\n\t-------------------------------------------------"
    print"\t1.4 Converting parameter to XML"
    print"\t-------------------------------------------------"
    parameter2xml(root,parameter,elec14)

#======================================================================================
# md2xml.read_topology
#======================================================================================
def read_topology(root,f):
    """
    reads CHARMM topology file and stores the information in the dictionaries
    mass (key=type) and residue(key=molname).
    """
    # reading mass section
    autogenerate_angle    = False
    autogenerate_dihedral = False
    mass = {}
    while True:
        # mass section ends with the first residue        
        rewind_position = f.tell()
        line = readline(f)
        if line[0:4]=='RESI':
            f.seek(rewind_position)
            break
        
        # creating mass dictionary
        elif line[0:4]=='MASS':
            line_element = line.split()
            mass[line_element[2]] = float(line_element[3])
            
        # checking for autogeneration
        elif line[0:4]=='AUTO':
            if 'ANGL' in line[5:]:
                autogenerate_angle = True
            if 'DIHE' in line[5:]:
                autogenerate_dihedral = True

        # anything else?
        else:
            print "?\t\t",line
                
    # Reading molecules            
    residue = {}
    while True:
        line = readline(f)
        # stops at end of topology file        
        if line[0:3]=='END':
            break

        # new residue        
        if line[0:4]=='RESI':
            name = (line.split())[1]
            if not name in residue:
                residue[name] = []
        residue[name].append(line.strip())
    return mass, residue, autogenerate_angle, autogenerate_dihedral

#======================================================================================
# md2xml.residue2xml
#======================================================================================
def residue2xml(root,mol,mass,residue,autogenerate_angle,autogenerate_dihedral):
    """
    decomposes the dictionary of one residue(key=molname) into atom, bond, angle, dihedral,
    improper and ic string arrays. These arrays are handle by the corresponding functions:

    -> residue_atom.read_atom()
    -> residue_bond.read_bond()
    -> residue_angle.read_angle()
    -> residue_dihedral.read_dihedral()
    -> residue_improper.read_improper()
    """
    
    atom = []
    bond = []
    angle = []
    dihedral = []
    improper = []
    lonepair = []
    ic = []
    for line in residue:
        if line[0:4] == "ATOM":
            atom.append(line)
        elif line[0:4] == "BOND": 
            bond.append(line)
        elif line[0:4] == "DOUB":
            bond.append(line)
        elif line[0:4] == "TRIP":
            bond.append(line)
        elif line[0:4] == "AROM":
            bond.append(line)
        elif line[0:4] == "ANGL":
            angle.append(line)
        elif line[0:4] == "THET":
            angle.append(line)
        elif line[0:4] == "DIHE":
            dihedral.append(line)
        elif line[0:4] == "IMPH":
            improper.append(line)
        elif line[0:4] == "IMPR":
            improper.append(line)
        elif line[0:8] == "LONEPAIR":
            lonepair.append(line)
        elif (line[0:2] == "IC"):
            ic.append(line)
        elif line[0:4] == "RESI":
            continue
        elif line[0:5] == "GROUP":
            continue
        else:
            print "?\t\t",line
    print
            
    import residue
    from ..md_xml import helpers

    # Atoms and virtual atoms
    (natom,nvirtual) = residue.read_atom(mol,atom,mass)
    residue.read_lonepair(mol,lonepair)    
    print'\t\tNumber of atoms         = %5s'%natom
    print'\t\tNumber of virtual atoms = %5s'%nvirtual
    
    # Bonds
    nbond = residue.read_bond(root,mol,bond)
    print'\t\tNumber of bonds         = %5s'%nbond

    # Angles
    if len(angle)==0:
        if not autogenerate_angle:
            helpers.warning("Molecule has no angles at all!")
        else:
            angle = ["AUTOGENERATE"]  
    nangle = residue.read_angle(root,mol,angle)
    print'\t\tNumber of angles        = %5s'%nangle,
    if autogenerate_angle:
        print' (autogenerated)',
    print

    # Dihedrals
    if len(dihedral)==0:
        if not autogenerate_dihedral:
            helpers.warning("Molecule has no dihedrals at all!")
        else:
            dihedral = ["AUTOGENERATE"]
    ndihedral = residue.read_dihedral(root,mol,dihedral,ic)
    print'\t\tNumber of dihedrals     = %5s'%ndihedral,
    if autogenerate_dihedral:
        print' (autogenerated)',
    print

    # Impropers
    nimproper = residue.read_improper(root,mol,improper)
    print'\t\tNumber of impropers     = %5s\n'%nimproper

#======================================================================================
# md2xml.read_parameter
#======================================================================================
def read_parameter(root,f):
    """
    decomposes parameter information into the following dictionary parameter with the
    following keys:
    * BONDS
    * ANGLES
    * DIHEDRALS
    * IMPROPERS
    * NONBONDED
    * NBFIX
    * DUMP    
    """
    from ..md_xml import helpers
    
    parameter = {}
    parameter["BONDS"] = []
    parameter["ANGLES"] = []
    parameter["DIHEDRALS"] = []
    parameter["IMPROPERS"] = []
    parameter["NONBONDED"] = []
    parameter["NBFIX"] = []
    parameter["DUMP"] = []
    key = "DUMP"
    rtf = True

    while True:
        line = readline(f)
        
        if len(parameter["BONDS"])>0:
            rtf = False
        if line[0:3]=="END" and not rtf:
            break    
        if line[0:5] == "BONDS":
            key = "BONDS"
            continue
        if line[0:6] == "ANGLES":
            key = "ANGLES"
            continue
        if line[0:5] == "DIHED":
            key = "DIHEDRALS"
            continue
        if line[0:5] == "IMPRO":
            key = "IMPROPERS"
            continue
        if line[0:5] == "NBFIX":
            key = "NBFIX"
            continue
        if line[0:5] == "NONBO":
            key = "NONBONDED"
            print "?\t",line
            if "SHIFT" in line:
                helpers.warning('SHIFT option will be disregarded due to transferability. This will result in a change of vdW energy!')
            line_element = line.split()
            try:
                elec14 = float(line_element[line_element.index("E14FAC")+1])
            except:
                elec14 = 1.0
            print "\n\t\telec14 = ",elec14
        if line[0:5]=="HBOND":
            key = "DUMP"
            print "?\t",line
        if line[0:5]=="THOLE":
           key = "DUMP"
           print "?\t",line
        parameter[key].append(line)

#    All force field lines containing wild cards are shifted to the end of the corresponding
#    parameter array, since CHARMM overwrites wild card containing potentials with respective
#    potentials containing particular atomtypes. Thus, wild card potentials are only taken into
#    account if no "simple" potential exists.    
    for key in parameter:
        new_parameter = []
        for line in parameter[key]:
            line_element = line.split()
            if "X" in line_element:
                new_parameter.insert(0,line)
            else:
                new_parameter.append(line)        
        parameter[key] = new_parameter[::-1]
    return parameter, elec14

#======================================================================================
# md2xml.parameter2xml
#======================================================================================
def parameter2xml(root,forcefield,elec14):
    """

    -> parameter.read_bonds()
    -> parameter.read_angles()
    -> parameter.read_dihedrals()
    -> parameter.read_impropers()
    -> parameter.read_nonbonded()
    """
    import parameter
    
    ibonds = parameter.read_bonds(root,forcefield["BONDS"])
    print "\t\tNumber of bond potentials       = %6s" %ibonds

    iangles = parameter.read_angles(root,forcefield["ANGLES"])
    print "\t\tNumber of angles potentials     = %6s" %iangles

    idihedrals = parameter.read_dihedrals(root,forcefield["DIHEDRALS"])
    print "\t\tNumber of dihedral potentials   = %6s" %idihedrals

    iimpropers = parameter.read_impropers(root,forcefield["IMPROPERS"])
    print "\t\tNumber of improper potentials   = %6s" %iimpropers

    (iatom,ivdw) = parameter.read_nonbonded(root,
                                            forcefield["NONBONDED"],
                                            forcefield["NBFIX"],
                                            elec14)
    print "\n\t\tNumber of atomic Lennard-Jones  = %6s" %iatom
    print "\t\tNumber of pair   Lennard-Jones  = %6s" %ivdw

#======================================================================================
# md2xml.readline
#======================================================================================
def readline(f):
    """
    reads one line from CHARMM force field file and
    * converts characters to upper case
    * cuts everything beyond "!" since it is a comment
    * skips blank lines
    """
    while True:
        buffer = ((f.readline()).upper().rstrip())
#       removing comments
        pos = buffer.find('!')
        if pos>0:
            buffer = (buffer[:pos-1]).rstrip()
        elif pos==0:
            continue
        else:
            buffer = buffer.strip()
        if len(buffer)<2:
            continue
#       lines ending with "-" are merged with the next line
        while True:
            if (buffer.strip())[-1]=="-":
                buffer = buffer[:-1].strip() +' '
                buffer += f.readline().strip()
            else:
                break
        break
    return buffer

