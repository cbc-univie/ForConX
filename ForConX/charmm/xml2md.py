import xml.etree.ElementTree as ET
#======================================================================================
# xml2md.write_topology
#======================================================================================
def write_topology(root):
    """
    converts XML structure to CHARMM topology file.

    -> input_output.header
    -> input_output.warranty
    -> molecule.moleculeElement

    -> residue.write_atom
    -> residue.write_bond
    -> residue.write_angle
    -> residue.write_dihedral
    -> residue.write_improper
    -> residue.write_lonepair

    """
    from ..md_xml import molecule
    from ..md_xml import input_output
    from ..md_xml import helpers
    
    pointer = root.find('./output/topology')
    try:
        topology_file = pointer.get('file')
        print"> Writing topology file %s ..."%topology_file
        f = open(topology_file,'w')
    except IOError:
        helpers.error('<output> Topology file cannot be opened.')
    except:
        helpers.error("<output> Can't find topology information.")

    current_output = input_output.mdElement(root,'output')
    distance_conversion = current_output.convert_distance()
    polarizability_conversion = current_output.convert_polarizability()
        
    input_output.header(f,"*.")
    f.write("*.\n")
    input_output.warranty(f,"*.")
    f.write("*"+"."*70+"\n")
    f.write("*\n")
    f.write("%5s %5s\n" %("99","1"))
    f.write(" \n")

    print"\n\t----------------------------------"
    print"\t3.1 Writing mass section"
    print"\t----------------------------------"
    type = []
    for i in root.findall('./molecule/atom'):
        line  = "%-6s  %8.3f" % (i.get("type"),float(i.get("mass")))
        type.append(line)
    for i in root.findall('./molecule/virtual'):
        line  = "%-6s  %8.3f" % (i.get("type"),0.0)
        type.append(line)
    type = sorted(set(type))
    
    polarizable = False
    if root.findall('./molecule/atom/[@alpha]'):
        line  = "%-6s  %8.3f" % ("DRUD",0.0)
        polarizable = True
        type.append(line)

    # generating mass section of topology file
    imass = 0
    mass = []
    for line in type:
        imass +=1
        mass.append("MASS %4s %s\n" % (str(imass),line))
    mass.append("\n")
    if polarizable:
        print "\t\tPolarizable force field ..."
        mass.append("AUTOGENERATE DRUDE\n\n")
    for line in mass:
        f.write(line)
        
    import residue
    print"\n\t----------------------------------"
    print"\t3.2 Writing residues"
    print"\t----------------------------------"
    print"!\tWarning: In the current ForConX version connection between molecules is not supported!"
    print"\t\t Thus, bonds, angles and dihedrals involving atoms of two different molecules,"
    print"\t\t e.g. amino acids, have to be added manually. Sorry for the inconvenience.\n"
    for mol in root.findall("./molecule"):
        molname = mol.get("name")
        print"\t\tRESI ",molname
        current_molecule = molecule.moleculeElement(mol,molname)
        molname = current_molecule.name
        current_residue = []
        current_residue.append("!"+"~"*60+"\n")
        current_residue.append("RESI %4s %8.5f\n" % (molname, current_molecule.charge))
        current_residue.append("!"+"~"*60+"\n")
        current_residue.append("GROUP\n")

        current_residue += residue.write_atom(current_molecule,polarizability_conversion)
        current_residue += residue.write_bond(current_molecule)
        current_residue += residue.write_angle(current_molecule)
        current_residue += residue.write_dihedral(current_molecule)
        current_residue += residue.write_improper(current_molecule)
        current_residue += residue.write_lonepair(current_molecule,distance_conversion)
        # This has to be changed in case of polymers and proteins!        
        current_residue.append("PATCHING FIRST NONE LAST NONE\n\n")
        
        for line in current_residue:
            f.write(line)
    f.write("END\n")
    print('\n\n')

#======================================================================================
# xml2md.write_parameter
#======================================================================================
def write_parameter(root):
    """
    converts XML structure to CHARMM topology file.

    -> input_output.header

    -> parameter_bonds.write_bonds
    -> parameter_angles.write_angles
    -> parameter_dihedrals.write_dihedrals
    -> parameter_impropers.write_impropers
    -> parameter_nonbonded.write_vdws
    -> parameter_nonbonded.write_exclusions

    """
    from ..md_xml import helpers
    pointer = root.find('./output/parameter')
    try:
        parameter_file = pointer.get('file')
        print"> Writing parameter file %s ..."%parameter_file
        f = open(parameter_file,'w')
    except IOError:
        helpers.error('<output> Parameter file cannot be opened.')
    except:
        helpers.error("<output> Can't find parameter information.")
        
    from ..md_xml import input_output
    input_output.conk(f,"*.")
    f.write("*\n\n")

    current_output = input_output.mdElement(root,'output')
    distance_conversion = current_output.convert_distance()
    energy_conversion = current_output.convert_energy()
        
    import parameter
    print"\n\t-----------------------------------------"
    print"\t3.3 Writing intramolecular potentials"
    print"\t-----------------------------------------"
    f.write("BONDS\n")
    f.write("!\n")
    f.write("! U_bond = k ( r - r0 )^2\n")
    f.write("!\n")
    f.write("!"+"~"*70+"\n")
    f.write("!TYPE1    TYPE2    k [kcal/mol Angstroem^2]     r0 [Angstroem]\n")
    f.write("!"+"~"*70+"\n")
    bonds = sorted(parameter.write_bonds(root,energy_conversion,distance_conversion))
    for line in bonds:
        f.write(line)
    f.write("\n\n")
    print"\t\tNumber of bond potentials     = ",len(bonds)

    f.write("ANGLES\n")
    f.write("!\n")
    f.write("! U_angle = k ( theta - theta0 )^2\n")
    f.write("!\n")
    f.write("!"+"~"*70+"\n")
    f.write("!TYPE1    TYPE2     TYPE3     k [kcal/mol rad^2]     theta0 [deg]\n")
    f.write("!"+"~"*70+"\n")
    angles =  sorted(parameter.write_angles(root,energy_conversion,distance_conversion))
    for line in angles:
        f.write(line)
    f.write("\n\n")
    print"\t\tNumber of angle potentials    = ",len(angles)

    f.write("DIHEDRALS\n")
    f.write("!\n")
    f.write("! U_dihedral = k ( 1 + Cos[n phi - delta] )\n")
    f.write("!\n")
    f.write("!"+"~"*70+"\n")
    f.write("!TYPE1    TYPE2     TYPE3     TYPE4    k [kcal/mol]    n     delta [deg]\n")
    f.write("!"+"~"*70+"\n")
    dihedrals = sorted(parameter.write_dihedrals(root,energy_conversion))
    for line in dihedrals:
        f.write(line)
    f.write("\n\n")
    print"\t\tNumber of dihedral potentials = ",len(dihedrals)

    f.write("IMPROPERS\n")
    f.write("! n > 0: \n")
    f.write("! U_improper = k ( 1 + Cos[n phi - delta] )\n")
    f.write("!\n")
    f.write("! n = 0: \n")
    f.write("! U_improper = k ( phi - delta )^2\n")
    f.write("!\n")
    f.write("!"+"~"*80+"\n")
    f.write("!TYPE1    TYPE2     TYPE3     TYPE4    k [kcal/mol]    n     delta [deg]\n")
    f.write("!"+"~"*80+"\n")
    impropers = sorted(parameter.write_impropers(root,energy_conversion))
    for line in impropers:
        f.write(line)
    f.write("\n\n")
    print"\t\tNumber of improper potentials = ",len(impropers)

    print"\n\t--------------------------------------"
    print"\t3.4 Writing nonbonded potentials"
    print"\t--------------------------------------"
    from ..md_xml import nonbonded
    current_nonbonded = nonbonded.nonbondedElement(root)
    e14fac = current_nonbonded.elec14_total()
    
    f.write("NONBONDED NBXMOD 5 ATOM CDIEL VATOM VDISTANCE VSWITCH -\n")
    f.write("CUTNB 16.0 CTOFNB 12.0 CTONNB 10.0 E14FAC %8.3f EPS 1.0 WMIN 1.5\n" % e14fac)
    f.write("!"+"~"*70+"\n")
    f.write("!TYPE1          epsilon    rmin/2       epsilon14   rmin14/2 \n")
    f.write("!           [kcal/mol] [Angstroem]    [kcal/mol] [Angstroem]\n")
    f.write("!"+"~"*70+"\n")
    vdws = sorted(parameter.write_vdws(root,energy_conversion,distance_conversion))
    for line in vdws:
        f.write(line)
    print"\t\tNumber of atomic Lennard-Jones= ",len(vdws)

    f.write("\nNBFIX\n")
    f.write("!"+"~"*70+"\n")
    f.write("!TYPE1    TYPE2   epsilon [kcal/mol] rmin [Angstroem]\n")
    f.write("!"+"~"*70+"\n")
    exclusions = sorted(parameter.write_exclusions(root,energy_conversion,distance_conversion))
    for line in exclusions:
        f.write(line)
    print"\t\tNumber of pair Lennard-Jones  = ",len(exclusions)
    f.write("END\n")
    f.close()

