import xml.etree.ElementTree as ET

#===============================================================================================
# xml2md.write_field
#===============================================================================================
def write_field(root):
    """
    """
    from ..md_xml import molecule
    from ..md_xml import input_output
    from ..md_xml import helpers

    # header
    pointer = root.find('./output/field')
    try:
        field_file = pointer.get('file')
        print"> Writing field file %s ..."%field_file
        f = open(field_file,'w')
    except IOError:
        helpers.error('Field file cannot be opened.')
    except:
        helpers.error("Can't find field information in <output>")

    current_output = input_output.mdElement(root,'output')
    energy_conversion   = current_output.convert_energy()
    distance_conversion = current_output.convert_distance()
    polarizability_conversion = current_output.convert_polarizability()

    input_output.header(f,"#")
    input_output.conk(f,"#")
    f.write("#"+"."*70+"\n")
    input_output.warranty(f,"#")
    f.write('units %10s \n'%current_output.energy_unit)

    molecules = root.findall('./molecule')
    f.write('molecules %i \n'%len(molecules))
    
    # topology    
    print"\n\t----------------------------------------------------------"
    print"\t3.1 Writing molecules including intramolecular interactions"
    print"\t----------------------------------------------------------"
    print"\n\tVirtual atoms are currently not supported for DLPOLY.\n"
    for mol in molecules:
        molname = mol.get('name')
        print'\t<molecule name="%s">' % molname
        current_molecule = molecule.moleculeElement(mol,molname)
        f.write('%-8s\n'%molname)
        f.write('nummols %10s\n'%mol.get('nmol'))

        # Atoms
        (atomlist,atoms) = write_atom(current_molecule)
        f.write("atoms   %10s\n"%len(atoms))
        for line in atoms:
            f.write(line)
        print "\t\tAtoms ...       %6s" % len(atoms)

        # Bonds
        (bonds,constraints,shells) = write_bonds(root,
                                                 current_molecule,
                                                 atomlist,
                                                 distance_conversion,
                                                 polarizability_conversion,
                                                 energy_conversion)
        nbond = len(bonds)
        print "\t\tBonds ...       %6s" % nbond
        if nbond>0:
            f.write("bonds %-10s\n" % nbond)
            for line in sorted(bonds):
                f.write(line)

        # Contraints        
        nconstraint = len(constraints)
        if nconstraint>0:
            print "\t\tConstraints ... %6s" % nconstraint
            f.write("constraints %-10s\n" % nconstraint)
            for line in sorted(constraints):
                f.write(line)

        # Shells        
        nshell = len(shells)
        if nshell>0:
            print '\t\tDrudes ...      %6s'%nshell
            f.write("shell %6s    1\n" % nshell)
            for line in shells:
                f.write(line)

        # Angles        
        angles = write_angles(root,
                              current_molecule,
                              atomlist,
                              energy_conversion)
        nangle = len(angles)
        print "\t\tAngles ...      %6s"%nangle
        if nangle>0:
            f.write("angles %-10s\n" % len(angles))
            for line in sorted(angles):
                f.write(line)

        # Dihedrals and impropers        
        (dihedrals,impropers) = write_dihedrals(root,
                                                current_molecule,
                                                atomlist,
                                                energy_conversion)
        ndihedral = len(dihedrals)
        print "\t\tDihedrals ...   %6s" %ndihedral
        nimproper = len(impropers)
        print "\t\tImpropers ...   %6s" % nimproper
        if ndihedral+nimproper>0:
            f.write("dihedrals %-10s\n" % (ndihedral+nimproper))
            for line in sorted(dihedrals):
                f.write(line)
            for line in sorted(impropers):
                f.write(line)
        f.write("finish\n")
        print
    
    print"\n\t--------------------------------------"
    print"\t3.2 Writing nonbonded potentials"
    print"\t--------------------------------------"
    vdw = write_vdw(root,distance_conversion,energy_conversion)
    nvdw = len(vdw)
    print "\t\tVan-der-Waals ... ",nvdw
    f.write("vdw %-10s\n" % (nvdw))
    for line in sorted(vdw):
        f.write(line)
    f.write('close')
    f.close()

#===============================================================================================
# dlpoly.write_atom
#===============================================================================================
def write_atom(current_molecule):
    """
    """
    from ..md_xml import molecule

    atomlist  = []
    parameter = []
    for i in  current_molecule.list('atom'):
        current_atom = molecule.atomClass(current_molecule.mol,i)
        if current_atom.alpha is None:
            line  = "%-8s "  % current_atom.type
            line += "%8.3f " % current_atom.mass
            line += "%8.3f " % current_atom.charge
            line += "1 \n"
            parameter.append(line)
            atomlist.append(current_atom.name)
        else:
            # Normal polarizable atom
            line  = "%-8s "  % current_atom.type
            line += "%8.3f " % (current_atom.mass-0.3)
            line += "%8.3f " % (current_atom.charge+2.0)
            line += "1 \n"
            parameter.append(line)
            atomlist.append(current_atom.name)

            # Followed by Drude particle
            drude = "D"+current_atom.type
            line  = "%-8s "  % drude
            line += "%8.3f " % 0.3
            line += "%8.3f " % -2.0
            line += "1 \n"
            parameter.append(line)
            atomlist.append(drude)
    return atomlist,parameter
        
#===============================================================================================
# dlpoly.write_bonds
#===============================================================================================
def write_bonds(root, current_molecule, atomlist,
                distance_conversion, polarizability_conversion, energy_conversion):
    """
    """
    from ..md_xml import molecule
    from ..md_xml import bonds

    current_bonds = bonds.bondsElement(root)
    bond_lines = []
    constraint_lines = []

    # Bond and constraint potentials    
    for i in current_molecule.list('bond'):
        bond = molecule.bondClass(current_molecule.mol,i)
        type = bond.type
        name = (bond.name).split()        
        (potential,pointer) = current_bonds.find_type(type)
        if potential=="HARM":
            harm = bonds.harmClass(root,type)
            r0 = harm.r0 * distance_conversion
            k  = harm.k   
            if k is None:
                line  = "     "
                line += "%8s " % str(atomlist.index(name[0])+1)
                line += "%8s " % str(atomlist.index(name[1])+1)
                line += "%10s" % " "
                line += "%8.3f\n" % r0
                constraint_lines.append(line)
            else:
                k = float(k)*energy_conversion/distance_conversion**2
                line = "harm "
                line += "%8s " % str(atomlist.index(name[0])+1)
                line += "%8s " % str(atomlist.index(name[1])+1)
                line += "%10.5f" % (k * 2.)
                line += "%8.3f\n" % r0
                bond_lines.append(line)
        if potential=="MORS":
            mors = bonds.morsClass(root,type)
            D0    = mors.D0   * energy_conversion
            alpha = mors.beta / distance_conversion
            r0    = mors.r0   * distance_conversion
            line = "mors "
            line += "%8s " % str(atomlist.index(name[0])+1)
            line += "%8s " % str(atomlist.index(name[1])+1)
            line += "%10.5f" % D0
            line += "%8.3f"  % r0
            line += "%10.5f\n" % alpha
            bond_lines.append(line)

    # Shell        
    shell_lines = []
    for i in range(1,len(atomlist)):
        if (atomlist[i])[0] =="D":
            atom_i = molecule.atomClass(current_molecule.mol,atomlist[i-1])
            alpha = atom_i.alpha * polarizability_conversion
            k = 1602.2 / (3.1415926*8.85) / alpha
            # additional factors!!!
            line  = "     "
            line += "%8s " % str(i)
            line += "%8s " % str(i+1)
            line += "%10.5f\n" % k
            shell_lines.append(line)
            # pseudo-Thole: induced dipolar interactions are excluded between atoms connected by
            # a bond or an angle
            for j in range(i+1,len(atomlist)):
                if (atomlist[j])[0] =="D":
                    path = current_molecule.find_path(atomlist[i-1],atomlist[j-1])
                    if len(path)<4:
                        line = "harm "
                        line += "%8s " % str(i+1)
                        line += "%8s " % str(j+1)
                        line += "%10.5f" % 0.0
                        line += "%8.3f\n" % 0.0
                        bond_lines.append(line)
    return bond_lines,constraint_lines,shell_lines
            
#===============================================================================================
# dlpoly.write_angles
#===============================================================================================
def write_angles(root,current_molecule,atomlist,energy_conversion):
    """
    """
    from ..md_xml import molecule
    from ..md_xml import angles
    from ..md_xml import bonds
    from ..md_xml import helpers
    import numpy
    current_angles = angles.anglesElement(root)
    tmp_list = {}
    
    for i in current_molecule.list('angle'):
        angle = molecule.angleClass(current_molecule.mol,i)
        type = angle.type
        name = (angle.name).split()
        key  = "%8s" % str(atomlist.index(name[0])+1)
        key += "%8s" % str(atomlist.index(name[1])+1)
        key += "%8s" % str(atomlist.index(name[2])+1)

        (potential,pointer) = current_angles.find_type(type)
        if potential=="HARM":
            harm   = angles.harmClass(root,type)
            theta0 = harm.theta0 
            k      = 2.*harm.k * energy_conversion
            tmp_list[key] = "%10.5f%8.3f\n"%(k,theta0)
        if potential=="UREY":
            # looking for corresponding harmonic potential
            try:
                harm = angles.harmClass(root,type)
            except:
                helpers.error('<urey="%s"> without corresponding harmonic potential'%type)
                
            # compare r0 and theta0    
            urey   = angles.ureyClass(root,type)
            (delta_k, urey_theta0) = urey.harm()
            if (urey_theta0/harm.theta0 >1.02) or (urey_theta0/harm.theta0 <0.98):
                helpers.error('<urey="%s"> Urey r0 and harmonic theta0 do not fit!')
            
            helpers.warning('Urey-Bradley potential will be applied to corresponding harmonic angle!')
            k = harm.k + delta_k
            tmp_list[key] = "%10.5f%8.3f\n"%(k,harm.theta0)

    angle_lines = []
    for i in tmp_list:
        line  = 'harm '+i 
        line += tmp_list[i]
        angle_lines.append(line)
    return angle_lines
            
#===============================================================================================
# dlpoly.write_dihedrals
#===============================================================================================
def write_dihedrals(root,current_molecule,atomlist,energy_conversion):
    """
    """
    from ..md_xml import molecule
    from ..md_xml import dihedrals
    from ..md_xml import impropers
    from ..md_xml import nonbonded
    from ..md_xml import helpers
    current_nonbonded = nonbonded.nonbondedElement(root)
    
#   Dihedrals    
    current_dihedrals = dihedrals.dihedralsElement(root)    
    dihedral_lines = []
    
    for i in current_molecule.list('dihedral'):
        dihedral = molecule.dihedralClass(current_molecule.mol,i)
        type     = dihedral.type
        tmp      = type.split()
        current_vdw = nonbonded.vdwClass(root,' '.join([tmp[0],tmp[3]]))
        elec14      = current_vdw.elec14
        vdw14       = current_vdw.vdw14
        (potential,pointer) = current_dihedrals.find_type(type)
        name     = (dihedral.name).split()
        if potential=="COS":
            cos = dihedrals.cosClass(root,type)
            n = cos.n
            k = cos.k
            delta =cos.delta
        if potential=="RYCK":
            ryck  = dihedrals.ryckClass(root,type)
            helpers.warning('<impropers/ryck> converted to cos.')
            (n,k,delta) = ryck.cos()
        for j in range(len(n)):
            line = "cos  "
            line += "%8s"   % str(atomlist.index(name[0])+1)
            line += "%8s"   % str(atomlist.index(name[1])+1)
            line += "%8s"   % str(atomlist.index(name[2])+1)
            line += "%8s"   % str(atomlist.index(name[3])+1)
            line += "%10.5f" % (k[j] * energy_conversion) 
            line += "%10.5f" % delta[j]
            line += "%8.3f"  % n[j]
            line += "%8.3f"  % elec14
            line += "%8.3f\n"% vdw14
            dihedral_lines.append(line)
                
#   Impropers
    current_impropers = impropers.impropersElement(root)
    improper_lines = []
    for i in current_molecule.list('improper'):
        improper = molecule.improperClass(current_molecule.mol,i)
        type     = improper.type
        (potential,pointer) = current_impropers.find_type(type)
        name = (improper.name).split()
        if potential=="COS":
            cos = impropers.cosClass(root,type)
            for j in range(len(cos.k)):
                line = "cos  "
                line += "%8s"   % str(atomlist.index(name[0])+1)
                line += "%8s"   % str(atomlist.index(name[1])+1)
                line += "%8s"   % str(atomlist.index(name[2])+1)
                line += "%8s"   % str(atomlist.index(name[3])+1)
                line += "%10.5f" % (cos.k[j] * energy_conversion)
                line += "%10.5f" % cos.delta[j]
                line += "%8.3f"  % cos.n[j]
                line += "%8.3f"  % 0.0
                line += "%8.3f\n"% 0.0
                improper_lines.append(line)
        elif potential=="RYCK":
            ryck  = impropers.ryckClass(root,type)
            helpers.warning('<impropers/ryck> converted to cos.')
            (n,k,delta) = ryck.cos()
            for j in range(len(n)):
                line = "cos  "
                line += "%8s"   % str(atomlist.index(name[0])+1)
                line += "%8s"   % str(atomlist.index(name[1])+1)
                line += "%8s"   % str(atomlist.index(name[2])+1)
                line += "%8s"   % str(atomlist.index(name[3])+1)
                line += "%10.5f" % (k[j] * energy_conversion) 
                line += "%10.5f" % delta[j]
                line += "%8.3f"  % n[j]
                line += "%8.3f"  % elec14
                line += "%8.3f\n"% vdw14
                improper_lines.append(line)
        elif potential=="HARM":
            harm = impropers.harmClass(root,type)
            line = "harm "
            line += "%8s"    % str(atomlist.index(name[0])+1)
            line += "%8s"    % str(atomlist.index(name[1])+1)
            line += "%8s"    % str(atomlist.index(name[2])+1)
            line += "%8s"    % str(atomlist.index(name[3])+1)
            line += "%10.5f"  % (2.0*harm.k* energy_conversion)
            line += "%10.5f"  % 0.0
            line += "%8.3f"   % 0.0
            line += "%8.3f"   % 0.0
            line += "%8.3f\n" % 0.0
            improper_lines.append(line)
    return dihedral_lines,improper_lines
            
#===============================================================================================
# dlpoly.write_vdw
#===============================================================================================
def write_vdw(root,distance_conversion,energy_conversion):
    from ..md_xml import nonbonded
    
    current_nonbonded = nonbonded.nonbondedElement(root)
    vdw_lines = []
    for i in current_nonbonded.list('vdw'):
        current_vdw = nonbonded.vdwClass(root,i)
        type = i.split()
        lj  = "%8s " % type[0]
        lj += "%8s " % type[1]
        lj += "%8s " % "lj"
        lj += "%10.7f "   % (current_vdw.epsilon * energy_conversion)
        lj += "%10.7f \n" % (current_vdw.sigma   * distance_conversion)
        vdw_lines.append(lj)
    return vdw_lines

        
