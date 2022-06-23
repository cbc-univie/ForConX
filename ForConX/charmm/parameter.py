#======================================================================================
# parameter.read_bonds
#======================================================================================
def read_bonds(root,parameter):
    """
    adds force constants and equilibrium distance to all entries in <bonds>, since
    the topology already created the necessary types.

    -> bonds.bondsElement
    -> bonds.harmClass
    -> bonds.sequence
    """
    from ..md_xml import bonds

    ibonds = 0
    current_bonds = bonds.bondsElement(root)
    for bond_xml in current_bonds.list("HARM"):
        for line in parameter:
            line_element = line.split()
            type = bonds.sequence(' '.join(line_element[0:2]))
            hit = match_type(type.split(),bond_xml.split())
            if hit:
                harm = bonds.harmClass(root,bond_xml)
#               If normal bond exists, wild card potential is disregarded                
                if "X" in line_element[0:2] and harm.get_k() is not None:                    
                    continue
                else:
                    harm.k  = float(line_element[2])
                    harm.r0 = float(line_element[3])
                    ibonds += 1
    return ibonds

#======================================================================================
# parameter.write_bonds
#======================================================================================
def write_bonds(root,energy_unit,distance_unit):
    """

    -> bonds.bondsElement
    -> bonds.harmClass
    -> bonds.morsClass

    -> helpers.warning
    -> helpers.error
    """
    from ..md_xml import bonds
    from ..md_xml import helpers

    current_bonds = bonds.bondsElement(root)
    parameter = []
    for i in current_bonds.list('HARM MORS'):
        (bondclass,pointer) = current_bonds.find_type(i)
        if bondclass=="HARM":
            harm = bonds.harmClass(root,i)
            type = i.split()
            r0   = harm.r0 * distance_unit
            k    = harm.k
            if k is None:
                helpers.warning('<bond type="%s"> is constant, but will be set to 750 kcal/mol Angstroem^2.'%i)
                k = 750.0
            else:
                k = float(k) * energy_unit / distance_unit / distance_unit
        elif bondclass=="MORS":
            helpers.warning('<mors type="%s"> is not available in CHARMM and consequently converted to an harmonic bond.'%i)
            mors = bonds.morsClass(root,i)
            type = i.split()
            (k,r) = mors.harm()
            k  *= energy_unit / distance_unit / distance_unit
            r0 *= distance_unit
        else:
            helpers.error('<bonds> %s potential is unknown' % bondclass)
        
        line =  "%-6s    " % type[0]
        line += "%-6s    " % type[1]
        line +="%6s %10.4f %10s" % (" ",k," ")
        line += "%10.4f\n" % r0
        parameter.append(line)
        
#   Adding constant Drude harmonic potential, i.e. the effective Drude charge depends on the
#   polarizability
    if root.findall('./molecule/atom/[@alpha]'):
        line =  "%-6s    " % ("DRUD")
        line += "%-6s    " % ("X")
        line += "%6s %10.4f %10s" % (" ",500.0," ")
        line += "%10.4f\n" % (0.0)
        parameter.append(line)
    return parameter

#======================================================================================
# parameter.read_angles
#======================================================================================
def read_angles(root,parameter):
    """
    adds force constants and equilibrium angles to all entries in <angles>, since
    the topology already created the necessary types.

    -> angles.anglesElement
    -> angles.harmClass
    -> angles.ureyClass
    -> angles.sequence
    """
    from ..md_xml import angles
    
    iangles = 0
    current_angles = angles.anglesElement(root)
    for angle_xml in current_angles.list("HARM"):
        for line in parameter:
            line_element = line.split()
            type = angles.sequence(' '.join(line_element[0:3]))
            hit = match_type(type.split(), angle_xml.split())
            if hit:
                harm = angles.harmClass(root,angle_xml)
#               If normal angle exists, wild card potential is disregarded                
                if "X" in line_element[0:3] and harm.k is not None:                    
                    continue
                else:
                    harm.k      = float(line_element[3])
                    harm.theta0 = float(line_element[4])

#                   Urey-Bradley
                    if len(line_element)>5:
                        urey = angles.ureyClass(root,angle_xml)
                        urey.k  = float(line_element[5])
                        urey.r0 = float(line_element[6])
                iangles += 1
    return iangles

#======================================================================================
# parameter.write_angles
#======================================================================================
def write_angles(root,energy_unit,distance_unit):
    """

    -> angles.anglesElement()
    -> angles.harmClass()
    -> angles.UreyClass()
    """
    from ..md_xml import angles

    current_angles = angles.anglesElement(root)
    tmp = {}
    for i in current_angles.list('HARM'):
        harm = angles.harmClass(root,i)
        type = i.split()
        line =  "%-6s    " % type[0]
        line += "%-6s    " % type[1]
        line += "%-6s    " % type[2]
        k = harm.k * energy_unit 
        line += "%2s %10.4f %9s" % (" ",k," ")
        theta0 = harm.theta0
        line += "%10.4f" % theta0
        tmp[i] = line
    for i in current_angles.list('UREY'):
        urey = angles.ureyClass(root,i)
        if i in tmp:
            line  = "%10.4f" % (float(urey.k)  * energy_unit)  
            line += "%10.4f" % (float(urey.r0) * distance_unit)
        else:
            type = i.split()
            line =  "%-6s    " % type[0]
            line += "%-6s    " % type[1]
            line += "%-6s    " % type[2]
            line += "%2s %10.4f %9s" % (" ",0.0," ")
            line += "%10.4f" % 0.0
            line  = "%10.4f" % (float(urey.k)  * energy_unit)  
            line += "%10.4f" % (float(urey.r0) * distance_unit)
        tmp[i] += line 

    parameter = []
    for i in tmp:
        parameter.append(tmp[i]+'\n')
    return parameter

#======================================================================================
# parameter.read_dihedrals
#======================================================================================
def read_dihedrals(root,parameter):
    """

    -> dihedrals.dihedralsElement
    -> dihedrals.cosClass
    -> dihedrals.sequence
    """
    from ..md_xml import dihedrals

    idihedrals = 0
    current_dihedrals = dihedrals.dihedralsElement(root)
    for dihedral_xml in current_dihedrals.list("COS"):
        for line in parameter:
            line_element = line.split()
            type = dihedrals.sequence(' '.join(line_element[0:4]))

#           Wildcards may change the alphabetic order of the middle atoms
            hit1 = match_type(type.split(), dihedral_xml.split())
            hit2 = match_type((type.split())[::-1], dihedral_xml.split())
            if hit1 or hit2:
                cos = dihedrals.cosClass(root,dihedral_xml)
                old_n = cos.n
                old_k = cos.k
                old_delta = cos.delta

                if "X" in line_element[0:4] and len(old_n)>0:
                    continue
                
                current_k = float(line_element[4])
                current_n = int(line_element[5])
                current_delta = float(line_element[6])

                if current_n in old_n:
                    continue
                else:
                    cos.n = old_n + [current_n]
                    cos.k = old_k + [current_k]
                    cos.delta = old_delta + [current_delta]
                idihedrals += 1
    return idihedrals

#======================================================================================
# parameter.write_dihedrals
#======================================================================================
def write_dihedrals(root,energy_unit):
    """
    
    -> dihedrals.dihedralsElement
    -> dihedrals.cosClass
    -> dihedrals.ryckClass

    """
    from ..md_xml import dihedrals

    current_dihedrals = dihedrals.dihedralsElement(root)
    parameter = []
    for i in current_dihedrals.list('COS RYCK'):
        (dihedralclass,pointer) = current_dihedrals.find_type(i)
        if dihedralclass=="COS":
            cos  = dihedrals.cosClass(root,i)
            type = i.split()
            n = cos.n
            k = cos.k
            delta =cos.delta
        elif dihedralclass=="RYCK":
            ryck  = dihedrals.ryckClass(root,i)
            type  = i.split()
            (n,k,delta) = ryck.cos()
        else:
            helpers.warning('<dihedrals> %s potential is unknown'%dihedralclass)

        first = True
        for i in range(len(n)):
            if not first and k[i]==0:
                continue
            line =  "%-6s    " % type[0]
            line += "%-6s    " % type[1]
            line += "%-6s    " % type[2]
            line += "%-6s    " % type[3]
            line += "%10.4f"   % (float(k[i]) * energy_unit)
            line += "%4s%2s"   % (" ",int(n[i]))
            line += "%4s%10.4f\n" % (" ",delta[i])
            parameter.append(line)
            first = False
    return parameter

#======================================================================================
# parameter.read_improper
#======================================================================================
def read_impropers(root,parameter):
    """

    -> impropers.impropersElement
    -> impropers.cosClass
    -> impropers.harmClass
    -> impropers.sequence 
    """
    from ..md_xml import impropers
    
    iimpropers = 0
    current_impropers = impropers.impropersElement(root)
    for improper_xml in current_impropers.list("COS"):
        for line in parameter:
            line_element = line.split()
            type = impropers.sequence(' '.join(line_element[0:4]))
            hit1 = match_type(type.split(), improper_xml.split())
            hit2 = match_type(type.split()[::-1], improper_xml.split())
            if hit1 or hit2:
                k = float(line_element[4])
                n = float(line_element[5])
                delta = float(line_element[6])
                cos = impropers.cosClass(root,improper_xml)
                if n==0.0:
                    cos.remove()
                    harm = impropers.harmClass(root,improper_xml)
                    if "X" in line_element[0:4] and harm.k is not None:
                        continue
                    harm.k = k
                    harm.theta0 = delta
                else:
                    old_n = cos.n
                    old_k = cos.k
                    old_delta = cos.delta
                    if "X" in line_element[0:4] and old_n is not None:
                        continue
                    cos.add(n,k,delta)
                iimpropers += 1
    return iimpropers

#======================================================================================
# parameter.write_impropers
#======================================================================================
def write_impropers(root,energy_unit):
    """

    -> impropers.impropersElement
    -> impropers.cosClass
    -> impropers.harmClass

    -> helpers.warning
    """
    from ..md_xml import molecule
    from ..md_xml import impropers
    from ..md_xml import helpers
    
    current_impropers = impropers.impropersElement(root)
    parameter = []
    
    for i in current_impropers.list('COS HARM RYCK'):
        (improperclass,pointer) = current_impropers.find_type(i)
        if improperclass=="COS":
            cos = impropers.cosClass(root,i)
            type = i.split()
            n = cos.n
            k = cos.k
            delta = cos.delta
            for i in range(len(n)):
                if float(k[i]) == 0.:
                    continue
                line =  "%-6s    " % type[0]
                line += "%-6s    " % type[1]
                line += "%-6s    " % type[2]
                line += "%-6s    " % type[3]
                line += "%10.4f"   % (float(k[i]) * energy_unit)
                line += "%4s%2s"   % (" ",int(n[i]))
                line += "%4s%10.4f\n" % (" ",delta[i])
                parameter.append(line)
        elif improperclass=="RYCK":
            ryck = impropers.ryckClass(root,i)
            helpers.warning("<impropers/ryck>  %s will be converted to <impropers/cos>."%i)
            (n,k,delta) = ryck.cos()
            type = i.split()
            for i in range(len(n)):
                if float(k[i]) == 0.:
                    continue
                line =  "%-6s    " % type[0]
                line += "%-6s    " % type[1]
                line += "%-6s    " % type[2]
                line += "%-6s    " % type[3]
                line += "%10.4f"   % (float(k[i]) * energy_unit)
                line += "%4s%2s"   % (" ",int(n[i]))
                line += "%4s%10.4f\n" % (" ",delta[i])
                parameter.append(line)
        elif improperclass=="HARM":
            harm   = impropers.harmClass(root,i)
            type   = i.split()
            k      = harm.k
            theta0 = harm.theta0
            line =  "%-6s    " % type[0]
            line += "%-6s    " % type[1]
            line += "%-6s    " % type[2]
            line += "%-6s    " % type[3]
            line += "%10.4f"   % (float(k) * energy_unit)
            line += "%5s%2s"   % (" ",0.)
            line += "%4s%10.4f\n" % (" ",theta0)
            parameter.append(line)
        else:
            helpers.warning('<impropers> %s potential is unknown'%improperclass)
    return parameter

#======================================================================================
# parameter.read_nonbonded
#======================================================================================
def read_nonbonded(root,parameter_atom,parameter_vdw,elec14):
    """

    -> nonbonded.nonbondedElement
    -> nonbonded.atomClass
    -> nonbonded.vdwClass

    -> helpers.warning
    """
    from ..md_xml import nonbonded
    from ..md_xml import molecule
    from ..md_xml import helpers
    
    iatom = 0
    current_nonbonded = nonbonded.nonbondedElement(root)
    current_nonbonded.mixing_epsilon="geometric"
    current_nonbonded.mixing_sigma="arithmetic"
    
#   Generating typelist    
    typelist = []
    for i in root.findall("molecule/atom"):
        typelist.append(i.get("type"))
    typelist = sorted(set(typelist))

#   Looking for Lennard-Jones parameters of these types    
    for line in parameter_atom:
        line_element = line.split()
        type = line_element[0]
        if not type in typelist:
            continue
        current_atom = nonbonded.atomClass(root,type)
        current_atom.epsilon = -float(line_element[2])
#       sigma = 2^(5/6) rmim/2 
        current_atom.sigma   = float(line_element[3])*1.781797
        current_atom.elec14  = float(elec14)

        if len(line_element)<5:
            current_atom.vdw14 = 1.0
        else:
            if line_element[3] != line_element[6]:
                helpers.warning('<atom type="%s"> sigma14 differs from sigma! Extrapolating vdw14 ...'%type)
                C12_14 = 4. * float(line_element[5])* float(line_element[6])**12
                C6_14  = 4. * float(line_element[5])* float(line_element[6])**6

                C12    = 4. * float(line_element[2])* float(line_element[3])**12
                C6     = 4. * float(line_element[2])* float(line_element[3])**6

                current_atom.vdw14 = 0.5 * (C12_14/C12 + C6_14/C6)
            else:
                try:
                    current_atom.vdw14 = float(line_element[5])/float(line_element[2])
                except:
                    current_atom.vdw14 = 1.0
        iatom +=1

        #   Pairwise interactions
    ivdw = 0
    for line in parameter_vdw:
        line_element = line.split()
        if line_element[0] not in typelist:
            continue
        if line_element[1] not in typelist:
            continue    
        current_vdw = nonbonded.vdwClass(root,' '.join([line_element[0],line_element[1]]))
        current_vdw.set_epsilon(-float(line_element[2]))
        current_vdw.set_sigma(float(line_element[3])*1.781797/2.)
        if len(line_element)<5:
            current_atom.vdw14 = 1.0
        else:
            try:
                current_vdw.vdw14 = float(line_element[4])/float(line_element[2])
            except:
                current_vdw.vdw14 = 1.0
        ivdw += 1
    return iatom, ivdw

#======================================================================================
# parameter.write_vdws
#======================================================================================
def write_vdws(root,energy_unit,distance_unit):
    """

    -> nonbonded.nonbondedElement
    -> nonbonded.atomClass

    -> input_output.energy_unit
    -> input_output.distance_unit
    """
    from ..md_xml import nonbonded
    current_nonbonded = nonbonded.nonbondedElement(root)
    current_nonbonded.vdw2atom()
    parameter = []
    for i in current_nonbonded.list('atom'):
        nb = nonbonded.atomClass(root,i)
        epsilon = -nb.epsilon * energy_unit
#       rmin/2 = 2^(-5/6) sigma
        rmin2 = nb.sigma * 0.561231 * distance_unit
        line =  "%-6s " % i
        line += "%5.2f %10.5f" % (0., epsilon)
        line += "%10.5f" % rmin2
#       1-4 interactions
        vdw14 = nb.vdw14
        if vdw14 is not None and vdw14 != 1.0:
            epsilon14 = vdw14*epsilon
            line += "%5.2f %10.5f" % (0.,epsilon14)
            line += "%10.5f" % rmin2
        line += "\n"
        parameter.append(line)
#   Polarizable Drude particles        
    if root.findall('./molecule/atom/[@alpha]'):
        line =  "%-6s %5.2f " % ("DRUD",0.)
        line += "%10.5f%10.5f" % (0.0,0.0)
        line += "%5.2f %10.5f%10.5f\n" % (0.0,0.0,0.0)
        parameter.append(line)
    return parameter

#======================================================================================
# parameter.write_exclusions
#======================================================================================
def write_exclusions(root,energy_unit,distance_unit):
    """

    -> nonbonded.nonbondedElement
    -> nonbonded.vdwClass
    """
    from ..md_xml import nonbonded

    current_nonbonded = nonbonded.nonbondedElement(root)
    parameter = []
    for i in current_nonbonded.list('vdw'):
        nb = nonbonded.vdwClass(root,i)
        epsilon = -nb.epsilon * energy_unit
#       rmin = 2^(1/6) sigma 
        rmin = nb.sigma * 1.122462 * distance_unit
        line =  "%-6s " % i
        line += "%10.5f %10.5f" % (epsilon,rmin)
#       1-4 interactions
        vdw14 = nb.vdw14
        if vdw14 is not None and vdw14 != 1.0:
            epsilon14 = vdw14*epsilon
            line += "%10.5f %10.5f" % (epsilon14,rmin)
        line += "\n"
        parameter.append(line)
    return parameter

#======================================================================================
# parameter.match_type
#======================================================================================
def match_type(type_ff,type_xml):
    """
    In CHARMM X is used as wildcard. 
    This routine checks if current xml type of a bond, angle, dihedral or improper
    coincides with one of the force field.
    """
    hit = type_ff.count('X') 
    for i in range(len(type_ff)):
        if type_ff[i] == type_xml[i]:
            hit +=1
    return hit==len(type_ff)

