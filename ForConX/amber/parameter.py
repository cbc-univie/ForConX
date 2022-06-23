#==================================================================================================
# amber.write_bonds
#==================================================================================================
def write_bonds(root,energy_unit,distance_unit):
    from ..md_xml import bonds
    from ..md_xml import helpers

    current_bonds = bonds.bondsElement(root)
    bondlist = []
    
    #   harmonic bonds    
    for i in current_bonds.list('HARM'):
        harm = bonds.harmClass(root,i)
        atomtypes = i.split()
        line =  "%-2s" % atomtypes[0]
        line += "-%-2s    " % atomtypes[1]
#       force constant
        k = harm.k
        if k is None:
            helpers.warning('Harmonic bond "%s %s" is probably constrained, but will be set to 750 kcal/mol Angstroem^2. for compatibility'%(atomtypes[0],atomtypes[1]))
            k = 750.0
        else:
            k = float(k) * energy_unit / distance_unit / distance_unit
        line +="%6s %10.4f %10s" % (" ",k," ")
#       equilibrium distance
        r0 = float(harm.r0) * distance_unit
        line += "%10.4f\n" % r0
        bondlist.append(line)

    #   Morse potentials
    for i in current_bonds.list('MORS'):
        helpers.warning('Morse potentials are not available in AMBER. Therefore, this potential is converted to an harmonic potential.')

        mors = bonds.morsClass(root,i)
        # CHECK if the next line is necessary
        (k,r0) = mors.harm()
        k  = mors.k * energy_unit / distance_unit / distance_unit
        r0 = mors.r0 * distance_unit
        atomtypes = atomtypes.split()
        line =  "%-2s    " % atomtypes[0]
        line += "-%-2s    " % atomtypes[1]
        line += "%6s %10.4f %10s" % (" ",k," ")
        line += "%10.4f\n" % r0
        bondlist.append(line)
    return bondlist


#==================================================================================================
# amber.write_angles
#==================================================================================================
def write_angles(root,energy_unit):
    from ..md_xml import angles
    from ..md_xml import helpers

    current_angles = angles.anglesElement(root)
    anglelist = []

    for i in current_angles.list('HARM'):
        harm = angles.harmClass(root,i)
        atomtypes = i.split()
        line =  "%-2s" % atomtypes[0]
        line += "-%-2s" % atomtypes[1]
        line += "-%-2s" % atomtypes[2]
#       force constant
        k = harm.k * energy_unit
        line += "%2s %10.4f %9s" % (" ",k," ")
#       equilibrium angle
        theta0 = harm.theta0
        line += "%10.4f\n" % theta0
        anglelist.append(line)
    return anglelist

#==================================================================================================
# amber.write_dihedrals
#==================================================================================================
def write_dihedrals(root,energy_unit):
    from ..md_xml import dihedrals
    from ..md_xml import nonbonded
    from ..md_xml import helpers

    dihedrallist = []
    # The current version of ForConX only supports a global scaling of 1-4 interactions. If you need special 1-4 scaling, please adjust the frcmod file.
    print "\t\t\tForConX uses the one scaling factor 1-4. ",
    print "If you want to use special 1-4 interactions, please edit amber_new.frcmod. "

    current_nonbonded = nonbonded.nonbondedElement(root)
    current_dihedrals = dihedrals.dihedralsElement(root)

    for i in current_dihedrals.list("cos"):
        cos = dihedrals.cosClass(root,i)
        atomtypes = i.split()

        elec14 = current_nonbonded.elec14([atomtypes[0],atomtypes[3]])
        vdw14  = current_nonbonded.vdw14([atomtypes[0],atomtypes[3]])
        k = cos.k
        n = cos.n
        delta = cos.delta
        for i in range(len(k)):
            line =  "%-2s" % atomtypes[0]
            line += "-%-2s" % atomtypes[1]
            line += "-%-2s" % atomtypes[2]
            line += "-%-2s" % atomtypes[3]
            line += "  1   %10.4f"   % (float(k[i]) * energy_unit)
            line += "%4s%10.4f" % (" ",delta[i])
            if (len(n)>1):
                if(i==(len(n)-1)):
                    line += "%5s%s"   % (" ",n[i])
                else:
                    line += "%5s%-s"   % (" -",n[i])
            else:
                line += "%5s%s"   % (" ",n[i])
            line += "   SCEE=%1.2f  " % (1.0/elec14)
            line += "   SCNB=%1.2f\n" % (1.0/vdw14)
                
            dihedrallist.append(line)
#   TO DO Ryckaert-Bellemans
    for i in current_dihedrals.list("ryck"):
        ryck        = dihedrals.ryckClass(root,i)
        atomtypes   = i.split()
        elec14 = current_nonbonded.elec14([atomtypes[0],atomtypes[3]])
        vdw14  = current_nonbonded.vdw14([atomtypes[0],atomtypes[3]])
        helpers.warning('The potential form Ryckaert-Bellmanns for dihedral "%s" is not supported by AMBER. It is converted to a cosine potential.' %(i))
        (k,n,delta) = ryck.cos()
        list_k_unequal_zero = []
        for i in range(len(k)):
            if ((float)(k[i]) != 0.0):
                list_k_unequal_zero.append(i)
        for m in range(len(list_k_unequal_zero)):
            i = list_k_unequal_zero[m]
            line =  "%-2s" % atomtypes[0]
            line += "-%-2s" % atomtypes[1]
            line += "-%-2s" % atomtypes[2]
            line += "-%-2s" % atomtypes[3]
            line += "  1   %10.4f"   % (float(k[i]) * energy_unit)
            line += "%4s%10.4f" % (" ",delta[i])
            if (len(list_k_unequal_zero)>1):
                if(m==(len(list_k_unequal_zero)-1)):
                    line += "%5s%s"   % (" ",n[i])
                else:
                    line += "%5s%s"   % (" -",n[i])
            else:
                line += "%5s%s"   % (" ",n[i])
            line += "   SCEE=%1.2f  " % (1.0/elec14)
            line += "   SCNB=%1.2f\n" % (1.0/vdw14)
            dihedrallist.append(line)
        if (len(list_k_unequal_zero)==0):
            line =  "%-2s" % atomtypes[0]
            line += "-%-2s" % atomtypes[1]
            line += "-%-2s" % atomtypes[2]
            line += "-%-2s" % atomtypes[3]
            line += "  1   %10.4f"   % (0.0)
            line += "%4s%10.4f" % (" ",delta[i])
            line += "%5s%s"   % (" ",1)
            line += "   SCEE=%1.2f  " % (1.0/elec14)
            line += "   SCNB=%1.2f\n" % (1.0/vdw14)
            dihedrallist.append(line)
    return dihedrallist

#==================================================================================================
# amber.write_impropers
#==================================================================================================
def write_impropers(root,energy_unit):
    from ..md_xml import impropers
    from ..md_xml import helpers
    import math

    current_impr = impropers.impropersElement(root)
    improperlist = []
    for i in current_impr.list('cos'):
        cos = impropers.cosClass(root,i)
        atomtypes = i.split()
        k = cos.k
        n = cos.n
        delta = cos.delta
        for i in range(len(k)):
            line =  "%-2s" % atomtypes[1]
            line += "-%-2s" % atomtypes[2]
            line += "-%-2s" % atomtypes[0]
            line += "-%-2s" % atomtypes[3]
            line += "  1   %10.4f"   % (float(k[i]) * energy_unit / 2.0)
            line += "%4s%10.4f" % (" ",delta[i])
            line += "%5s%2s\n"   % (" ",n[i])
            improperlist.append(line)
    for i in current_impr.list('harm'):
        harm = impropers.harmClass(root,i)
        atomtypes = i.split()
        k = harm.k
        phi0 = harm.theta0
        line = '%-2s-%2s-%2s-%2s  1  %10.4f   180.000   2\n'%(atomtypes[1],atomtypes[2],atomtypes[0],atomtypes[3],(float(k) * energy_unit /(math.pi*math.pi)))
        helpers.warning('AMBER does not support harmonic impropers for "%s". Therefore, the potential is converted to a cosine function.'%i)
        improperlist.append(line)
    return improperlist



#==================================================================================================
# amber.write_nonbonded
#==================================================================================================
def write_vdw(root,energy_unit,distance_unit):
    from ..md_xml import nonbonded

    current_nonbonded = nonbonded.nonbondedElement(root)
    current_nonbonded.vdw2atom()
    vdws = []
    
    for i in current_nonbonded.list('atom'):
        nb = nonbonded.atomClass(root,i)
        epsilon = -nb.epsilon * energy_unit
        rmin2 = nb.sigma * 0.561231 * distance_unit
        line =  "%-6s" % (i)
        line += "%10.5f" % rmin2
        line += "%10.5f" % epsilon
        line += "\n"
        vdws.append(line)
    return vdws
