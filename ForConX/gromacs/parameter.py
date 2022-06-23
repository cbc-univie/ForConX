import sys
this = sys.modules[__name__]


def unique_interaction_type(root, inter):
    from ..md_xml import bonds
    from ..md_xml import angles
    from ..md_xml import dihedrals
    itype=dict()
    ifunc=dict()
    if inter == 'bonds':
        ifunc['HARM']=1
        ifunc['MORS']=3
        current_bonds = bonds.bondsElement(root)
        for i in current_bonds.list('HARM MORS'):
              (bondclass,pointer) = current_bonds.find_type(i)
        if len(itype) == 1 :
              return ifunc[itype.keys()[0]]
    if inter == 'angles':
        ifunc['HARM']=1
        ifunc['UREY']=5
        current_angles = angles.anglesElement(root)
        for i in current_angles.list('HARM UREY'):
              (angleclass,pointer) = current_angles.find_type(i)
        if len(itype) == 1 :
              return ifunc[itype.keys()[0]]
    if inter == 'dihedrals':
        ifunc['COS']=9 ; # this is the multiple one
        ifunc['RYCK']=3
        current_dihedrals = dihedrals.dihedralsElement(root)
        for i in current_dihedrals.list('COS RYCK'):
              (dihedralclass,pointer) = current_dihedrals.find_type(i)
        if len(itype) == 1 :
              return ifunc[itype.keys()[0]]

    #this is the same for all interactions
    if len(itype) == 0 : # default
        return 1
    if len(itype) > 1 :
        print "More than a single type of interactions in the same topology. This is incompatible with pdb2gmx. Conversion of interactions will be added."
        exit()


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

#==============================================================================================================================================
# gromacs.mass2pse
# Conversion of mass to PSE number - important for gromacs itp files
#==============================================================================================================================================
def mass2pse(mass):
    if mass < 38:
        pse = round(0.481*mass+0.149)
    elif mass < 42:
        pse = round(16.63*mass**2-1315.5*mass+26030.9)
    elif mass < 67:
        pse = round(0.458*mass+0.16)
        if (mass>56.) and (mass<58.89):
#           Nickel            
            pse = 28
    elif mass < 88:
        pse = round(0.39*mass+3.6)
    elif mass < 115:
        pse = round(0.39*mass+4.5)
    elif mass < 139:
        pse = round(0.33*mass+10.7)
        if (mass >127.2) and (mass < 130):
#           Tellur            
            pse = 52
    elif mass<208.981:
        pse = round(0.357*mass+8.3)
    else:
        pse = round(0.124*mass+58.5)
    return int(pse)

#==================================================================================================
# parameter.write_forcefielditp(dirname)
#==================================================================================================
def create_forcefield_directory(dirname):
   import os
   this.dirname=dirname
   try: 
	os.mkdir(dirname)
   except:
	pass 	
#==================================================================================================
# parameter.write_forcefielditp 
#==================================================================================================
def write_forcefielditp(root):
    """

    -> nonbonded.nonbondedElement
    """
    from ..md_xml import nonbonded

    f = open(this.dirname+'/'+"forcefield.itp",'w')
    f.write("#define _FF_GROMACS\n")
    f.write("#define _FF_s\n\n")

    f.write("[ defaults ]\n")
    f.write(";nbfunc comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
    # Check handling of 1-4 interactions
    current_nonbonded = nonbonded.nonbondedElement(root)
    elec14 = current_nonbonded.elec14_total()
    vdw14 = current_nonbonded.vdw14_total()
    f.write("1     2     yes    %10.2f      %10.2f\n\n"%(elec14,vdw14))
    f.write("; Include bonded and nonbonded parameters\n")
    f.write("#include \"./ffnonbonded.itp\"\n")
    f.write("#include \"./ffbonded.itp\"\n")

#==================================================================================================
# parameter.write_ffnonbonded 
#==================================================================================================
def write_ffnonbonded(root,energy_unit,distance_unit):
    """
    
    -> nonbonded.nonbondedElement
    -> nonbonded.atomClass

    -> molecule.molecule Element
    -> molecule.dihedralClass
    -> molecule.bondClass
    -> molecule.angleClass

    -> input_output.energy_unit
    -> input_output.distance_unit
    """
    from ..md_xml import nonbonded


    f = open(this.dirname+'/'+"ffnonbonded.itp",'w')
    f.write("[ atomtypes ] \n; name  bond_type       mass    charge  ptype   sigma   epsilon\n")

    current_nonbonded = nonbonded.nonbondedElement(root)
    current_nonbonded.vdw2atom()
    polarizable = False
    for i in current_nonbonded.list('atom'):
        nb = nonbonded.atomClass(root,i)
        for iatom in root.findall('./molecule/atom'):
            if i == iatom.get('type'):
                mass = float(iatom.get('mass'))
                pse_number = mass2pse(mass)
                charge = float(iatom.get('charge'))
                alpha = iatom.get('alpha')
                if alpha:
                    # drude mass = 0.3
                   mass = mass - 0.3
                   polarizable = True
        f.write("%5s %5s %8.3f %8.3f  A  %8.3f %8.3f\n"%(i,pse_number,mass,charge,nb.sigma*distance_unit,nb.epsilon*energy_unit))
    if polarizable:
       f.write("%5s %5s %8.3f %8.3f  S  %8.3f %8.3f\n "%("D",1,0.3,-2,0.000,0.0000))

#===== write pair types ===========================================================================
    from ..md_xml import molecule
    from ..md_xml import dihedrals
    from ..md_xml import angles
    import math

    f.write('\n[ pairtypes ]\n')

    pairtypes = []
    for imol in root.findall("molecule"):
        molname = imol.get('name')
        current_molecule = molecule.moleculeElement(imol,molname)
        for i in current_molecule.list('dihedral'):
            dihedral = molecule.dihedralClass(current_molecule.mol,i)
            name     = (dihedral.name).split()
            type     = dihedral.type
            tmp      = type.split()
            write    = True
        # Check for 1-2 or 1-3 interaction
            for i in current_molecule.list('bond'):
                bond   = molecule.bondClass(current_molecule.mol,i)
                name_b = (bond.name).split()
            if ((name_b[0] == name[0]) and (name_b[1]==name[3])):
                write = False
            for i in current_molecule.list('angle'):
                angle  = molecule.angleClass(current_molecule.mol,i)
                name_a = (angle.name).split()
                if ((name_a[0] == name[0]) and (name_a[2]==name[3])):
                   write = False
            if (write):
               elec14 = current_nonbonded.elec14([tmp[0],tmp[3]]) 
               vdw14  = current_nonbonded.vdw14([tmp[0],tmp[3]])  
               nb1 = nonbonded.atomClass(root,tmp[0])
               nb2 = nonbonded.atomClass(root,tmp[3])
               line = ("%4s %4s   1  %8.3f %8.3f\n"%(tmp[0],tmp[3],0.5*(nb1.sigma+nb2.sigma)*distance_unit,(vdw14*math.sqrt(nb1.epsilon*nb2.epsilon)*energy_unit)))
               pairtypes.append(line)
    pairtypes = sorted(set(pairtypes))
    for line in pairtypes:
        f.write(line)

#==================================================================================================
# parameter.write_atomtypes
#==================================================================================================
def write_atomtypes(root):
    """

    -> molecule.moleculeElement
    -> molecule.atomClass
    """
    from ..md_xml import molecule

    f = open(this.dirname+'/'+'atomtypes.atp','w')
    f.write(';atomtypes\n')
    f.write(';type  mass\n')
    type = []
    polarizable = False
    for imol in root.findall("molecule"):
        molname = imol.get('name')
        current_molecule = molecule.moleculeElement(imol,molname)
        for i in current_molecule.list('atom'):
            current_atom = molecule.atomClass(current_molecule.mol,i)
            if current_atom.alpha is None:
                line = '%6s  %8s\n'%(current_atom.type,current_atom.mass)
                type.append(line)
            else:
                line = '%6s  %8s\n'%(current_atom.type,current_atom.mass-0.3)
                polarizable = True
                type.append(line)
    type = sorted(set(type))
    for line in type:
        f.write(line)
    if polarizable:
        f.write('%6s  %8s\n'%('D',0.3))

#==================================================================================================
# parameter.write_ffbonded
#==================================================================================================
def write_ffbonded(root,energy_unit,distance_unit):
    """
    -> bonds.bondsElement
    -> bonds.harmClass
    -> bonds.morsClass

    -> angles.anglesElement
    -> angles.harmClass
    -> angles.ureyClass

    -> dihedrals.dihedralsElement
    -> dihedrals.cosClass
    -> dihedrals.ryckClass

    -> impropers.impropersElement
    -> impropers.cosClass
    -> impropers.harmClass

    -> helpers.warning
    -> helpers.error
    """
    f = open(this.dirname+'/'+"ffbonded.itp",'w')
    import sys
    from ..md_xml import molecule
    from ..md_xml import helpers
    from ..md_xml import input_output

    input_output.conk(f,";")
    input_output.warranty(f,";")
    
#===== write bond types ===========================================================================
    from ..md_xml import bonds

    f.write("[ bondtypes ]\n")
    f.write("; ai    aj  funct  r  k\n")

    current_bonds = bonds.bondsElement(root)
    for i in current_bonds.list('HARM MORS'):
        (bondclass,pointer) = current_bonds.find_type(i)
        if bondclass=="HARM":
           harm = bonds.harmClass(root,i)
           type = i.split()
           r0   = harm.r0 * distance_unit
           k    = harm.k
           if k is None:
              helpers.warning('<bond type="%s"> is constant, but will be set to 750 kcal/mol Angstroem^2.'%i)
              k = 750.0 * energy_unit / distance_unit**2
           else:
              k = 2*float(k) * energy_unit / distance_unit**2
           f.write('%4s %4s   1  %8.4f  %10.7f\n'%(type[0],type[1],r0,k))
        elif bondclass=="MORS":
           mors = bonds.morsClass(root,type)
           f.write('%4s %4s   3  %8.4f  %8.3f  %8.3\n'%(type[0],type[1],float(bond.r0)*distance_unit, float(bond.D0)*energy_unit, float(bond.alpha)/distance_unit))
        else:
           helpers.error('Bond potential %s unknown'%(bondclass))

    polarizable = False
    for iatom in root.findall('./molecule/atom'):
        alpha = iatom.get('alpha')
        if alpha:
           polarizable = True

    if polarizable:
        f.write('%4s %4s   1  %8.4f  %10.7f\n'%("X","D",0.00,418680.00))

#===== write angle types ===========================================================================
    from ..md_xml import angles

    f.write("\n[ angletypes ]\n")
    f.write("; ai    aj   ak    funct  theta0  k\n")

    current_angles = angles.anglesElement(root)
    for i in current_angles.list('HARM UREY'):
        (angleclass,pointer) = current_angles.find_type(i)
        if angleclass=="HARM":
           harm   = angles.harmClass(root,i)
           type   = i.split()
           theta0 = harm.theta0
           k      = 2. * harm.k * energy_unit
           f.write('%4s %4s %4s   1  %s  %8.3f\n'%(type[0],type[1],type[2],theta0,k))
        elif angleclass=="UREY":
           try:
                harm = angles.harmClass(root,type)
           except:
                helpers.error('<urey="%s"> without corresponding harmonic potentia'%type)
           urey = angles.ureyClass(root,type)
           (delta_k, urey_theta0) = urey.harm()
           if (urey_theta0/harm.theta0 >1.02) or (urey_theta0/harm.theta0 <0.98):
               helpers.error('<urey="%s"> Urey r0 and harmonic theta0 do not fit!')
           helpers.warning('Urey-Bradley potential will be applied to corresponding harmonic angle!')
           k = 2*(harm.k + delta_k)*energy_unit
           f.write("   1  %s  %8.3f\n"%(harm.theta0,k))

#===== write dihedral types ===========================================================================
    from ..md_xml import dihedrals

    f.write("\n[ dihedraltypes ]\n")
    f.write("; ai   aj   ak   al   func\n")

    current_dihedrals = dihedrals.dihedralsElement(root)
    for i in current_dihedrals.list("COS RYCK"):
        (dihedralclass,pointer) = current_dihedrals.find_type(i)
        if dihedralclass=="COS":
            cos   = dihedrals.cosClass(root,i)
            type  = i.split()
            n     = cos.n
            k     = cos.k
            delta = cos.delta
        elif dihedralclass=="RYCK":
            ryck        = dihedrals.ryckClass(root,type)
            (n,k,delta) = ryck.cos()
        else:
            helpers.error('The potential from %s is currently not supported by ForConX and Gromacs.'%dihedralclass)
        for j in range(len(k)):
            f.write('%4s %4s %4s %4s   1   %s %8.3f %s\n'%(type[0],type[1],type[2],type[3],delta[j],k[j]*energy_unit,n[j]))

#===== write improper types ===========================================================================
    from ..md_xml import impropers

    f.write(";impropers\n")
    f.write("; ai   aj   ak   al   func\n")

    current_impropers = impropers.impropersElement(root)
    for i in current_impropers.list("COS HARM"):
        (improperclass,pointer) = current_impropers.find_type(i)
        if improperclass=="COS":
            cos   = impropers.cosClass(root,i)
            type  = i.split()
            n     = cos.n
            k     = cos.k
            delta = cos.delta
            for j in range(len(k)):
                f.write('%4s %4s %4s %4s   4   %5.2f %8.3f %4.1d\n'%(type[0],type[1],type[2],type[3],delta[j],k[j]*energy_unit,n[j]))   
        elif improperclass=="HARM":
            harm   = impropers.harmClass(root,i)
            type   = i.split()
            k      = harm.k * 2.0 * energy_unit
            theta0 = harm.theta0
            f.write('%4s %4s %4s %4s   2   %5.3f %8.3f\n'%(type[0],type[1],type[2],type[3],theta0,k))

#==================================================================================================
# parameter.write_rtp
#==================================================================================================
def write_rtp(root,energy_unit,distance_unit):
    """

    -> molecule.moleculeElement
    -> molecule.atomClass
    -> molecule.bondClass
    -> input_output.mdElement
    """
    from ..md_xml import molecule
    from ..md_xml import input_output
    import math

    f = open(this.dirname+'/'+'molecule.rtp','w')
    f.write('[ bondedtypes ]\n')

    polarizability = False
    for i in root.findall('molecule'):
        name = i.get('name')
        current_molecule = molecule.moleculeElement(i,name)
        for i in current_molecule.list('atom'):
            current_atom = molecule.atomClass(current_molecule.mol,i)
            if current_atom.alpha:
                polarizability = True
    # ; Col 1: Type of bond
    # ; Col 2: Type of angles
    # ; Col 3: Type of proper dihedrals
    # ; Col 4: Type of improper dihedrals
    # ; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.
    # ; Col 6: Number of excluded neighbors for nonbonded interactions
    # ; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1
    # ; Col 8: Remove propers over the same bond as an improper if it is 1    
    ifunc_string='   '+str(unique_interaction_type(root,'bonds'))+\
                 '   '+str(unique_interaction_type(root,'angles'))+\
                 '   '+str(unique_interaction_type(root,'dihedrals'))+\
                 '    2        1           3      1     0        '
    if polarizability:
        f.write('; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih Polarizable\n')
        f.write(ifunc_string+'     1\n')
    else:
        f.write('; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih\n')
        f.write(ifunc_string)
    unique_interaction_type(root,'bonds')
    current_output = input_output.mdElement(root,'output')
    polarizability_conversion = current_output.convert_polarizability()
    ctr = 1
    for imol in root.findall("molecule"):
        drudelist = []
        molname = imol.get('name')
        f.write('\n[ %s ]\n'%molname)
        f.write('   [ atoms ]\n')
        f.write(';name  type  charge  cgnr  alpha  thole\n')
        current_molecule = molecule.moleculeElement(imol,molname)
        for i in current_molecule.list('atom'):
            current_atom = molecule.atomClass(current_molecule.mol,i)
            if current_atom.alpha is None:
                f.write('%6s  %6s%6.2f%4d\n'%(current_atom.name,current_atom.type,current_atom.charge,ctr))
                ctr += 1
            else:
                alpha = current_atom.alpha * polarizability_conversion
                q     = math.sqrt(3013.49 * alpha)
                f.write('%6s  %6s%6.2f%4d%10.7f%6s\n'%(current_atom.name,current_atom.type,current_atom.charge+q,ctr,alpha,0.0))
                drudelist.append(current_atom.name)
                ctr += 1
                f.write('%6s  %6s%6.2f%4d\n'%('D'+current_atom.name,'D',-q,ctr))
                ctr += 1

        f.write('\n   [ bonds ]\n')
        for i in current_molecule.list('bond'):
            current_bond = molecule.bondClass(current_molecule.mol,i)
            name         = current_bond.name.split()
            f.write('%6s %6s\n'%(name[0],name[1]))
        for i in drudelist:
            f.write('%6s %6s\n'%(i,'D'+i))

#==================================================================================================
# parameter.write_forcefielddoc
#==================================================================================================
def write_forcefielddoc(root):
    """

    -> molecule.moleculeElement
    -> molecule.atomClass
    """
    from ..md_xml import molecule
    from ..md_xml import input_output

    f = open(this.dirname+'/'+'forcefield.doc','w')
    f.write('\n\nForcefield for')

    polarizability = False
    for i in root.findall('molecule'):
        name = i.get('name')
        f.write(' %s'%name)
        current_molecule = molecule.moleculeElement(i,name)
        for i in current_molecule.list('atom'):
            current_atom = molecule.atomClass(current_molecule.mol,i)
            if current_atom.alpha:
                polarizability = True
    f.write(' converted with ForConX\n\n')

    input_output.header(f,'*.')
    input_output.conk(f,'*.')

    if polarizability:
        f.write('\n\nPolarizable Force Field with Drude Particle, Drude Parameter:\n')
        f.write('mass = 0.3\n')
        f.write('k    = 418680.00\n')
        f.write('q    =\n')

    f.write('\nPlease move all files (atomtypes.atp, ffbonded.itp. ffnonbonded.itp, forcefield.itp, forcefield.doc and ion.rtp) in a new folder Force_Field.ff')
    f.write('Use pdb2gmx to build your simulaiton, the first line in forcefield.doc is the name of the force field in pdb2gmx')

