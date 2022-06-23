import xml.etree.ElementTree as ET
import sys

#==================================================================================================
# amber_write_lib
# Lib files contain the topology of each molecule and partial charges
# Subroutine:
#    - mass2pse
#==================================================================================================
def write_lib(root):

    print"\n\t----------------------------------"
    print"\t3.1 Writing lib files"
    print"\t----------------------------------"
    print"\t\tFor more information on the file format see http://ambermd.org/formats.html"
    print"\t\tIf you want to connect two residues you have to modify section residueconnect in the lib file.\n"

    for molecule in root.findall('molecule'):
        molname = molecule.get('name')[:3].rstrip()
        filename = ''.join([molname,".lib"])
        print"\t\tMolecule ",molname," -> ",filename
        f = open(filename,'w')
        f.write('!!index array str\n')
        f.write(' "%s"\n'%molname)
        f.write('!entry.%s.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n'%molname)
        atomlist = []
        index = 1
        # Via the bond section the program checks which atoms are present in the molecule
        for i in molecule.getchildren():
            if i.tag == 'bond':
                atomnames = i.get("name").split()
                if atomnames[0] not in atomlist:
                    atomlist.append(atomnames[0])
                    sorted(atomlist)
                if atomnames[1] not in atomlist:
                    atomlist.append(atomnames[1])
                    sorted(atomlist)
        indices = {}
        print "\t\t\t",
        count = 0
        for i in molecule.getchildren():
            if i.tag == 'atom':
                count +=1
        if count ==1:
            for i in molecule.getchildren():
                atomlist.append(i.get('name'))
        for child in molecule.getchildren():
            if child.tag == 'atom':
                if child.get('name') in atomlist:
                    pse_number = mass2pse(float(child.get('mass')))
                    print "%5s, " % child.get('name'),
                    if (index%5==0):
                        print "\n\t\t\t",
                    indices[child.get('name')] = index
                    line  = ' "%s" "%s" ' % (child.get('name'), child.get('type'))
                    line += '0 1  17956867 %i %i ' % (index, pse_number)
                    line += '%s \n' % child.get('charge')
                    f.write(line)
                    index +=1
        print
        f.write('!entry.%s.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n'%molname)
        for child in molecule.getchildren():
            if child.tag == 'atom':
                    if child.get('name') in atomlist:
                        f.write('"%s" "%s" 0 -1 0.0\n' % (child.get('name'), child.get('type')))
        f.write('!entry.%s.unit.boundbox array dbl\n -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n'%molname)
        f.write('!entry.%s.unit.childsequence single int\n 2\n'%molname)
        f.write('!entry.%s.unit.connect array int\n 0\n 0\n'%molname)
        f.write('!entry.%s.unit.connectivity table  int atom1x  int atom2x  int flags\n'%molname)
        start = 0
        index = 0
        for i in molecule.getchildren():
            if i.tag == 'bond':
                # start variable is necessary because the bonds assignment in the lib file is via integers
                atomnames = i.get('name').split()
                line  = ' %i' % indices[atomnames[0]]
                line += ' %i' % indices[atomnames[1]]
                line += ' 1\n'
                f.write(line)
                index+=1
        f.write('!entry.%s.unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx\n'%molname)
        index = 1
        f.write(' "U" 0 "R" %i\n'%index)
        for i in range(len(atomlist)):
            f.write(' "R" 1 "A" %i\n'%index)
            index+=1
        f.write('!entry.%s.unit.name single str\n "%s"\n' %(molname,molname))
        f.write('!entry.%s.unit.positions table  dbl x  dbl y  dbl z\n'%molname)
        for i in range(len(atomlist)):
            f.write(' 0.0 0.0 0.0\n')
        f.write('!entry.%s.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x\n 0 0 0 0 0 0\n'%molname)
        f.write('!entry.%s.unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx\n'%molname)
        f.write(' "%s" 1 %i 1 "?" 0\n'%(molname,len(atomlist)+2))
        f.write('!entry.%s.unit.residuesPdbSequenceNumber array int\n 1\n'%molname)
        f.write('!entry.%s.unit.solventcap array dbl\n -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n'%molname)
        f.write('!entry.%s.unit.velocities table  dbl x  dbl y  dbl z\n')
        for i in range(len(atomlist)):
            f.write(' 0.0 0.0 0.0\n')
    return

#==============================================================================================================================================
# write_amber.mass2pse
# The lib files need the pse number for identification
# This routine finds the pse number based on the mass
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
# amber_write_frcmod
# frcmod file contains all force constant, equilibrium distances,...
# Subroutines:
#    - write_bonds
#    - write_angles
#    - write_dihedrals
#    - write_impropers
#    - write_nonbonded
#==================================================================================================
def write_frcmod(root,distance_unit,energy_unit,polarizability_unit):
    import math
    import xml.etree.ElementTree as ET
    import parameter

    tree = ET.parse('forconx.xml')
    root = tree.getroot()

    f = open('amber_new.frcmod', 'w')
    
    iter = root.getiterator()

    print"\n\t----------------------------------"
    print"\t3.2 Writing frcmod file"
    print"\t----------------------------------"
    print"\t\tFor more information on the file format see http://ambermd.org/formats.html"

    #    output.conk(f,"")
    # Add reference or website for more information about ForConX
    f.write('This frcmod file is generated by ForConX\n')
    f.write('MASS\n')
    atomtype_list = []

    for molecule in root.findall('molecule'):
        for i in molecule.getchildren():
            if i.tag == 'atom':
                if i.get('type') not in atomtype_list:
                    if(str(i.get('polarizability')) == "None"):
                        string = ''.join([i.get('type')," ",i.get('mass'),"\n"])
                    else:
                        string = ''.join([i.get('type')," ",i.get('mass'),"   ", str(i.get('polarizability')),"\n"])
                    atomtype_list.append(i.get('type'))
                    f.write(string)

#   Writing bonds
    print"\t\tWriting <bonds> ..." 
    string = ''.join(["\n\nBOND\n"])
    f.write(string)
    bonds = sorted(parameter.write_bonds(root,energy_unit,distance_unit))
    for line in bonds:
        f.write(line)
    f.write("\n")
    print"\t\t\tNumber of bond potentials = ",len(bonds)
#   Writing angles
    print"\t\tWriting <angles> ..."
    string = ''.join(["\nANGL\n"])
    f.write(string)    
    anglelist = sorted(parameter.write_angles(root,energy_unit))
    for line in anglelist:
        f.write(line)
    f.write("\n")
    print"\t\t\tNumber of angle potentials = ",len(anglelist)
#   Writing dihedrals
    print"\t\tWriting <dihedrals> ..."
    string = ''.join(["\nDIHE\n"])
    f.write(string)
    dihedrals = parameter.write_dihedrals(root,energy_unit)
    for line in dihedrals:
        f.write(line)
    f.write("\n")
    print"\t\t\tNumber of dihedral potentials = ",len(dihedrals)
#   Writing impropers    
    print"\t\tWriting <impropers> ..."
    string = ''.join(["\nIMPR\n"])
    f.write(string)
    impropers = parameter.write_impropers(root,energy_unit)
    for line in impropers:
        f.write(line)
    f.write("\n")
    print"\t\t\tNumber of improper potentials = ",len(impropers)
#   Writing nonbonded
    print"\t\tWriting <nonbonded> ..."
    string = ''.join(["\nNONB\n"])
    f.write(string)
    vdws = sorted(set(parameter.write_vdw(root,energy_unit,distance_unit)))
    for line in vdws:
        f.write(line)
    f.write("\n\n")
    print"\t\t\tNumber of vdw types = ",len(vdws)
    return
