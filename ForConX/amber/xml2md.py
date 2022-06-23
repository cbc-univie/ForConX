import xml.etree.ElementTree as ET
import sys

#==================================================================================================
# amber_write_lib
# Lib files contain the topology of each molecule and partial charges
# Subroutine:
#    - mass2pse
#==================================================================================================
def write_lib(root):
    from ..md_xml import molecule
    
    print"\n\t----------------------------------"
    print"\t3.1 Writing lib files"
    print"\t----------------------------------"
    print"\t\tFor more information on the file format see http://ambermd.org/formats.html"
    print"\t\tIf you want to connect two residues you have to modify section residueconnect in the lib file.\n"

    for mol in root.findall('molecule'):
        molname = mol.get('name')
        current_molecule = molecule.moleculeElement(mol,molname)
        filename = ''.join([molname,".lib"])

        print '\n\t\t<molecule ="%s">   --> %s' % (molname,filename)
        f = open(filename,'w')
        f.write('!!index array str\n')
        f.write(' "%s"\n'%molname)
        # Atomtype definition per molecule - charges are assigned
        f.write('!entry.%s.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n'%molname)

        index = 1
        indices = {}
        print "\t\t\t",
        for i in current_molecule.list('ATOM'):
            print "%5s, " %i,
            if (index%5==0):
                print "\n\t\t\t",
            
            current_atom = molecule.atomClass(mol,i)
            pse_number = mass2pse(current_atom.mass)
            line  = ' "%s"' % current_atom.name
            line += ' "%s"' % current_atom.type
            line += ' 0 1  17956867 %i %i ' % (index, pse_number)
            line += '%s\n' % current_atom.charge
            f.write(line)
            indices[current_atom.name] = index
            index +=1
            
        print
        f.write('!entry.%s.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n'%molname)
        for i in current_molecule.list('ATOM'):
            current_atom = molecule.atomClass(mol,i)
            line  = ' "%s"' % current_atom.name
            line += ' "%s"' % current_atom.type
            line += ' 0 -1 0.0\n'
            f.write(line)
        f.write('!entry.%s.unit.boundbox array dbl\n -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n'%molname)
        f.write('!entry.%s.unit.childsequence single int\n 2\n'%molname)
        f.write('!entry.%s.unit.connect array int\n 0\n 0\n'%molname)
        f.write('!entry.%s.unit.connectivity table  int atom1x  int atom2x  int flags\n'%molname)
        index = 0

        for i in current_molecule.list('BOND'):
            atomnames = i.split()
            line  = ' %i' % indices[atomnames[0]]
            line += ' %i' % indices[atomnames[1]]
            line += ' 1\n'
            f.write(line)
            index+=1
            
        f.write('!entry.%s.unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx\n'%molname)
        index = 1
        f.write(' "U" 0 "R" %i\n'%index)
        natom = len(current_molecule.list('ATOM'))
        for i in range(natom):
            f.write(' "R" 1 "A" %i\n'%index)
            index+=1
        f.write('!entry.%s.unit.name single str\n "%s"\n' %(molname,molname))
        f.write('!entry.%s.unit.positions table  dbl x  dbl y  dbl z\n'%molname)
        for i in range(natom):
            f.write(' 0.0 0.0 0.0\n')
        f.write('!entry.%s.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x\n 0 0 0 0 0 0\n'%molname)
        f.write('!entry.%s.unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx\n'%molname)
        f.write(' "%s" 1 %i 1 "?" 0\n'%(molname,natom+2))
        f.write('!entry.%s.unit.residuesPdbSequenceNumber array int\n 1\n'%molname)
        f.write('!entry.%s.unit.solventcap array dbl\n -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n'%molname)
        f.write('!entry.%s.unit.velocities table  dbl x  dbl y  dbl z\n'%molname)
        for i in range(natom):
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
def write_frcmod(root):
    import math
    import xml.etree.ElementTree as ET
    import parameter

    from ..md_xml import input_output
    current_output = input_output.mdElement(root,'output')
    energy_conversion   = current_output.convert_energy()
    distance_conversion = current_output.convert_distance()
    polarizability_conversion = current_output.convert_polarizability()
 
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
                        f.write("%s  %s\n"%(i.get('type'),i.get('mass')))
                    else:
                        f.write("%s  %s   %s\n"%(i.get('type'),i.get('mass'),str(i.get('polarizability'))))
                    atomtype_list.append(i.get('type'))

#   Writing bonds
    print"\t\tWriting <bonds> ..." 
    f.write('\n\nBOND\n')
    bonds = sorted(parameter.write_bonds(root,energy_conversion,distance_conversion))
    for line in bonds:
        f.write(line)
    print"\t\t\tNumber of bond potentials = ",len(bonds)
#   Writing angles
    print"\t\tWriting <angles> ..."
    f.write('\n\nANGL\n')    
    anglelist = sorted(parameter.write_angles(root,energy_conversion))
    for line in anglelist:
        f.write(line)
    print"\t\t\tNumber of angle potentials = ",len(anglelist)
#   Writing dihedrals
    print"\t\tWriting <dihedrals> ..."
    f.write('\n\nDIHE\n')
    dihedrals = parameter.write_dihedrals(root,energy_conversion)
    for line in dihedrals:
        f.write(line)
    print"\t\t\tNumber of dihedral potentials = ",len(dihedrals)
#   Writing impropers    
    print"\t\tWriting <impropers> ..."
    f.write('\n\nIMPR\n')
    impropers = parameter.write_impropers(root,energy_conversion)
    for line in impropers:
        f.write(line)
    print"\t\t\tNumber of improper potentials = ",len(impropers)
#   Writing nonbonded
    print"\t\tWriting <nonbonded> ..."
    f.write('\n\nNONB\n')
    vdws = sorted(set(parameter.write_vdw(root,energy_conversion,distance_conversion)))
    for line in vdws:
        f.write(line)
    print"\t\t\tNumber of vdw types = ",len(vdws)
    return
