import xml.etree.ElementTree as ET
#===============================================================================================
# helpers.check_xml
#===============================================================================================
def check_xml(root):
    print"\n\t----------------------------------"
    print"\t2.1 Checking <molecule> ..."
    print"\t----------------------------------"
    import molecule
    for mol in root.findall('molecule'):
        molname = mol.get('name')
        print '\t<molecule name="%s"> ...' % molname
        md = molecule.moleculeElement(mol,molname)
        polarizable = md.check()
        print '\t</molecule> ',
        if polarizable is not None:
            print 'is polarizable.',
        print '\n'

    print"\t----------------------------------"
    print"\t2.2 Checking <bonds> ..."
    print"\t----------------------------------"
    import bonds
    current_bonds = bonds.bondsElement(root)
    current_bonds.check()

    print"\n\t----------------------------------"
    print"\t2.3 Checking <angles> ..."
    print"\t----------------------------------"
    import angles
    current_angles = angles.anglesElement(root)
    current_angles.check()

    print"\n\t----------------------------------"
    print"\t2.4 Checking <dihedrals> ..."
    print"\t----------------------------------"
    import dihedrals 
    current_dihedrals = dihedrals.dihedralsElement(root)
    current_dihedrals.check()

    print"\n\t----------------------------------"
    print"\t2.5 Checking <impropers> ..."
    print"\t----------------------------------"
    import impropers 
    current_impropers = impropers.impropersElement(root)
    current_impropers.check()
    
    print"\n\t----------------------------------"
    print"\t2.6 Checking <nonbonded> ..."
    print"\t----------------------------------"
    import nonbonded 
    current_nonbonded = nonbonded.nonbondedElement(root)
    current_nonbonded.check()
    print
    
#===============================================================================================
# helpers.warning
#===============================================================================================
def warning(message):
    print '\n!~~~~~~~~~!'
    print '! Warning ! %s' %message
    print '!~~~~~~~~~!'

#===============================================================================================
# helpers.error
#===============================================================================================
def error(message):
    import sys
    print '\n!~~~~~~~!'
    print '! Error !',message
    print '!~~~~~~~!'
    sys.exit()

#======================================================================================
# helpers.write_xml
#======================================================================================
def write_xml(root,filename):
    """
    prints the current XML-tree.
    """
    import re
    from xml.dom import minidom
    print'> Printing current xml-tree',filename,' ...'
    f=open(filename,'w')
    xmltree=re.sub( '\s+', ' ', ET.tostring(root) ).strip().replace('> <','><')
    f.write(minidom.parseString(xmltree).toprettyxml() )
    f.close()
    print

    
