import xml.etree.ElementTree as ET
#======================================================================================
# xyz.pdb2config
#======================================================================================
def pdb2config(root,configfile):
    """
    """
    from ..md_xml import helpers
    from ..md_xml import molecule
    
    print"\n\t-------------------------------------------------"
    print"\t3.3 Converting coordinate file to DLPOLY CONFIG"
    print"\t-------------------------------------------------"
    print '\t<pdb file ="forconx.pdb">'
    try:
        f = open('forconx.pdb','r')
    except IOError:
        helpers.error('forconx.pdb was not found!')
    coor = []
    atoms = []
    residues = []
    last_index = 0
    for line in f:
        if line[0:4] == 'ATOM':
            atoms.append(line[12:16].strip())
            current_index = int(line[23:26])
            if current_index != last_index:
                last_index = current_index
                residues.append(line[17:22].strip())

            coor.append( [ line[31:38], line[39:46], line[47:54] ] ) 
    f.close()
    print '\t\tNumber of atoms    = ',len(atoms)
    print '\t\tNumber of residues = ',len(residues)
    print '\n\t<config file ="%s">' % configfile
    try:
        f = open(configfile,'w')
    except IOError:
        helpers.error('%s cannot be opened!'%configfile)
    f.write('Converted from forconx.pdb\n')
    f.write('%8s%8s\n'%(0,0))
    print '\t\tPlease note, that no periodic boundary conditions are applied!'
    print '\t\tAll virtual atoms are removed since current ForConX version does not support the conversion to DLPOLY.'

    polarizable = root.find('./molecule/atom/[@alpha]')
    mol_index  = 0
    atom_index = 0
    for i in residues:
        mol = root.find('./molecule/[@name="%s"]'%i)
        if mol is None:
            helpers.error('Molecule mismatch: %s not found in XML!'%i)
        mol_index += 1
        
        current_molecule = molecule.moleculeElement(mol,i)
        current_atomlist = current_molecule.list('ATOM')
        current_virtuallist = current_molecule.list('VIRTUAL')
        for atom in current_atomlist:
            pdb_atom = atoms.pop(0)
            if pdb_atom in current_virtuallist:
                continue
            if pdb_atom != atom:
                helpers.error('Atom mismatch in molecule %s: %s not equal to %s'%(i,pdb_atom,atom))
            atom_index += 1
            current_atom = molecule.atomClass(mol,atom)
            current_xyz = coor.pop(0)
            line  = '%-8s%8s\n' %(current_atom.type,atom_index)
            line += '%20.10f'   % float(current_xyz[0])
            line += '%20.10f'   % float(current_xyz[1])
            line += '%20.10f\n' % float(current_xyz[2])
            f.write(line)
            if polarizable:
                atom_index += 1
                line  = '%-8s%8s\n' %('D'+atom,atom_index)
                line += '%20.10f'   % float(current_xyz[0])
                line += '%20.10f'   % float(current_xyz[1])
                line += '%20.10f\n' % float(current_xyz[2])
                f.write(line)
    print '\n'
                
#======================================================================================
# xyz.config2pdb
#======================================================================================
def config2pdb(root,configfile):
    """
    """
    from ..md_xml import helpers
    from ..md_xml import molecule
    
    print"\n\t-------------------------------------------------"
    print"\t1.5 Converting coordinate file to XML"
    print"\t-------------------------------------------------"
    helpers.write_xml(root,'forconx.xml')
    print '\t< config file ="%s">' % configfile
    try:
        f = open(configfile,'r')
    except IOError:
        helpers.error('File not found!')

    # header and box information
    f.readline()
    line_element = (f.readline()).split()
    levcfg = int(line_element[0])
    imcon  = int(line_element[1])
    if imcon >0:
        for i in range(3):
            line = f.readline()
            line_element = line.split()
            for j in range(3):
                if j!=i and abs(float(line_element[j])) >0.001:
                    helpers.warning('No cubic box!')
    coor = []
    atoms = []
    for line in f:
        line_element = line.split()
        atoms.append(line_element[0])

        # xyz
        line = f.next()
        line_element = line.split()
        coor.append( line_element[0:3] )
        if levcfg >0:
            # velocities
            line = f.next()
            if levcfg >1:
                # forces
                line = f.next()
    f.close()
    print '\n'
    
    # Writing pdb coordinate file
    pdbfile = 'forconx.pdb'
    print '> Printing pdb file %s ...' % pdbfile 
    try:
        f = open(pdbfile,'w')
    except IOError:
        helpers.error('PDB file cannot be opened!')
    f.write('REMARK   Converted from %s\n'%configfile)
    mol_index  = 0
    atom_index = 0
    for mol in root.findall('./molecule'):
        molname = mol.get('name')
        current_molecule = molecule.moleculeElement(mol,molname)
        try:
            nmol = int(mol.get('nmol'))
        except:
            helpers.error('Number of molecules of <molecule name="%s"> is not defined!'%molname)
        
        for i in range(nmol):
            mol_index +=1
            for atom in current_molecule.list('ATOM VIRTUAL'):
                current_atom = molecule.atomClass(mol,atom)
                crd_atom = atoms.pop(0)
                if crd_atom != current_atom.type:
                    helpers.error('Atom mismatch in <molecule name="%s">: %s not equal to %s'%(molname,crd_atom,current_atom.type))
                atom_index += 1
                current_xyz = coor.pop(0)
                line  = 'ATOM  '
                line += '%5i '    % atom_index
                line += '%4s '    % current_atom.name
                line += '%4s '    % molname[0:4]
                line += '%4s    ' % mol_index
                line += '%8.3f'   % float(current_xyz[0])
                line += '%8.3f'   % float(current_xyz[1])
                line += '%8.3f'   % float(current_xyz[2])
                line += '  0.00  0.00\n'
                f.write(line)
    print '\n'    

