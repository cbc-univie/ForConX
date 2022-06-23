import xml.etree.ElementTree as ET
#======================================================================================
# xyz.pdb2crd
#======================================================================================
def pdb2crd(root,crdfile):
    """
    """
    from ..md_xml import helpers
    from ..md_xml import molecule
    
    print"\n\t-------------------------------------------------"
    print"\t3.5 Converting coordinate file to CHARMM crd"
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
    print '\n\t<crd file ="%s">' % crdfile
    try:
        f = open(crdfile,'w')
    except IOError:
        helpers.error('%s cannot be opened!'%crdfile)
    f.write('* Converted from forconx.pdb\n')
    f.write('*\n')
    f.write('%10s%5s\n'%(len(atoms),'EXT'))

    polarizable = root.find('./molecule/atom/[@alpha]')
    mol_index  = 0
    atom_index = 0

    # first 4 letters of <molecule/name> for each molecule
    residuelist = {}
    for i in root.findall('./molecule'):
        molname = i.get('name')
        residuelist[molname[0:4]] = molname 
    
    for i in residues:
        if i not in residuelist.keys():
            helpers.error('Molecule mismatch: %s not found in XML!'%i)
        else:
            mol = root.find('./molecule/[@name="%s"]'%residuelist[i])
        mol_index += 1
        
        current_molecule = molecule.moleculeElement(mol,i)
        for atom in current_molecule.list('ATOM VIRTUAL'):
            pdb_atom = atoms.pop(0)
            if pdb_atom != atom:
                helpers.error('Atom mismatch in molecule %s: %s not equal to %s'%(i,pdb_atom,atom))
            atom_index += 1
            current_xyz = coor.pop(0)
            line  = '%10s'      % atom_index
            line += '%10s  '    % mol_index
            line += '%-8s  '    % residuelist[i]
            line += '%-8s'      % atom
            line += '%20.10f'   % float(current_xyz[0])
            line += '%20.10f'   % float(current_xyz[1])
            line += '%20.10f  ' % float(current_xyz[2])
            line += '%-8s  '    % i
            line += '%-8s'      % mol_index
            line += '%20.10f\n' % 0
            f.write(line)
            if polarizable:
                atom_index += 1
                line  = '%10s'      % atom_index
                line += '%10s  '    % mol_index
                line += '%-8s  '    % i
                line += '%-8s'      % ('D'+atom)
                line += '%20.10f'   % float(current_xyz[0])
                line += '%20.10f'   % float(current_xyz[1])
                line += '%20.10f  ' % float(current_xyz[2])
                line += '%-8s  '    % i
                line += '%-8s'      % mol_index
                line += '%20.10f\n' % 0
                f.write(line)
    print '\n'
    
#======================================================================================
# xyz.crd2pdb
#======================================================================================
def crd2pdb(root,crdfile):
    """
    """
    from ..md_xml import helpers
    from ..md_xml import molecule
    
    print"\n\t-------------------------------------------------"
    print"\t1.5 Converting coordinate file to XML"
    print"\t-------------------------------------------------"
    print '\t< crd file ="%s">' % crdfile
    try:
        f = open(crdfile,'r')
    except IOError:
        helpers.error('File not found!')

    # reading header
    extended = False
    print '\t\t',
    while True:
        line = f.readline()
        if not line[0:1]=="*":
            # extended file format
            if 'EXT' in line:
                extended = True
            break
    # reading coordinates
    coor = []
    last_index = 0
    atoms = []
    residues = []    
    for line in f:
        if extended:
            current_index   = int(line[11:20])
            current_residue = line[21:29].strip()
            current_atom    = line[32:40].strip()
            current_xyz     = [ line[41:60],line[61:80],line[81:100] ]
        else:
            current_index   = int(line[6:10])
            current_residue = line[11:15].strip()
            current_atom    = line[16:20].strip()
            current_xyz     = [ line[21:30],line[31:40],line[41:50] ]

        # not storing Drude particles
        if current_atom[0:1]=="D":
                continue
        atoms.append(current_atom)
        coor.append(current_xyz)
        
        # new residue?
        if current_index != last_index:
            last_index = current_index
            residues.append(current_residue)
            print current_residue,', ',
            if (len(residues)%10==0):
                print '\n\t\t',
    f.close()
    print '\n'
    
    # Writing pdb coordinate file
    from ..md_xml import molecule
    pdbfile = 'forconx.pdb'
    print '> Printing pdb file %s ...' % pdbfile 
    try:
        f = open(pdbfile,'w')
    except IOError:
        helpers.error('PDB file cannot be opened!')
    f.write('REMARK   Converted from %s\n'%crdfile)
    mol_index  = 0
    atom_index = 0
    for i in residues:
        mol = root.find('./molecule/[@name="%s"]'%i)
        if mol is None:
            helpers.error('Molecule mismatch: %s not found in XML!'%i)
        mol_index += 1
        
        current_molecule = molecule.moleculeElement(mol,i)
        for atom in current_molecule.list('ATOM VIRTUAL'):
            crd_atom = atoms.pop(0)
            if crd_atom != atom:
                helpers.error('Atom mismatch in molecule %s: %s not equal to %s'%(i,crd_atom,atom))
            atom_index += 1
            current_xyz = coor.pop(0)
            line  = 'ATOM  '
            line += '%5i '    % atom_index
            line += '%4s '    % atom
            line += '%4s '    % i
            line += '%4s    ' % mol_index
            line += '%8.3f'   % float(current_xyz[0])
            line += '%8.3f'   % float(current_xyz[1])
            line += '%8.3f'   % float(current_xyz[2])
            line += '  0.00  0.00\n'
            f.write(line)
    print '\n'
