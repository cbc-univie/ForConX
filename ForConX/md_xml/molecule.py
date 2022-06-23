from xml.etree   import ElementTree as ET
#===============================================================================================
# molecule.moleculeElement                                                                         
#===============================================================================================
class moleculeElement:
    """
    This object handles <molecule>:
      <molecule name="..." nmol="..." >
        <atom     ... />
        <virtual  ... />
        <bond     ... />
        <angle    ... />
        <dihedral ... />
        <improper ... />

        <coordinate name="..." imol="..." xyz="... ... ..." />
      </molecule>

                               << interface >>
    ______________________________________________________________________________________
    + mol       memory address of the corresponding molecule
    + nmol      number of the current molecular species
    + name      name of the current molecule
    + charge    molecular charge of the current molecule
    + mass      molecular mass of the current molecule
    ______________________________________________________________________________________
    + unique_atomname(atomname)           makes the atomname unique
    + adjacency_list()                    creates connectivity matrix
    + generate_angle()                    autogenerates all angles
    + generate_dihedral()                 autogenerates all dihedrals
    + find_path(first_atom,last_atom)     returns the shortest path between two atoms
    + find_name(name)                     returns memory address of the subelement defined by name
    + list(tag)                           returns all subelements with name tag
    + check()                             checks the XML sanity of moleculeElement
    ______________________________________________________________________________________
    """
    def __init__(self,mol,name):
        self.mol = mol
        self.name = name
        try:
            self.nmol = int(self.mol.get('nmol'))
        except:
            self.nmol = 1
        self.charge = 0.0
        self.mass   = 0.0
        for atom in mol.findall('atom'):
            self.charge += float(atom.get('charge'))
            self.mass   += float(atom.get('mass'))
        for virtual in mol.findall('virtual'):
            self.charge += float(virtual.get('charge'))
            
    def remove(self):
        self.mol.clear()
        del self

    def unique_atomname(self,atomname):
        atomlist = self.list("ATOM VIRTUAL") 
        while atomname in atomlist:
            last_letter = atomname[-1]
            ascii = ord(last_letter)
            if (len(atomname)==1) or (ascii==90):
                atomname += "1"
                continue
            if ascii==57:
                ascii=64
            atomname=atomname[:-1]+chr(ascii+1)
        return atomname
        
    def adjacency_list(self):
        adjacency = {}
        for i in self.mol.findall('bond'):
            atom = (i.get('name')).split()
            if not adjacency.has_key(atom[0]):
                adjacency[atom[0]] = []
            if not adjacency.has_key(atom[1]):
                adjacency[atom[1]] = []
            adjacency[atom[0]].append(atom[1])
            adjacency[atom[1]].append(atom[0])
        return adjacency
            
    def generate_angle(self):
        import angles
        list = []
        adjacency = self.adjacency_list()
        for i in adjacency.keys():
            for j in adjacency[i]:
                for k in adjacency[j]:
                    if k == i:
                        continue
                    list.append(angles.sequence(' '.join([i,j,k])))
        return sorted(set(list))

    def generate_dihedral(self):
        import dihedrals
        list = []
        adjacency = self.adjacency_list()
        for i in adjacency.keys():
            for j in adjacency[i]:
                for k in adjacency[j]:
                    if k == i:
                        continue
                    for l in adjacency[k]:
                        if l == j or l == i:
                            continue
                        list.append(dihedrals.sequence(' '.join([i,j,k,l])))
        return sorted(set(list))

    def generate_improper(self):
        import impropers
        list = []
        adjacency = self.adjacency_list()
        for i in adjacency.keys():
            for j in adjacency[i]:
                for k in adjacency[i]:
                    if k == j:
                        continue
                    for l in adjacency[i]:
                        if l == j or l == k:
                            continue
                        list.append(impropers.sequence(' '.join([i,j,k,l])))
        return sorted(set(list))
    
    def find_path(self, first_atom, last_atom, path=[]):
        adjacency = self.adjacency_list()
        path = path + [first_atom]
        if first_atom == last_atom:
            return path
        if not adjacency.has_key(first_atom):
            return None
        shortest = None
        for node in adjacency[first_atom]:
            if node not in path:
                newpath = self.find_path(node,last_atom,path)
                if newpath:
                    if not shortest or len(newpath) < len(shortest):
                        shortest = newpath
        return shortest

    def find_name(self,name):
        pointer = self.mol.find('./*/[@name="%s"]'%name)
        try:
            return pointer.tag.upper(), pointer
        except:
            return None, None
    
    def list(self,tag):
        list = []
        for i in self.mol:
            if i.tag.upper() in tag.upper():
                list.append(i.get('name'))
        return list

    def check(self):
        # print additional tags
        for i in self.mol.findall('*'):
            if i.tag not in 'atom virtual bond angle dihedral improper coordinates':
                print '?\t\t<%s>\n'%i.tag

        print '\t\tAtoms ...        ',
        number = 0
        for i in self.mol.findall('atom'):
            number += 1
            name = i.get('name')
            atom = atomClass(self.mol,name)
            atom.check()
        print '%5s'%number

        print '\t\tVirtual atoms ...',
        number = 0
        for i in self.mol.findall('virtual'):
            number += 1
            name = i.get('name')
            virtual = virtualClass(self.mol,name)
            virtual.check()
        print '%5s'%number

        print '\t\tBonds ...        ',
        number = 0
        for i in self.mol.findall('bond'):
            number += 1
            name = i.get('name')
            bond = bondClass(self.mol,name)
            bond.check()
        print '%5s'%number

        print '\t\tAngles ...       ',
        number = 0
        for i in self.mol.findall('angle'):
            number += 1
            name = i.get('name')
            angle = angleClass(self.mol,name)
            angle.check()
        print '%5s'%number

        print '\t\tDihedrals ...    ',
        number = 0
        for i in self.mol.findall('dihedral'):
            number += 1
            name = i.get('name')
            dihedral = dihedralClass(self.mol,name)
            dihedral.check()
        print '%5s'%number

        print '\t\tImpropers ...    ',
        number = 0
        for i in self.mol.findall('improper'):
            number += 1
            name = i.get('name')
            improper = improperClass(self.mol,name)
            improper.check()
        print '%5s'%number
        
        
#===============================================================================================
# molecule.atomClass                                                                         
#===============================================================================================
class atomClass:
    """
    This object handles <atom>-SubElements of <molecule>:
      <atom name="..." type="..." charge="..." alpha="..." mass="..." />

    ______________________________________________________________________________________
    + mol       memory address of the corresponding molecule
    + pointer   memory address of the current atom
    + name      name of the current atom
    + type      type of the current atom
    + charge    partial charge of the current atom
    + mass      mass of the current atom
    + alpha     polarizability if this atom is polarizable)
    ______________________________________________________________________________________
    - __init__(mol,name)           initializes atomClass object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes atomClass object (including XML entry)
    + check()                      checks the XML sanity of the atomClass object 
    ______________________________________________________________________________________
    """
    def __init__(self,mol,name):
        self.mol     = mol 
        self.pointer = self.mol.find('./atom/[@name="%s"]'%name)
        if self.pointer is None:
            self.pointer = ET.SubElement(mol,'atom')
            self.type   = None
            self.mass   = None
            self.charge = None
            self.alpha  = None
        else:
            self.type   = self.get('type')
            self.mass   = self.get('mass')
            self.charge = self.get('charge')
            self.alpha  = self.get('alpha')
        self.name = name
 
    def __setattr__(self,attribute,value):
        if attribute in 'pointer mol':
            self.__dict__[attribute] = value
            
        elif attribute in 'name type':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.pointer.set(attribute,value)
            except:
                self.__dict__[attribute] = None
                
        elif attribute in 'mass':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,3)))

        elif attribute in 'alpha':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,6)))
                
        elif attribute in 'charge':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,6)))

    def get(self,attribute):
        try:
            value = self.pointer.get(attribute)
            if attribute in 'mass charge alpha':
                value = float(value)
        except:
            value = None
        return value

    def remove(self):
        self.pointer.clear()
        self.mol.remove(self.pointer)
        self.__dict__[pointer] = None
        del self
        
    def check(self):
        import helpers
        if not self.pointer.get("type"):
            helpers.error('<atom name="%s"> has no atomtype!'%self.name)
        if not self.pointer.get("charge"):
            helpers.error('<atom name="%s"> has no charge!'%self.name)
        if not self.pointer.get("mass"):
            helpers.error('<atom name="%s"> has no mass!'%self.name)

            
#===============================================================================================
# molecule.virtualClass                                                                         
#===============================================================================================
class virtualClass:
    """
    This object handles <virtual>-SubElements of <molecule>:
      <virtual name="..." type="..." charge="..." 
               zmatrix="... ... ..." r="..." theta="..." phi="..." />

    ______________________________________________________________________________________
    + mol       memory address of the corresponding molecule
    + pointer   memory address of the current virtual atom
    + name      name of the current virtual atom
    + type      type of the current virtual atom
    + charge    partial charge of the current virtual atom
    + zmatrix   sequence of real atom names to define position of current virtual
    + r         distance between virtual--1
    + theta     angle between    virtual--1--2
    + phi       dihedral between virtual--1--2--3
    ______________________________________________________________________________________
    - __init__(mol,name) initializes virtualClass object
    - __setattr__()      overloads "=" operation
    + remove()           removes virtualClass object (including XML entry)
    + check()            checks the XML sanity of the virtualClass object 
    ______________________________________________________________________________________
    """
    def __init__(self,mol,name):
        self.mol = mol 
        self.pointer = self.mol.find('./virtual/[@name="%s"]'%name)
        if self.pointer is None:
            self.pointer = ET.SubElement(mol,'virtual')
            self.type    = None
            self.charge  = None
            self.zmatrix = None
            self.dist    = None
            self.theta   = None
            self.phi     = None
        else:
            self.type    = self.get('type')
            self.charge  = self.get('charge')
            self.zmatrix = self.get('zmatrix')
            self.dist    = self.get('dist')
            self.theta   = self.get('theta')
            self.phi     = self.get('phi')
        self.name = name 

    def __setattr__(self,attribute,value):
        if attribute in 'pointer mol':
            self.__dict__[attribute] = value
            
        elif attribute in 'name type zmatrix':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.pointer.set(attribute,value)
            except:
                self.__dict__[attribute] = None
                
        elif attribute in 'charge dist theta phi':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,5)))
            
    def get(self,attribute):
        try:
            value = self.pointer.get(attribute)
            if attribute in 'charge dist theta phi':
                value = float(value)
        except:
            value = None
        return value
        
    def remove(self):
        self.pointer.clear()
        self.mol.remove(self.pointer)
        self.__dict__[pointer] = None
        del self

    def check(self):
        import helpers
        if not self.pointer.get("type"):
            helpers.error('<virtual name="%s"> has no atomtype!'%self.name)
        if not self.pointer.get("zmatrix"):
            helpers.error('<virtual name="%s"> has no zmatrix!'%self.name)
        else:
            zmatrix = self.zmatrix.split()
            if len(zmatrix) != 3:
                helpers.error('<virtual name="%s"> has incorrect number of zmatrix atoms!'%self.name)
        if not self.pointer.get("dist"):
            helpers.error('<virtual name="%s"> has no distance dist!'%self.name)
        if not self.pointer.get("theta"):
            helpers.error('<virtual name="%s"> has no angle theta!'%self.name)
        if not self.pointer.get("phi"):
            helpers.error('<virtual name="%s"> has no dihedral phi!'%self.name)


#===============================================================================================
# molecule.bondClass                                                                         
#===============================================================================================
class bondClass:
    """
    This object handles <bond>-SubElements of <molecule>:
      <bond name="... ..."  />

    ______________________________________________________________________________________
    + mol       memory address of the corresponding molecule
    + pointer   memory address of the current bond
    + name      atom names of the current bond
    + type      atom types of the current bond
    ______________________________________________________________________________________
    - __init__(mol,name) initializes bondClass object
    - __setattr__()      overloads "=" operation
    + remove()           removes bondClass object (including XML entry)
    + check()            checks the XML sanity of the bondClass object 
    ______________________________________________________________________________________
    """
    def __init__(self,mol,name):
        self.mol = mol 
        self.pointer = self.mol.find('./bond/[@name="%s"]'%name)
        if self.pointer is None:
            self.pointer = ET.SubElement(self.mol,'bond')
            self.type = None
        else:
            import bonds
            atoms  = name.split()
            atom_1 = atomClass(self.mol,atoms[0])
            atom_2 = atomClass(self.mol,atoms[1])
            self.type = bonds.sequence(' '.join([atom_1.type,atom_2.type]))
        self.name = name 

    def __setattr__(self,attribute,value):
        if attribute in 'pointer mol':
            self.__dict__[attribute] = value
            
        elif attribute in 'name type':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.pointer.set(attribute,value)
            except:
                self.__dict__[attribute] = None
                
    def remove(self):
        self.pointer.clear()
        self.mol.remove(self.pointer)
        self.__dict__[pointer] = None
        del self

    def check(self):
        import helpers
        name = self.name.split()
        if len(name) != 2:
            helpers.error('<bond name="%s"> has incorrect number of atom names!'%self.name)


#===============================================================================================
# molecule.angleClass                                                                         
#===============================================================================================
class angleClass:
    """
    This object handles <angle>-SubElements of <molecule>:
      <angle name="... ... ..."  />

    ______________________________________________________________________________________
    + mol       memory address of the corresponding molecule
    + pointer   memory address of the current angle
    + name      atom names of the current angle
    + type      atom types of the current angle
    ______________________________________________________________________________________
    - __init__(mol,name) initializes angleClass object
    - __setattr__()      overloads "=" operation
    + remove()           removes angleClass object (including XML entry)
    + check()            checks the XML sanity of the angleClass object 
    ______________________________________________________________________________________
    """
    def __init__(self,mol,name):
        self.mol = mol 
        self.pointer = self.mol.find('./angle/[@name="%s"]'%name)
        if self.pointer is None:
            self.pointer = ET.SubElement(self.mol,'angle')
            self.type = None
        else:
            import angles
            atoms  = name.split()
            atom_1 = atomClass(self.mol,atoms[0])
            atom_2 = atomClass(self.mol,atoms[1])
            atom_3 = atomClass(self.mol,atoms[2])
            self.type = angles.sequence(' '.join([atom_1.type,atom_2.type,atom_3.type]))
        self.name = name 

    def __setattr__(self,attribute,value):
        if attribute in 'pointer mol':
            self.__dict__[attribute] = value
            
        elif attribute in 'name type':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.pointer.set(attribute,value)
            except:
                self.__dict__[attribute] = None
                
    def remove(self):
        self.pointer.clear()
        self.mol.remove(self.pointer)
        self.__dict__[pointer] = None
        del self

    def check(self):
        import helpers
        name = self.name.split()
        if len(name) != 3:
            helpers.error('<angle name="%s"> has incorrect number of atom names!'%self.name)

        current_molecule = moleculeElement(self.mol,self.mol.get('name'))
        for i in range(1):
            path = current_molecule.find_path(name[i],name[i+1])
            if not path or len(path) != 2:
                helpers.error('<angle name="%s"> has no bond between %s and %s!'%(self.name,name[i],name[i+1]))
            

#===============================================================================================
# molecule.dihedralClass                                                                         
#===============================================================================================
class dihedralClass:
    """
    This object handles <dihedral>-SubElements of <molecule>:
      <dihedral name="... ... ..."  />

    ______________________________________________________________________________________
    + mol       memory address of the corresponding molecule
    + pointer   memory address of the current dihedral
    + name      atom names of the current dihedral
    + type      atom types of the current dihedral
    ______________________________________________________________________________________
    - __init__(mol,name) initializes dihedralClass object
    - __setattr__()      overloads "=" operation
    + remove()           removes dihedralClass object (including XML entry)
    + check()            checks the XML sanity of the dihedralClass object 
    ______________________________________________________________________________________
    """
    def __init__(self,mol,name):
        self.mol = mol 
        self.pointer = self.mol.find('./dihedral/[@name="%s"]'%name)
        if self.pointer is None:
            self.pointer = ET.SubElement(mol,'dihedral')
            self.type = None
        else:
            import dihedrals
            atoms  = name.split()
            atom_1 = atomClass(self.mol,atoms[0])
            atom_2 = atomClass(self.mol,atoms[1])
            atom_3 = atomClass(self.mol,atoms[2])
            atom_4 = atomClass(self.mol,atoms[3])
            self.type = dihedrals.sequence(' '.join([atom_1.type,atom_2.type,atom_3.type,atom_4.type]))
        self.name = name 

    def __setattr__(self,attribute,value):
        if attribute in 'pointer mol linear':
            self.__dict__[attribute] = value            
        elif attribute in 'name type':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.pointer.set(attribute,value)
            except:
                self.__dict__[attribute] = None
                
    def remove(self):
        self.pointer.clear()
        self.mol.remove(self.pointer)
        self.__dict__[pointer] = None
        del self

    def check(self):
        import helpers
        name = self.name.split()
        if len(name) != 4:
            helpers.error('<dihedral name="%s"> has incorrect number of atom names!'%self.name)

        current_molecule = moleculeElement(self.mol,self.mol.get('name'))
        for i in range(2):
            path = current_molecule.find_path(name[i],name[i+1])
            if not path or len(path) != 2:
                helpers.error('<dihedral name="%s"> has no bond between %s and %s!'%(self.name,name[i],name[i+1]))


#===============================================================================================
# molecule.improperClass                                                                         
#===============================================================================================
class improperClass:
    """
    This object handles <improper>-SubElements of <molecule>:
      <improper name="... ... ..." central="..."  />

    ______________________________________________________________________________________
    + mol       memory address of the corresponding molecule
    + pointer   memory address of the current improper
    + name      atom names of the current improper
    + type      atom types of the current improper
    + central   name of the central atom of the current improper
    ______________________________________________________________________________________
    - __init__(mol,name) initializes improperClass object
    - __setattr__()      overloads "=" operation
    + remove()           removes improperClass object (including XML entry)
    + check()            checks the XML sanity of the improperClass object 
    ______________________________________________________________________________________
    """
    def __init__(self,mol,name):
        self.mol = mol 
        self.pointer = self.mol.find('./improper/[@name="%s"]'%name)
        if self.pointer is None:
            self.pointer = ET.SubElement(mol,'improper')
            self.type    = None
            self.central = None
        else:
            import impropers
            atoms  = name.split()
            atom_1 = atomClass(self.mol,atoms[0])
            atom_2 = atomClass(self.mol,atoms[1])
            atom_3 = atomClass(self.mol,atoms[2])
            atom_4 = atomClass(self.mol,atoms[3])
            self.type = impropers.sequence(' '.join([atom_1.type,atom_2.type,atom_3.type,atom_4.type]))
            self.central = self.pointer.get('central')
        self.name = name 
            
    def __setattr__(self,attribute,value):
        if attribute in 'pointer mol':
            self.__dict__[attribute] = value
            
        elif attribute in 'name type central':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.pointer.set(attribute,value)
            except:
                self.__dict__[attribute] = None
                
    def remove(self):
        self.pointer.clear()
        self.mol.remove(self.pointer)
        self.__dict__[pointer] = None
        del self

    def check(self):
        import helpers
        name = self.name.split()
        if len(name) != 4:
            helpers.error('<improper name="%s"> has incorrect number of atom names!'%self.name)

        if not self.central:
            helpers.error('<improper name="%s"> has no central atom!'%self.name)
            
        current_molecule = moleculeElement(self.mol,self.mol.get('name'))
        for i in range(3):
            if name[i]==self.central:
                continue
            path = current_molecule.find_path(self.central,name[i])
            if not path or len(path) != 2:
                helpers.error('<improper name="%s"> has no bond between %s and %s!'%(self.name,name[i],name[i+1]))



