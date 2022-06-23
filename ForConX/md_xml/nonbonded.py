import xml.etree.ElementTree as ET
import helpers 
#===============================================================================================
# Currently the following potentials are supported.
# -lj    V = 4*epsilon ( [sigma/r]^{12} - [sigma/r]^6 )
#===============================================================================================
from ..potentials import lennard_jones 

#===============================================================================================
# nonbonded.average
#===============================================================================================
def average(method,value_i,value_j):
    if method=='geometric':
        try:
            return (value_i*value_j)**0.5
        except:
            helpers.error('Geometric averaging was unsuccessful!')
    elif method=='arithmetic':
        try:
            return (value_i+value_j)*0.5
        except:
            helpers.error('Arithmetic averaging was unsuccessful!')
    else:
        helpers.error('Method %s unknown.'%method)

        
#===============================================================================================
# nonbondedElement                                                                       
#===============================================================================================
class nonbondedElement:
    """
    An XML entry handled by this class should look like this:
    <nonbonded>
       <mixing_rules epsilon="geometric" sigma="arithmetic" />
       <atom type="..." epsilon="..." sigma="..." vdw14="..." elec14="..." />
             ...
       <atom type="..." epsilon="..." sigma="..." vdw14="..." elec14="..." />
       <vdw type="... ..." epsilon="..." sigma="..." vdw14="..." elec14="..." />
             ...
       <vdw type="... ..." epsilon="..." sigma="..." vdw14="..." elec14="..." />
    </nonbonded>

                                    << interface >>
    ______________________________________________________________________________________
    + root                     memory address of root
    + pointer                  memory address of <dihedrals>
    + mixing_epsilon           mixing rule for epsilons
    + mixing_sigma             mixing rule for sigmas
    + elec14_total             total (?= average) value of 1-4 Coulomb scaling factors
    + vdw14_total              total (?= average) value of 1-4 van-der-Waals 
                               scaling factors
    ______________________________________________________________________________________
    - __init__(root)               initializes <nonbonded> object
    - __setattr__(attribute,value) 
    + remove()                     removes <nonbonded> object
    + find_type(type)              returns memory address of the nonbonded with the 
                                   atom types "... ..." or "..."
    + list(tag)                    returns a list of all types of a particular 
                                   subelement of <nonbonded>
    + check()                      checks the XML sanity of the <nonbonded> object
    + epsilon(type_i,type_j)       computes the epsilon_ij according to the mixing rules
    + sigma(type_i,type_j)         computes the sigma_ij according to the mixing rules
    + atom2vdw()                   replaces all <atom> entries in <nonbonded> with
                                   corresponding <vdw>s
    + vdw2atom()                   reduces the number of <vdw> as many as possible by
                                   creating <atom> and relying on mixing rules
    ______________________________________________________________________________________
    """
    def __init__(self,root):
        self.root = root
        self.pointer = self.root.find('nonbonded')
        if self.pointer is None:
            self.pointer = ET.SubElement(self.root,'nonbonded')        
        pointer = self.root.find('nonbonded/mixing_rules')
        if pointer is None:
            ET.SubElement(self.pointer,'mixing_rules')
            self.mixing_epsilon = 'geometric'
            self.mixing_sigma   = 'arithmetic'
        else:
            self.mixing_epsilon = pointer.get('epsilon')
            self.mixing_sigma   = pointer.get('sigma')

    def __setattr__(self,attribute,value):
        if attribute == 'mixing_epsilon':
            value = (value.lower()).strip()
            self.__dict__[attribute] = value
            pointer = self.root.find('nonbonded/mixing_rules')
            pointer.set('epsilon',str(value))
        elif attribute == 'mixing_sigma':
            value = (value.lower()).strip()
            self.__dict__[attribute] = value
            pointer = self.root.find('nonbonded/mixing_rules')
            pointer.set('sigma',str(value))
        elif attribute in 'root pointer':
            self.__dict__[attribute] = value
        elif attribute in 'elec14 vdw14 elec14_total vdw14_total':
            self.__dict__[attribute] = value
            
    def remove(self):
        self.pointer.clear()
        del self

    def list(self,tag):
        list = []
        if tag =="14":
            import molecule
            for mol in self.root.findall('molecule'):
                for dihedral in mol.findall('dihedral'):
                    current_dihedral = molecule.dihedralClass(mol,dihedral.get('name'))
                    type = (current_dihedral.type).split()
                    list.append(' '.join([type[0],type[3]]))
        else:
            for i in self.pointer:
                if i.tag.upper() in tag.upper():
                    list.append(i.get("type"))
        return list

    def find_type(self,type):
        try:
            pointer = self.pointer.find('./*/[@type="%s"]'%type)
            return pointer.tag.upper(),pointer
        except:
            return None, None
        
    def epsilon(self,type):
        try:
            pointer_i = self.root.find('nonbonded/atom/[@type="%s"]'%type[0])
            epsilon_i = float(pointer_i.get('epsilon'))
            pointer_j = self.root.find('nonbonded/atom/[@type="%s"]'%type[1])
            epsilon_j = float(pointer_j.get('epsilon'))
            epsilon_ij= average(self.mixing_epsilon,epsilon_i,epsilon_j)
            return epsilon_ij
        except:
            try:
                import bonds
                pointer = self.root.find('nonbonded/vdw/[@type="%s"]'%bonds.sequence(' '.join(type)))
                epsilon_ij = float(pointer.get('epsilon'))
                return epsilon_ij
            except:
                return None
    
    def sigma(self,type):
       try:
            pointer_i = self.root.find('nonbonded/atom/[@type="%s"]'%type[0])
            sigma_i   = float(pointer_i.get('sigma'))
            pointer_j = self.root.find('nonbonded/atom/[@type="%s"]'%type[1])
            sigma_j   = float(pointer_j.get('sigma'))
            sigma_ij  = average(self.mixing_sigma,sigma_i,sigma_j)
            return sigma_ij
       except:
           try:
               import bonds
               pointer    = self.root.find('nonbonded/vdw/[@type="%s"]'%bonds.sequence(' '.join(type)))
               sigma_ij   = float(pointer.get('sigma'))
               return sigma_ij
           except:
               return None

    def elec14(self,type):
        import bonds  
        current_type = bonds.sequence(' '.join(type))
        if current_type not in self.list("14"):
            return None
        try:
            pointer = self.root.find('nonbonded/vdw/[@type="%s"]'%current_type)
            return float(pointer.get('elec14'))
        except:
            atom_1 = atomClass(self.root,type[0])
            atom_4 = atomClass(self.root,type[1])
            elec14_1 = atom_1.elec14
            elec14_4 = atom_4.elec14
            if elec14_1 is None or elec14_4 is None:
                return None
            elif elec14_1 != elec14_4:
                helpers.warning('Atom types %s and %s have different elec14!'%(type[0],type[1]))
                return average('geometric',elec14_1,elec14_4)
            return elec14_1

    def vdw14(self,type):
        import bonds
        current_type = bonds.sequence(' '.join(type))
        if current_type not in self.list("14"):
            return None
        try:
            pointer = self.root.find('nonbonded/vdw/[@type="%s"]'%current_type)
            return float(pointer.get('vdw14'))
        except:
            atom_1 = atomClass(self.root,type[0])
            atom_4 = atomClass(self.root,type[1])
            vdw14_1 = atom_1.vdw14
            vdw14_4 = atom_4.vdw14
            if vdw14_1 is None or vdw14_4 is None:
                return None
            elif vdw14_1 != vdw14_4:
                helpers.warning('Atom types %s and %s have different vdw14!'%(type[0],type[1]))
                return average('geometric',vdw14_1,vdw14_4)
            return vdw14_1

    def elec14_total(self):
        elec14 = []
        for i in self.list("14"):
            type = i.split()
            elec14.append(self.elec14(type))
        import numpy
        if len(elec14)>0:
            return numpy.mean(elec14)
        return 1.0
        
    def vdw14_total(self):
        vdw14  = []
        for i in self.list("14"):
            type = i.split()        
            vdw14.append(self.vdw14(type))
        import numpy
        if len(vdw14)>0:
            return numpy.mean(vdw14)
        return 1.0
        
    def atom2vdw(self):
        import bonds
        for i in self.root.findall('./nonbonded/atom'):
            type_i = i.get('type')
            for j in self.root.findall('./nonbonded/atom'):
                type_j = j.get('type')
                epsilon = self.epsilon([type_i,type_j])
                sigma   = self.sigma([type_i,type_j])
                elec14  = self.elec14([type_i,type_j])
                vdw14   = self.vdw14([type_i,type_j])
                type = bonds.sequence(' '.join([type_i,type_j]))
                vdw = vdwClass(self.root,type)
                vdw.epsilon = epsilon
                vdw.sigma   = sigma

                if elec14 is not None:
                    vdw.elec14 = elec14
                if vdw14 is not None:
                    vdw.vdw14 = vdw14
            self.pointer.remove(i)

    def vdw2atom(self):
#       Computing atomic values from self interaction        
        for i in self.root.findall('./nonbonded/vdw'):
            vdw = vdwClass(self.root,i.get('type'))
            type = i.get('type').split()
            if type[0]==type[1]:
                atom = atomClass(self.root,type[0])
                atom.epsilon = vdw.epsilon 
                atom.sigma   = vdw.sigma
                
#       Remove redundant vdw pairs
        for i in self.root.findall('./nonbonded/atom'):
            type_i = i.get('type')
            for j in self.root.findall('./nonbonded/atom'):
                type_j = j.get('type')
                if type_j<type_i:
                    continue
                epsilon = self.epsilon([type_i,type_j])
                sigma   = self.sigma([type_i,type_j])
                type = ' '.join([type_i,type_j])

                current_vdw  = vdwClass(self.root,type)
                epsilon_xml = current_vdw.epsilon
                sigma_xml   = current_vdw.sigma
                redundant   = True
                if current_vdw.elec14 is not None:
                    redundant = False
                if current_vdw.vdw14 is not None:
                    redundant = False
                try:
                    ratio = epsilon / epsilon_xml
                    if (ratio>1.005) or (ratio<0.995):
                        redundant = False
                    ratio = sigma / sigma_xml
                    if (ratio>1.005) or (ratio<0.995):
                        redundant = False
                except:
                    pass
                if redundant:
                    current_vdw.remove()

    def check(self):
        print "\t<atom> ...\t",
        number = 0
        for i in self.pointer.findall('atom'):
            number += 1
            type = i.get("type")
            atom = atomClass(self.root,type)
            atom.check()
        print "%5s potentials checked."%number
        
        print "\t<vdw>  ...\t",
        number = 0
        for i in self.pointer.findall('vdw'):
            number += 1
            type = i.get("type")
            vdw = vdwClass(self.root,type)
            vdw.check()
        print "%5s potentials checked."%number

        for i in self.root.findall('./nonbonded/*'):
            if i.tag not in 'atom vdw mixing_rules':
                print "?\t<%s>\n"%i.tag
        list14 = self.list("14")
        for i in set(list14):
            type = i.split()
            if self.vdw14(type) is None:
                helpers.error('<vdw14 type="%s"> does not exist!'%i)
            if self.elec14(type) is None:
                helpers.error('<elec14 type="%s"> does not exist!'%i)
        print "\n\t%s 1-4 interactions checked."%len(list14)
                

#===============================================================================================
# nonbonded.atomClass                                                                         
#===============================================================================================
class atomClass(lennard_jones.ljClass):
    """
    An XML entry handled by this class should look like this:
       <atom type="..." epsilon="..." sigma="..." vdw14="..." elec14="..." />

    ______________________________________________________________________________________
    + type       atom types of involved atom
    + sigma      LJ - sigma
    + epsilon    LJ - epsilon
    + elec14     scaling factor for 1-4 Coulomb interactions
    + vdw14      scaling factor for 1-4 van-der-Waals interactions
    + root       memory address of root
    + pointer    memory address of current Lennard-Jones
    ______________________________________________________________________________________
    - __init__(root,element,type)  initializes LJ object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes LJ object
    + check()                      checks the XML-sanity of <nonbonded/atom>
    ______________________________________________________________________________________
    """
    def __init__(self,root,type):
        lennard_jones.ljClass.__init__(self,root,'./nonbonded/atom',type )

    def check(self):
        if self.epsilon is None:
            helpers.error('<atom type="%s"> has no epsilon!'%self.type)
        if self.sigma is None:
            helpers.error('<atom type="%s"> has no sigma!'%self.type)
        

#===============================================================================================
# nonbonded.vdwClass                                                                          
#===============================================================================================
class vdwClass(lennard_jones.ljClass):
    """
    An XML entry handled by this class should look like this:
       <vdw type="... ..." epsilon="..." sigma="..." vdw14="..." elec14="..." />

    ______________________________________________________________________________________
    + type       atom types of involved atoms
    + sigma      LJ - sigma
    + epsilon    LJ - epsilon
    + elec14     scaling factor for 1-4 Coulomb interactions
    + vdw14      scaling factor for 1-4 van-der-Waals interactions
    + root       memory address of root
    + pointer    memory address of current Lennard-Jones
    ______________________________________________________________________________________
    - __init__(root,element,type)  initializes LJ object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes LJ object
    + check()                      checks the XML-sanity of <nonbonded/atom>
    ______________________________________________________________________________________
    """
    def __init__(self,root,type):
        import bonds
        lennard_jones.ljClass.__init__(self,root,'./nonbonded/vdw',bonds.sequence(type) )

    def check(self):
        if self.epsilon is None:
            helpers.error('<vdw type="%s"> has no epsilon!'%self.type)
        if self.sigma is None:
            helpers.error('<vdw type="%s"> has no sigma!'%self.type)
        current_nonbonded = nonbondedElement(self.root)
        if self.type in current_nonbonded.list("14"):
            if self.elec14 is None:
                helpers.error('<vdw type="%s"> has no elec14.'%self.type)    
            if self.vdw14 is None:
                helpers.error('<vdw type="%s"> has no vdw14.'%self.type)    

                
