import xml.etree.ElementTree as ET

#=============================================================================================
# Currently supported potentials:
# -harm  V = k  ( r - r0 )^2
# -mors  V = D0 ( 1 - exp[-beta (r-r0)] )^2
#=============================================================================================
from ..potentials import harmonic_bond 
from ..potentials import morse         

#=============================================================================================
# bonds.sequence
#=============================================================================================
def sequence(type):
    """
    sorts alphabetically the bond atom types to facilitate the search for duplicates:
    {A-B, B-A} --> A-B
    """
    tmp = type.split()
    if tmp[0] >tmp[1]:
        tmp = tmp[::-1]
    return ' '.join(tmp)

#=============================================================================================
# bondsElement                                                                          
#=============================================================================================
class bondsElement:
    """    
    An XML entry handled by this class should look like this:
    <bonds>
       <harm type="..." r0="..." k="..." />
             ...
       <harm type="..." r0="..." k="..." />
       <mors type="..." r0="..." D0="..." beta="..." />
             ...
       <mors type="..." r0="..." D0="..." beta="..." />
    </bonds>

                                << interface >>
    ______________________________________________________________________________________
    + root            memory address of root
    + pointer         memory address of <bonds>
    ______________________________________________________________________________________
    - __init__(root)  initializes <bonds> object
    + remove()        removes <bonds> object
    + find_type(type) returns memory address of the bond with the atom types "... ..."
    + list(tag)       returns a list of all types of a particular subelement of <bonds>
    + check()         checks the XML sanity of the <bonds> object
    ______________________________________________________________________________________
    """
    def __init__(self,root):
        self.root = root
        self.pointer = root.find('bonds')
        if self.pointer is None:
            self.pointer = ET.SubElement(root,'bonds')

    def remove(self):
        self.pointer.clear()
        del self
        
    def find_type(self,type):
        pointer = self.pointer.find('./*/[@type="%s"]'%type)
        try:
            return pointer.tag.upper(), pointer
        except:
            return None, None
        
    def list(self,tag):
        list = []
        for i in self.pointer:
            if i.tag.upper() in tag.upper():
                list.append(i.get("type"))
        return list

    def check(self):
        import helpers
        import molecule
        print "\t<molecule> ...",
        number = 0
        for mol in self.root.findall('./molecule'):
            for i in mol.findall('bond'):
                current_bond = molecule.bondClass(mol,i.get('name'))
                existing = False
                number += 1
                if self.find_type(current_bond.type) is not None:
                    existing = True
                if not existing:
                    helpers.error('<bond name="%s"> in molecule %s has no corresponding bond potential'%(i.get("name"),mol.get("name")))
        print "%5s potentials checked."%number
        
        print "\t<harm> ..."
        number = 0
        for i in self.pointer.findall('harm'):
            number += 1
            type = i.get("type")
            harm = harmClass(self.root,type)
            harm.check()
        print "\n\t\t\t%5s potentials checked."%number
        
        print "\t<mors> ..."
        number = 0
        for i in self.pointer.findall('mors'):
            number += 1
            type = i.get("type")
            mors = morsClass(self.root,type)
            mors.check()
        print "\n\t\t\t%5s potentials checked."%number

        for i in self.root.findall('./bonds/*'):
            if i.tag not in 'harm mors':
                print "?\t<%s>\n"%i.tag

                
#=============================================================================================
# bonds.harmClass                                                                             
#=============================================================================================
class harmClass(harmonic_bond.hbondClass):
    """
    harmClass is an object of a particular ./bonds/harm which should look like this:
        <harm k="..." r0="..." type="... ..." />

    ______________________________________________________________________________________
    + k         force constant
    + r0        equilibrium distance
    + type      type of the involved atoms
    + root      memory address of root
    ______________________________________________________________________________________
    - __init__(root,type)          initializes harmClass object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes current harmmClass object (including XML entry)
    + check()                      checks the XML sanity of the harmClass object
    ______________________________________________________________________________________
    """
    def __init__(self,root,type):
        if root.find('bonds') is None:
            ET.SubElement(root,'bonds')
        harmonic_bond.hbondClass.__init__(self,root,'./bonds/harm',sequence(type) )

    def check(self):
        import helpers
        import molecule
        existing = False
        for mol in self.root.findall('molecule'):
            current_molecule = molecule.moleculeElement(mol,mol.get('name'))
            for i in current_molecule.list('bond'):
                name = i.split()
                atom_i = molecule.atomClass(mol,name[0])
                atom_j = molecule.atomClass(mol,name[1])
                if self.type == sequence(' '.join([atom_i.type,atom_j.type])):
                    existing = True
        if not existing:
            helpers.warning('<harm type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return

        if self.k is None:
            helpers.warning('<harm type="%s"> has no force constant! Constant bond?'%self.type)
        if self.r0 is None:
            helpers.error('<harm type="%s"> has no equilibrium distance!'%self.type)

            
#=============================================================================================
# bonds.morsClass                                                                             
#=============================================================================================
class morsClass(morse.morseClass):
    """
    morsClass is an object of a particular ./bonds/mors which should look like this:
        <mors beta="..." D0="..." r0="..." type="... ..." />

    ______________________________________________________________________________________
    + D0        Dissociation energy
    + r0        equilibrium distance
    + beta      stretching parameter
    + type      type of the involved atoms
    + root      memory address of root
    ______________________________________________________________________________________
    - __init__(root,type)          initializes harmClass object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes current harmmClass object (including XML entry)
    + check()                      checks the XML sanity of the harmClass object
    + harm()                       returns harmonic force constants and equilibrium
                                   distance
    ______________________________________________________________________________________
    """
    def __init__(self,root,type):
        if root.find('bonds') is None:
            ET.SubElement(root,'bonds')
        morse.morseClass.__init__(self,root,'./bonds/mors',sequence(type) )

    def harm(self):
        k = self.D0*self.beta*self.beta
        return k, self.r0
        
    def check(self):
        import helpers
        import molecule
        existing = False
        for mol in self.root.findall('molecule'):
            current_molecule = molecule.moleculeElement(mol,mol.get('name'))
            for i in current_molecule.list('bond'):
                name = i.split()
                atom_i = molecule.atomClass(mol,name[0])
                atom_j = molecule.atomClass(mol,name[1])
                type = sequence(' '.join([atom_i.type,atom_j.type]))
                if type == self.type:
                    existing = True
        if not existing:
            helpers.warning('<mors type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return
        if self.k is None:
            helpers.error('<mors type="%s"> has no force constant!'%self.type)
        if self.r0 is None:
            helpers.error('<mors type="%s"> has no equilibrium distance!'%self.type)
        if self.D0 is None:
            helpers.error('<mors type="%s"> has no dissociation energy!'%self.type)
        if self.beta is None:
            helpers.error('<mors type="%s"> has no stretching constant beta!'%self.type)

        
