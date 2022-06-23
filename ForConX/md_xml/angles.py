import xml.etree.ElementTree as ET

#===============================================================================================
# Currently supported potentials:
# -harm  V = k ( theta - theta0 )^2
# -urey  V = k ( r13   - r13_0  )^2
#===============================================================================================
from ..potentials import harmonic_angle 
from ..potentials import harmonic_bond  

#===============================================================================================
# angles.sequence
#===============================================================================================
def sequence(type):
    """
    sorts alphabetically the angle potential types to facilitate the search for duplicates.
    {A-B-C, C-B-A} --> A-B-C
    """
    tmp = type.split()
    if tmp[0] > tmp[2]:
#       reverting sequence
        tmp = tmp[::-1]
    return ' '.join(tmp)


#===============================================================================================
# angles.anglesElement                                                                         
#===============================================================================================
class anglesElement:
    """
    An XML entry handled by this class should look like this:
    <angles>
       <harm type="..." k="..." theta0="..." />
             ...
       <harm type="..." k="..." theta0="..." />
       <urey type="..." k="..." r0="..." />
             ...
       <urey type="..." k="..." r0="..." />
    </angles>


                                << interface >>
    ______________________________________________________________________________________
    + root           memory address of root
    + pointer        memory address of <angles>
    ______________________________________________________________________________________
    - __init__(root) initializes <angles> object
    + remove()       removes <angles> object
    + find_type(type) returns memory address of the angle with the atom types "... ..."
    + list(tag)       returns a list of all types of a particular subelement of <angles>
    + check()         checks the XML sanity of the <angles> object
    ______________________________________________________________________________________
    """
    def __init__(self,root):
        self.root = root
        self.pointer = root.find('angles')
        if self.pointer is None:
            self.pointer = ET.SubElement(root,'angles')

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
        print "\t<molecule> ...\t",
        number = 0
        for mol in self.root.findall('./molecule'):
            for i in mol.findall('angle'):
                current_angle = molecule.angleClass(mol,i.get('name'))
                existing = False
                number += 1
                if self.find_type(current_angle.type) is not None:
                    existing = True
                if not existing:
                    helpers.error('<angle name="%s"> in molecule %s has no corresponding angle potential'%(i.get("name"),mol.get("name")))
        print "%5s potentials checked."%number
        
        print "\t<harm> ..."
        number = 0
        for i in self.pointer.findall('harm'):
            number += 1
            type = i.get("type")
            harm = harmClass(self.root,type)
            harm.check()
        print "\n\t\t\t%5s potentials checked."%number
        
        print "\t<urey> ...\t"
        number = 0
        for i in self.pointer.findall('urey'):
            number += 1
            type = i.get("type")
            urey = ureyClass(self.root,type)
            urey.check()
        print "\n\t\t\t%5s potentials checked."%number

        for i in self.root.findall('./angles/*'):
            if i.tag not in 'harm urey':
                print "?\t<%s>\n"%i.tag

                
#===============================================================================================
# angles.harmClass                                                                            
#===============================================================================================
class harmClass(harmonic_angle.hangleClass):
    """
    harmClass is an object of a particular ./angles/harm which should look like this:
        <harm k="..." theta0="..." type="... ..." />

    ______________________________________________________________________________________
    + k         force constant
    + theta0    equilibrium angle
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
        if root.find('angles') is None:
            ET.SubElement(root,'angles')
        harmonic_angle.hangleClass.__init__(self,root,'./angles/harm',sequence(type) )

    def check(self):
        import helpers 
        import molecule 
        existing = False
        for mol in self.root.findall("molecule"):
            current_molecule = molecule.moleculeElement(mol,mol.get('name'))
            for i in current_molecule.list("angle"):
                name = i.split()
                atom_i = molecule.atomClass(mol,name[0])
                atom_j = molecule.atomClass(mol,name[1])
                atom_k = molecule.atomClass(mol,name[2])
                type = sequence(' '.join([atom_i.type,atom_j.type,atom_k.type]))
                if type == self.type:
                    existing = True
        if not existing:
            helpers.warning('<harm type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return
        if self.k is None:
            helpers.error('<harm type="%s"> has no force constant!'%self.type)
        if self.theta0 is None:
            helpers.error('<harm type="%s"> has no equilibrium angle!'%self.type)

            
#===============================================================================================
# angles.ureyClass                                                                            
#===============================================================================================
class ureyClass(harmonic_bond.hbondClass):
    """
    ureyClass is an object of a particular ./angles/urey which should look like this:
        <urey k="..." r0="..." type="... ..." />

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
        if root.find('angles') is None:
            ET.SubElement(root,'angles')
        harmonic_bond.hbondClass.__init__(self,root,'./angles/urey',sequence(type) )

    def harm(self):
        import bonds
        import molecule
        import numpy

        type = self.type.split()
        
        # looking for participating bonds
        try:
            bond_ab = bonds.harmClass(root,'%s %s'%(type[0],type[1]))
            ab  = bond_ab.r0
            ab2 = bond_ab.r0 * bond_ab.r0
            
            bond_ac = bonds.harmClass(root,'%s %s'%(type[1],type[2]))
            ac  = bond_ac.r0
            ac2 = bond_ac.r0 * bond_ac.r0
        except:
            helpers.error('<urey="%s"> without existing bonds!'%self.type)

        cos_theta = (self.r0*self.r0 - ab2 - ac2)/(-2.0*ab*ac)
        sin_theta2 = 1 - cos_theta*cos_theta

        harm_k = ab2 * ac2 / self.r0 / self.r0 * sin_theta2 * self.k
        harm_theta0 = numpy.arccos(cos_theta)*180/3.14159265

        return harm_k, harm_theta0
            
    def check(self):
        import helpers 
        import molecule 
        existing = False
        for mol in self.root.findall("molecule"):
            current_molecule = molecule.moleculeElement(mol,mol.get('name'))
            for i in current_molecule.list("angle"):
                name = i.split()
                atom_i = molecule.atomClass(mol,name[0])
                atom_j = molecule.atomClass(mol,name[1])
                atom_k = molecule.atomClass(mol,name[2])
                if self.type == sequence(' '.join([atom_i.type,atom_j.type,atom_k.type])):
                    existing = True
        if not existing:
            helpers.warning('<urey type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return
        if self.k is None:
            helpers.error('<urey type="%s"> has no force constant!'%self.type)
        if self.r0 is None:
            helpers.error('<urey type="%s"> has no equilibrium distance!'%self.type)

            




