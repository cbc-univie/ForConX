import xml.etree.ElementTree as ET

#===============================================================================================
# Currently supported potentials:
# - cos      V = k ( 1 + cos[n*phi - delta] )
# -harm      V = k ( theta - theta0)^2
# -ryck      V = sum_{i=0}^5 k_i (cos phi)^i
#===============================================================================================
from ..potentials import cosine_torsion 
from ..potentials import harmonic_angle 
from ..potentials import ryckaert_bellemans

#===============================================================================================
# impropers.sequence
# MD program      central atom
# AMBER               3.
# CHARMM              1. or 4.
# DLPOLY              2.
# GROMACS             1. or 4.
#===============================================================================================
def sequence(type):
    """
    Currently disabled !!!

    sorts the improper atoms alphabetically. 
    The first atom should always be the central atom.
    The last  atom should always be the atom moving out-of-plane.
    {A-B-C-D, A-C-B-D}  ->  A-B-C-D
    """
    tmp = type.split()
#    if tmp[1] > tmp[2]:
#        tmp1 = tmp[1]
#        tmp[1] = tmp[2]
#        tmp[2] = tmp1
    return ' '.join(tmp)


#===============================================================================================
# impropersElement                                                                       
#===============================================================================================
class impropersElement:
    """
    An XML entry handled by this class should look like this:
    <impropers>
       <cos delta="... ..." k="... ..." n="... ..." type="..." />
             ...
       <cos delta="..." k="..." n="..."  type="..." />
       <harm k="..." phi0="..." type="..." />
             ...
       <harm k="..." phi0="..." type="..." />
    </impropers>

                                << interface >>
    ______________________________________________________________________________________
    + root                     memory address of root
    + pointer                  memory address of <impropers>
    ______________________________________________________________________________________
    - __init__(root)           initializes <impropers> object
    + remove()                 removes <impropers> object
    + find_type(type)          returns memory address of the impropers with the 
                               atom types "... ..."
    + list(tag)                returns a list of all types of a particular 
                               subelement of <impropers>
    + list_molecule_improper() returns a list of all impropers of <molecule>
    + check()                  checks the XML sanity of the <impropers> object
    ______________________________________________________________________________________
    """
    def __init__(self,root):
        self.root = root
        self.pointer = root.find('impropers')
        if self.pointer is None:
            self.pointer = ET.SubElement(root,'impropers')
            
    def remove(self):
        self.pointer.clear()
        del self

    def find_type(self,type):
        pointer = self.pointer.find('./*/[@type="%s"]'%type)
        try:
            return pointer.tag.upper(),pointer
        except:
            return None, None
        
    def list(self,tag):
        list = []
        for i in self.pointer:
            if i.tag.upper() in tag.upper():
                list.append(i.get('type'))
        return list
    
    def list_molecule_improper(self):
        import molecule
        list = []
        for mol in self.root.findall("molecule"):
            current_molecule = molecule.moleculeElement(mol,mol.get('name'))
            for i in current_molecule.list('improper'):
                name = i.split()
                atom_i = molecule.atomClass(mol,name[0])
                atom_j = molecule.atomClass(mol,name[1])
                atom_k = molecule.atomClass(mol,name[2])
                atom_l = molecule.atomClass(mol,name[3])
                type = sequence(' '.join([atom_i.type,
                                          atom_j.type,
                                          atom_k.type,
                                          atom_l.type]))
                list.append(type)
        return list
    
    def check(self):
        import helpers
        import molecule
        print "\t<molecule> ...\t",
        number = 0
        for mol in self.root.findall('./molecule'):
            for i in mol.findall('improper'):
                current_improper = molecule.improperClass(mol,i.get('name'))
                existing = False
                number += 1
                if self.find_type(current_improper.type) is not None:
                    existing = True
                if not existing:
                    helpers.error('<improper name="%s"> in molecule %s has no corresponding improper potential'%(i.get("name"),mol.get("name")))
        print "%5s potentials checked."%number
        
        print "\t<cos>  ..."
        number = 0
        for i in self.pointer.findall('cos'):
            number += 1
            type = i.get("type")
            cos = cosClass(self.root,type)
            cos.check()
        print "\n\t\t\t%5s potentials checked."%number

        print "\t<ryck> ..."
        number = 0
        for i in self.pointer.findall('ryck'):
            number += 1
            type = i.get("type")
            ryck = ryckClass(self.root,type)
            ryck.check()
        print "\n\t\t\t%5s potentials checked." %number
        
        print "\t<harm> ..."
        number = 0
        for i in self.pointer.findall('harm'):
            number += 1
            type = i.get("type")
            harm = harmClass(self.root,type)
            harm.check()
        print "\n\t\t\t%5s potentials checked."%number

        for i in self.root.findall('./impropers/*'):
            if i.tag not in 'harm cos ryck':
                print "?\t<%s>\n"%i.tag

                
#===============================================================================================
# impropers.cosClass
#===============================================================================================
class cosClass(cosine_torsion.cosineClass):
    """
    An XML entry handled by this class should look like this:
       <cos delta="... ..." k="... ..." n="... ..." type="..." />
             ...
       <cos delta="..." k="..." n="..."  type="..." />

    ______________________________________________________________________________________
    + k         force constant
    + n         multiplicity
    + delta     phase shift
    + type      type of the involved atoms
    + root      memory address of root
    ______________________________________________________________________________________
    - __init__(root,type)          initializes cosClass object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes current cosmClass object (including XML entry)
    + check()                      checks the XML sanity of the cosClass object
    ______________________________________________________________________________________
    """
    def __init__(self,root,type):
        if root.find('impropers') is None:
            ET.SubElement(root,'impropers')
        cosine_torsion.cosineClass.__init__(self,root,'./impropers/cos',sequence(type))

    def check(self):
        import helpers
        if not self.pointer.get("type"):
            helpers.error('<cos type="%s"> has no atomtypes!'%self.type)

#       Testing existence of parameters and their completeness
        len_k     = len(self.pointer.get("k").split())
        len_n     = len(self.pointer.get("n").split())
        len_delta = len(self.pointer.get("delta").split())                                                        
        if not len_k>0:
            helpers.error('<cos type="%s"> has no force constant k'%self.type)
        if not len_n:
            helpers.error('<cos type="%s"> has no multiplicity n'%self.type)
        if not len_delta:
            helpers.error('<cos type="%s"> has no phase shift delta'%self.type)
        if (len_k != len_n) or (len_k != len_delta):
            helpers.error('<cos type="%s"> Cosine potential is not complete!' %self.type)

#       Testing if this potential is necessary for the molecules under investigation            
        current_impropers = impropersElement(self.root)
        if not self.type in current_impropers.list_molecule_improper():
            helpers.warning('<cos type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return

#===============================================================================================
# impropers.ryckClass                                                                         
#===============================================================================================
class ryckClass(ryckaert_bellemans.ryckaertClass):
    """
    An XML entry handled by this class should look like this:
       <ryck k="... ... ... ... ... ..." type="..." />
             ...
       <ryck k="... ... ... ... ... ..." type="..." />

    ______________________________________________________________________________________
    + k         force constant
    + type      type of the involved atoms
    + root      memory address of root
    ______________________________________________________________________________________
    - __init__(root,type)          initializes cosClass object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes current cosmClass object (including XML entry)
    + cos()                        converts Ryckaert-Bellemans to cosine torsion  
    + check()                      checks the XML sanity of the cosClass object
    ______________________________________________________________________________________
    """
    def __init__(self,root,type):
        if root.find('impropers') is None:
            ET.SubElement(root,'impropers')
        ryckaert_bellemans.ryckaertClass.__init__(self,root,'./impropers/ryck',sequence(type))

    def check(self):
        import helpers
        if self.type is None:
            helpers.error('<ryck> has no type')
        if self.k is None:
            helpers.error('<ryck type="%s"> has no force constants!'%self.type)
        current_impropers = impropersElement(self.root)
        if not self.type in current_impropers.list_molecule_improper():
            helpers.warning('<ryck type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return
        
#===============================================================================================
# impropers.harmClass                                                                         
#===============================================================================================
class harmClass(harmonic_angle.hangleClass):
    """
    harmClass is an object of a particular ./impropers/harm which should look like this:
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
        if root.find('impropers') is None:
            ET.SubElement(root,'impropers')
        harmonic_angle.hangleClass.__init__(self,root,'./impropers/harm',sequence(type))

    def check(self):
        import helpers
        if self.type is None:
            helpers.error('<cos> has no atomtypes!')
        if self.k is None:
            helpers.error('<cos type="%s"> has no force constants k!'%self.type)
        if self.theta0 is None:
            helpers.error('<cos type="%s"> has no equilibrium angle theta0!'%self.type)
        current_impropers = impropersElement(self.root)
        if not self.type in current_impropers.list_molecule_improper():
            helpers.warning('<cos type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return

        
