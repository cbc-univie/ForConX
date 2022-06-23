import xml.etree.ElementTree as ET

#===============================================================================================
# Currently supported potentials:
# -cos   V = k ( 1 + cos[n*phi - delta] )
# -ryck  V = sum_{i=0}^5 k_i (cos phi)^i
#===============================================================================================
from ..potentials import cosine_torsion     
from ..potentials import ryckaert_bellemans 

#===============================================================================================
# dihedrals.sequence
#===============================================================================================
def sequence(type):
    """
    sorts alphabetically the atomtypes in dihedrals to facilitate the search for duplicates, i.e. 
    {A-B-C-D, D-C-B-A} --> A-B-C-D
    """
    tmp = type.split()
    if tmp[0]>tmp[3]:
        tmp = tmp[::-1]
    elif tmp[0]==tmp[3] and tmp[1]>tmp[2]:
        tmp = tmp[::-1]
    return ' '.join(tmp)

#===============================================================================================
# dihedralsElement                                                                      
#===============================================================================================
class dihedralsElement:
    """
    An XML entry handled by this class should look like this:
    <dihedrals>
       <cos delta="... ..." k="... ..." n="... ..." type="..." />
             ...
       <cos delta="..." k="..." n="..."  type="..." />
       <ryck k="... ... ... ... ... ..." type="..." />
             ...
       <ryck k="... ... ... ... ... ..." type="..." />
    </dihedrals>
    
                                << interface >>
    ______________________________________________________________________________________
    + root                     memory address of root
    + pointer                  memory address of <dihedrals>
    ______________________________________________________________________________________
    - __init__(root)           initializes <dihedrals> object
    + remove()                 removes <dihedrals> object
    + find_type(type)          returns memory address of the dihedral with the 
                               atom types "... ..."
    + list(tag)                returns a list of all types of a particular 
                               subelement of <dihedrals>
    + list_molecule_dihedral() returns a list of all dihedrals of <molecule>
    + check()                  checks the XML sanity of the <dihedrals> object
    ______________________________________________________________________________________
    """
    def __init__(self,root):
        self.root = root
        self.pointer = root.find('dihedrals')
        if self.pointer is None:
            self.pointer = ET.SubElement(root,'dihedrals')
        
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
                list.append(i.get("type"))
        return list

    def list_molecule_dihedral(self):
        import molecule
        list = []
        for mol in self.root.findall("molecule"):
            current_molecule = molecule.moleculeElement(mol,mol.get('name'))
            for i in current_molecule.list('dihedral'):
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
            for i in mol.findall('dihedral'):
                current_dihedral = molecule.dihedralClass(mol,i.get('name'))
                existing = False
                number += 1
                if self.find_type(current_dihedral.type) is not None:
                    existing = True
                if not existing:
                    helpers.error('<dihedral name="%s"> in molecule %s has no corresponding dihedral potential'%(i.get("name"),mol.get("name")))
        print "%5s potentials checked."%number
        
        print "\t<cos>  ..."
        number = 0
        for i in self.pointer.findall('cos'):
            number += 1
            type = i.get("type")
            cos = cosClass(self.root,type)
            cos.check()
        print "\n\t\t\t%5s potentials checked." %number
        
        print "\t<ryck> ..."
        number = 0
        for i in self.pointer.findall('ryck'):
            number += 1
            type = i.get("type")
            ryck = ryckClass(self.root,type)
            ryck.check()
        print "\n\t\t\t%5s potentials checked." %number

        for i in self.root.findall('./dihedrals/*'):
            if i.tag not in 'cos ryck':
                print "?\t<%s>\n"%i.tag

#===============================================================================================
# dihedrals.cosClass                                                                          
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
    + ryck()                       converts cosine torsion to Ryckaert-Bellemans
    + check()                      checks the XML sanity of the cosClass object
    ______________________________________________________________________________________

    """
    def __init__(self,root,type):
        if root.find('dihedrals') is None:
            ET.SubElement(root,'dihedrals')
        cosine_torsion.cosineClass.__init__(self,root,'./dihedrals/cos',sequence(type))
        
    def check(self):
        import helpers
        if self.type is None:
            error('<cos> has no atomtypes!')

        # Testing existence of parameters and their completeness
        try:
            k = self.pointer.get("k").split()
            n = self.pointer.get("n").split()
            delta = self.pointer.get("delta").split()
            len_k     = len(k)
            len_n     = len(n)
            len_delta = len(delta)
        except:
            len_k     = 0
        if not len_k>0:
            helpers.error('<cos type="%s"> has no force constant k'%self.type)
        if not len_n:
            helpers.error('<cos type="%s"> has no multiplicity n'%self.type)
        if not len_delta:
            helpers.error('<cos type="%s"> has no phase shift delta'%self.type)
        if (len_k != len_n) or (len_k != len_delta):
            helpers.error('<cos type="%s"> Cosine potential is not complete!' %self.type)

        # Testing if this potential is necessary for the molecules under investigation            
        current_dihedrals = dihedralsElement(self.root)
        if not self.type in current_dihedrals.list_molecule_dihedral():
            warning('<cos type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return

        # sorting dihedrals
        tmp_k = {}
        tmp_delta = {}
        for i in range(len_n):
            tmp_k[n[i]] = k[i]
            tmp_delta[n[i]] = delta[i]

        new_k = []
        new_n = []
        new_delta = []
        for i in sorted(set(n)):
            new_k.append(float(tmp_k[i]))
            new_n.append(int(i))
            new_delta.append(float(tmp_delta[i]))
        self.k = new_k
        self.n = new_n
        self.delta = new_delta
            
    # Conversion to Ryckaert-Bellemans        
    def ryck(self):
        import helpers
        kcos = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        for i in range(len(self.n)):
            if (self.n[i]<1) or (self.n[i]>5): 
                return None
            if not (self.delta[i]==0.):
                return None
            kcos[self.n[i]] = self.k[i]
        kryck = [-kcos[1] + 2*kcos[2] - kcos[3] - kcos[5],
                 kcos[1] - 3*kcos[3]+ 5*kcos[5],
                 -2*kcos[2] + 8*kcos[4],
                 4*kcos[3] - 20*kcos[5],
                 -8*kcos[4],
                 16*kcos[5]]
        shift = kryck[0]+0.5*kryck[2]+0.375*kryck[4]
        helpers.warning('Ryckaert-Bellemans conversion has a constant energy shift %8.3f!'%shift)
        return kryck

    
#===============================================================================================
# dihedrals.ryckClass                                                                         
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
        if root.find('dihedrals') is None:
            ET.SubElement(root,'dihedrals')
        ryckaert_bellemans.ryckaertClass.__init__(self,root,'./dihedrals/ryck',sequence(type))

    def check(self):
        import helpers
        if self.type is None:
            helpers.error('<ryck> has no type')
        if self.k is None:
            helpers.error('<ryck type="%s"> has no force constants!'%self.type)
        current_dihedrals = dihedralsElement(self.root)
        if not self.type in current_dihedrals.list_molecule_dihedral():
            helpers.warning('<ryck type="%s"> is not part of a molecule and therefore removed.' %self.type)
            self.remove()
            return

