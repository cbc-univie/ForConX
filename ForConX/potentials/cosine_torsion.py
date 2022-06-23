from xml.etree import ElementTree as ET

#===============================================================================================
# cosineClass
#===============================================================================================
class cosineClass:
    """
    The cosine torsion potential V = sum_i k_i [ 1 - cos(n_i phi -delta_i) ] may contain
    contributions of various multiplicities n concerning the very same torsion.

    ______________________________________________________________________________________
    +  k       force constant
    +  n       multiplicity
    +  delta   phase shift
    +  type    atom types of involved atoms
    +  root    memory address of root
    +  pointer memory address of this cosine potential
    +  element "dihedral" or "improper"
    _______________________________________________________________________________________
    - __init__(root,element,type)  initializes cosine object
    - __setattr__(attribute,value) overloads "=" operation
    + sort()
    ______________________________________________________________________________________
    """
    def __init__(self,root,element,type):
        self.root = root
        self.type = type
        tmp = element.split('/')
        self.element = self.root.find(tmp[1])
        self.pointer = root.find('%s/[@type="%s"]'%(element,type))
        if self.pointer is None:
            self.pointer = ET.SubElement(self.element,tmp[2])
            self.pointer.set('type',type)
            self.k = []
            self.n = []
            self.delta = []
        else:
            self.k = map(float,(self.pointer.get('k')).split())
            self.n = map(int,(self.pointer.get('n')).split())
            self.delta = map(float,(self.pointer.get('delta')).split())
            
    def __setattr__(self,attribute,value):
        if attribute =='n':
            self.__dict__[attribute] = value
            try:
                string = ''
                for i in value:
                    string += str(int(i)) + ' '
                self.pointer.set(attribute,string.strip())
            except TypeError:
                from ..md_xml import helpers
                helpers.error('<%s/cos> n must be an array.'%self.element.tag)
                
        elif attribute in 'k delta':
            self.__dict__[attribute] = value
            try:
                string = ''
                for i in value:
                    string += str(round(i,4)) + ' '
                self.pointer.set(attribute,string.strip())
            except TypeError:
                from ..md_xml import helpers
                helpers.error('<%s/cos> k and delta must be arrays.'%self.element.tag)
        elif attribute in 'root pointer element':
            self.__dict__[attribute] = value
        elif attribute in 'type':
            value = (value.upper()).strip()
            self.__dict__[attribute] = value
            
    def remove(self):
        self.pointer.clear()
        self.element.remove(self.pointer)
        self.pointer = None
        del self

