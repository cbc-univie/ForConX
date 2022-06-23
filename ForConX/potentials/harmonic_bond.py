from xml.etree import ElementTree as ET
#===============================================================================================
# hbondClass
#===============================================================================================
class hbondClass:
    """
    The harmonic bond potential V = k ( r-r0 )^2 is defined by the 
    force constant k and the equilibrium distance r0.

    ______________________________________________________________________________________
    + k             force constant
    + r0            equilibrium distance
    + type          atom types of involved atoms
    + element       "bond" or "urey"
    + root          memory address of root
    + pointer       memory address of current harmonic bond
    ______________________________________________________________________________________
    - __init__(root,element,type)  initializes harmonic bond object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes harmonic bond object
    ______________________________________________________________________________________
    """
    def __init__(self,root,element,type):
        self.root    = root
        self.type    = type
        tmp = element.split('/')
        self.element = self.root.find(tmp[1])
        self.pointer = root.find('%s/[@type="%s"]'%(element,type))
        if self.pointer is None:
            self.pointer = ET.SubElement(self.element,tmp[2])
            self.pointer.set('type',type)
            self.k  = None
            self.r0 = None
        else:
            self.k  = self.get('k')
            self.r0 = self.get('r0')
                
    def __setattr__(self,attribute,value):
        if attribute in 'root pointer element':
            self.__dict__[attribute] = value
        elif attribute in 'r0 k':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,6)))
        elif attribute in 'type':
            value = (value.upper()).strip()
            self.__dict__[attribute] = value

    def get(self,attribute):
        try:
            value = float(self.pointer.get(attribute))
        except:
            value = None
        return value
            
    def remove(self):
        self.pointer.clear()
        self.element.remove(self.pointer)
        self.pointer = None
        del self

    def equal(self, k, r0):
        result = True
        if self.k != k or self.r0 != r0:
            result = False
        return result
    
            

