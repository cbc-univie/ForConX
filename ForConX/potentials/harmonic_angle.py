from xml.etree import ElementTree as ET
#=============================================================================================
# hangleClass
#=============================================================================================
class hangleClass:
    """
    The harmonic angle potential V = k ( theta-theta0 )^2 is defined by the 
    force constant k and the equilibrium angle theta0.

    ______________________________________________________________________________________
    +  k                    force constant
    +  theta0               equilibrium angle
    +  type                 atom types of the involved atoms
    +  element              "angle" or "improper"
    +  root                 memory address of root
    +  pointer              memory address of current harmonic angle
    ______________________________________________________________________________________
    - __init__(root,element,type)  initializes harmonic angle object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes harmonic angle object
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
            self.k = None
            self.theta0 = None
        else:
            self.k = self.get('k')
            self.theta0 = self.get('theta0')
           
    def __setattr__(self,attribute,value):
        if attribute in 'root pointer element':
            self.__dict__[attribute] = value
        elif attribute in 'type':
            value = (value.upper()).strip()
            self.__dict__[attribute] = value
        elif attribute in 'k':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,6)))
        elif attribute in 'theta0':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,6)))

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

    def equal(self,k,theta0):
        result = True
        if self.k != k or self.theta0 != theta0:
            result = False
        return result
