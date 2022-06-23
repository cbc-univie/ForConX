from xml.etree import ElementTree as ET

#===============================================================================================
# morseClass
#===============================================================================================
class morseClass:
    """
    manages dissociation energy, stretching constant and equilibrium distance of 
    a Morse potential V = D0 ( 1- exp[ -beta*(r-r0) ])^2.

    ______________________________________________________________________________________

    +  D0            dissociation energy
    +  r0            equilibrium distance
    +  beta          stretching constant
    +  root          memory address of root
    +  type          atom types of connected atoms    
    ______________________________________________________________________________________

    - __init__(root,pointer,type)  initializes Morse object 
    - __setattr__(attribute,value) overloads "=" operation
    ______________________________________________________________________________________
    """
    def __init__(self,root,pointer,type):
        self.root    = root
        self.type    = type
        self.pointer = root.find('./bonds/mors/[@type="%s"]'%type)
        if self.pointer is None:
            self.pointer = ET.SubElement(root.find('bonds'),'mors')
            self.pointer.set('type',type)
            self.D0 = None
            self.r0 = None
            self.beta = None
        else:
            self.D0 = self.get('D0')
            self.r0 = self.get('r0')
            self.beta = self.get('beta')

    def __setattr__(self,attribute,value):
        if attribute in 'root pointer':
            self.__dict__[attribute] = value
        elif attribute in 'type':
            value = (value.upper()).strip()
            self.__dict__[attribute] = value
            self.pointer.set(attribute,value)
        elif attribute in 'D0 r0 beta':
            try:
                self.__dict__[attribute] = float(value)
                self.pointer.set(attribute,str(round(value,5)))
            except TypeError:
                self.__dict__[attribute] = None

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

    def equal(self, D0, beta, r0):
        result = True
        if self.D0 != D0 :
            result = False
        if self.beta != beta :
            result = False
        if self.r0 != r0:
            result = False
        return result
