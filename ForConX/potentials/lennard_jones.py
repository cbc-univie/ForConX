from xml.etree import ElementTree as ET
#=============================================================================================
# ljClass
#=============================================================================================
class ljClass:
    """
    manages epsilon and sigma of the Lennard-Jones potential

    V = 4 epsilon ( (sigma/r)^12 - (sigma/r)^6 )

    ______________________________________________________________________________________
    + sigma      LJ - sigma
    + epsilon    LJ - epsilon
    + type       atom types of involved atoms
    + element    "atom" or "vdw" 
    + root       memory address of root
    + pointer    memory address of current Lennard-Jones
    ______________________________________________________________________________________
    - __init__(root,element,type)  initializes LJ object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes LJ object
    ______________________________________________________________________________________
    """
    def __init__(self,root,element,type):
        self.root = root
        self.type = type
        self.nonbonded = self.root.find('./nonbonded')
        if self.nonbonded is None:
            self.nonbonded = ET.SubElement(root,'nonbonded')
        tmp = element.split('/')
        self.element = self.root.find(tmp[1])
        self.pointer = root.find('%s/[@type="%s"]'%(element,type))
        if self.pointer is None:
            self.pointer = ET.SubElement(self.element,tmp[2])
            self.pointer.set('type',type)
            self.sigma   = None
            self.epsilon = None
            self.elec14  = None
            self.vdw14   = None
        else:
            self.sigma   = self.get('sigma')
            self.epsilon = self.get('epsilon')
            self.elec14  = self.get('elec14')
            self.vdw14   = self.get('vdw14')
                
    def __setattr__(self,attribute,value):
        if attribute in 'root pointer element nonbonded':
            self.__dict__[attribute] = value
        elif attribute in 'type':
            value = (value.upper()).strip()
            self.__dict__[attribute] = value
        elif attribute in 'sigma epsilon elec14 vdw14':
            self.__dict__[attribute] = value
            if value is not None:
                self.pointer.set(attribute,str(round(value,5)))

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

    
    
