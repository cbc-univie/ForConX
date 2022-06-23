from xml.etree import ElementTree as ET
#===============================================================================================
# ryckaertClass
#===============================================================================================
class ryckaertClass:
    """
    manages force constants k of a Ryckaert-Bellemans potential:
    V = sum_i k_i cos(psi)^i with phi = 180 deg-psi.

    ______________________________________________________________________________________
    + k       force constants
    + type    atom types of involved atoms
    + root    memory address of root
    + pointer memory address of current Ryckaert Bellemans potential
    ______________________________________________________________________________________
    - __init__(root,element,type)  initializes ryckaert object
    - __setattr__(attribute,value) overloads "=" operation
    + remove()                     removes ryckaert object
    ______________________________________________________________________________________
    """
    def __init__(self,root,element,type):
        self.root = root
        self.type = type
        tmp = element.split('/')
        self.element = self.root.find(tmp[1])
        self.pointer = root.find('%s/[@type="%s"]'%(element,type))
        if self.pointer is None:
            self.pointer = ET.SubElement(self.element,'ryck')
            self.pointer.set('type',type)
            self.k = []
        else:
            self.k = map(float,(self.pointer.get('k')).split())
            
    def __setattr__(self,attribute,value):
        if attribute in 'k':
            self.__dict__[attribute] = value
            try:
                string = ''
                for i in value:
                    string += str(round(i,4)) + ' '
                self.pointer.set(attribute,string.strip())
            except TypeError:
                from ..md_xml import helpers
                helpers.error('<%s/ryck> k must be an array'%self.element.tag)
        elif attribute in 'root pointer element':
            self.__dict__[attribute] = value
        elif attribute in 'type':
            value = (value.upper()).strip()
            self.__dict__[attribute] = value
                
    def remove(self):
        self.pointer.clear()
        self.dihedrals.remove(self.pointer)
        self.pointer = None
        del self
    
    # Conversion to cosine potentials        
    def cos(self):
        from ..md_xml import helpers 
        kryck = self.k
        try:
            kcos = [-kryck[1]-0.75*kryck[3]-0.625*kryck[5],
                    -0.5*(kryck[2]+kryck[4]),
                    -0.25*kryck[3]-0.3125*kryck[5],
                    -0.125*kryck[4],
                    -0.0625*kryck[5]]
            n     = [   1,   2,   3,   4,   5 ]
            delta = [ 180.0, 180.0, 180.0, 180.0, 180.0 ]
            shift =  kryck[0]-kryck[1]+kryck[2]-kryck[3]+kryck[4]-kryck[5]
            if abs(shift)>0.01:
                helpers.warning('<ryck> conversion has a constant energy shift %8.3f!'%shift)
            return n,kcos,delta
        except:
            return None, None, None

    # Conversion from Fourier series (see Gromacs manual bonded interactions)
    def fourier(self,k):
        rb = []
        rb.append( k[1] + 0.5 * ( k[0] + k[2] ) )
        rb.append( 0.5 * ( -k[0] + 3 * k[2] ) )
        rb.append( - k[1] + 4 * k[3] )
        rb.append( -2 * k[2] )
        rb.append( -4 * k[3] )
        rb.append( 0.0 )
        self.k = rb
