import xml.etree.ElementTree as ET
#======================================================================================
# input_output.mdElement
#======================================================================================
class mdElement:
    """
    This object handles <input> or <output>:
      <input md="..." >
        < ... />
        <energy         unit="..." />
        <distance       unit="..." />
        <polarizability unit="..." />
    </input>

    <output md="..." >
        < ... />
        <energy         unit="..." />
        <distance       unit="..." />
        <polarizability unit="..." />
    </output>
    ______________________________________________________________________________________
    + md           Name of MD program
    + energy       Pointer to <energy>
    + energy_unit  Energy unit
    + distance     Pointer to <distance>
    + distance_unit Distance unit
    + polarizability Pointer to <polarizability>
    + polarizability_unit Polarizability unit
    ______________________________________________________________________________________
    - __init__(root,inout)                 initializes <input> or <output> object
    - __setattr__(attribute,value)         overloads "=" operation
    + convert_energy(input,output)         returns corresponding conversion factor
    + convert_distance(input,output)       returns corresponding conversion factor
    + convert_polarizability(input,output) returns corresponding conversion factor
    + read_md()                            read force field of MD program
    + write_md()                           write force field of MD program
    ______________________________________________________________________________________
    """
    def __init__(self,root,inout):
        self.root = root
        if inout not in 'input output':
            import helpers
            helpers.error('inout tag is unknown!')
        self.inout = root.find(inout)
        if self.inout is None:
            self.inout = ET.SubElement(root,inout)
        self.md = self.inout.get('md')
        if self.md is None:
            self.md = 'XML'
            
        self.energy = self.inout.find('energy')
        if self.energy is None:
            self.energy = ET.SubElement(self.inout,'energy')
        self.energy_unit = self.energy.get('unit')

        self.distance = self.inout.find('distance')
        if self.distance is None:
            self.distance = ET.SubElement(self.inout,'distance')
        self.distance_unit = self.distance.get('unit')

        self.polarizability = self.inout.find('polarizability')
        if self.polarizability is None:
            self.polarizability = ET.SubElement(self.inout,'polarizability')
        self.polarizability_unit = self.polarizability.get('unit')

    def __setattr__(self,attribute,value):
        if attribute in 'root inout energy distance polarizability':
            self.__dict__[attribute] = value
        elif attribute in 'md':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.inout.set('md',value)
            except AttributeError:
                self.__dict__[attribute] = None
        elif attribute in 'energy_unit':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.energy.set('unit',value)
            except AttributeError:
                self.__dict__[attribute] = None
        elif attribute in 'distance_unit':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.distance.set('unit',value)
            except AttributeError:
                self.__dict__[attribute] = None
        elif attribute in 'polarizability_unit':
            try:
                value = (value.upper()).strip()
                self.__dict__[attribute] = value
                self.polarizability.set('unit',value)
            except AttributeError:
                self.__dict__[attribute] = None
                
    def convert_energy(self,input=None,output=None):
        """
        computes the factor arising from different energy units in 
        input and output file.
        """
        factor = 1.0
        if input is None:
            pointer = self.root.find('./input/energy')
            input = pointer.get('unit')
        if input[0:2]=="KJ":
            factor *=0.239
        elif input[0:2]=="EV":
            factor *= 23.0609
        elif input=="K":
            factor *= 0.001987912

        if output is None:
            pointer = self.root.find("./output/energy")
            output  = pointer.get('unit')
        if output[0:2]=="KJ":
            factor *= 4.184
        elif output[0:2]=="EV":
            factor *= 0.043363
        elif output=="K":
            factor *= 503.2227
        return factor

    def convert_distance(self,input=None,output=None):
        """
        computes the factor arising from different distance units in 
        input and output file.
        """
        factor = 1.0

        if input is None:
            pointer  = self.root.find('input/distance')
            input = pointer.get('unit')
        if input=="NM":
            factor *= 10

        if output is None:
            pointer = self.root.find('output/distance')
            output = pointer.get('unit')
        if output=="NM":
            factor *= 0.1
        return factor

    def convert_polarizability(self,input=None,output=None):
        """
        polarizability_unit computes the factor arising from different polarizability units in the 
        input and output file. The default energy unit is ANGSTROEM^3.
        """
        factor = 1.0

        if input is None:
            pointer = self.root.find('input/polarizability')
            input = pointer.get('unit')
        if input == "BOHR":
            factor *= 6.748342563
	if input == "NM":
	    factor *= 1000

        if output is None:
            pointer = self.root.find('output/polarizability')
            output = pointer.get('unit')
        if output == "BOHR":
            factor *= 0.148184534
	if output == "NM":
	    factor *= 0.001
        return factor

    def read_md(self):
        print '> Input force field from %s ...'%self.md
        current_input = mdElement(self.root,'input')
        
        # AMBER
        if self.md == 'AMBER':
            from ..amber import md2xml
            current_input.energy_unit = "KCAL"
            current_input.distance_unit = "ANGSTROEM"
            current_input.polarizability_unit = "ANGSTROEM"
            md2xml.read(self.root)
# TODO:            pointer_xyz = self.root.find('./input/rst')
#            if pointer_xyz is not None:
#                xyz_file = pointer_xyz.get('file')
#                xyz.rst2pdb(self.root,xyz_file)

        # CHARMM
        elif self.md == 'CHARMM':
            from ..charmm import md2xml
            from ..charmm import xyz
            current_input.energy_unit = "KCAL"
            current_input.distance_unit = "ANGSTROEM"
            current_input.polarizability_unit = "ANGSTROEM"
            
            md2xml.read(self.root)
            pointer_xyz = self.root.find('./input/crd')
            if pointer_xyz is not None:
                xyz_file = pointer_xyz.get('file')
                xyz.crd2pdb(self.root,xyz_file)

        # DLPOLY
        elif self.md == 'DLPOLY':
            from ..dlpoly import md2xml
            from ..dlpoly import xyz
            # current_input.energy_unit is set in the FIELD file
            current_input.distance_unit = "ANGSTROEM"
            current_input.polarizability_unit = "ANGSTROEM"
            
            md2xml.read(self.root)
            pointer_xyz = self.root.find('./input/config')
            if pointer_xyz is not None:
                xyz_file = pointer_xyz.get('file')
                xyz.config2pdb(self.root,xyz_file)

        # GROMACS
        elif self.md == 'GROMACS':
            from ..gromacs import md2xml
            current_input.energy_unit = "KJ"
            current_input.distance_unit = "NM"
	    current_input.polarizability_unit="NM"

            topfile=self.root.find('./input/top').get('file')
            lines=md2xml.preprocess(topfile)
            md2xml.read_top(self.root,lines)

        # LAMMPS
        elif self.md == 'LAMMPS':
            from ..lammps import md2xml
            md2xml.read(self.root)

        # XML
        elif self.md == 'XML':
            pass

        # UNKNOWN
        else:
            import helpers 
            helpers.error('MD program is unknown.')
        
    def write_md(self):
        print '> Output force field to %s ...'%self.md
        current_output = mdElement(self.root,'output')

        # AMBER
        if self.md=='AMBER':
            from ..amber import xml2md
            current_output.energy_unit = "KCAL"
            current_output.distance_unit = "ANGSTROEM"
            current_output.polarizability_unit = "ANGSTROEM"
            xml2md.write_lib(self.root)
            xml2md.write_frcmod(self.root)

        # CHARMM
        elif self.md=='CHARMM':
            from ..charmm import xml2md
            from ..charmm import xyz

            current_output.energy_unit = "KCAL"
            current_output.distance_unit = "ANGSTROEM"
            current_output.polarizability_unit = "ANGSTROEM"
            xml2md.write_topology(self.root)
            xml2md.write_parameter(self.root)
            pointer_xyz = self.root.find('./output/crd')
            if pointer_xyz is not None:
                xyz_file = pointer_xyz.get('file')
                xyz.pdb2crd(self.root,xyz_file)

        # DLPOLY
        elif self.md=='DLPOLY':
            from ..dlpoly import xml2md
            from ..dlpoly import xyz
            import nonbonded
            import helpers

            if current_output.energy_unit is None:
                current_output.energy_unit = "KCAL"
            current_output.distance_unit = "ANGSTROEM"
            current_output.polarizability_unit = "ANGSTROEM"
            
            print '\t Converting atom-based van-der-Waals XML information to XML vdw pairs ...'
            current_nonbonded = nonbonded.nonbondedElement(self.root)
            current_nonbonded.atom2vdw()
            helpers.write_xml(self.root,'forconx.xml')

            xml2md.write_field(self.root)

            pointer_xyz = self.root.find('./output/config')
            if pointer_xyz is not None:
                xyz_file = pointer_xyz.get('file')
                xyz.pdb2config(self.root, xyz_file)
            
        # GROMACS
        elif self.md=='GROMACS':
            from ..gromacs import xml2md
            current_output.energy_unit = "KJ"
            current_output.distance_unit = "NM"
            current_output.polarizability_unit = "NM"  
            xml2md.write(self.root)

        # LAMMPS
        elif self.md=='LAMMPS':
            from ..lammps import xml2md
            xml2md.write(self.root)

        # XML
        elif self.md=='XML':
            pass

        # UNKNOWN
        else:
            helpers.error('MD program is unfortunately unknown.')

#===============================================================================================
# input_output.header
#===============================================================================================
def header(f,comment):
    """
    comment is the symbol to indicate a comment in a particular md force field parameter file. 
    This comment symbol varies within the MD packages.
    """
    f.write(comment+"."*70+"\n")
    f.write(comment+" ForConX version 0.1\n")
    f.write(comment+" \n")
    f.write(comment+" Authors in alphabetical order:\n")
    f.write(comment+" Alain Dequidt, Diddo Diddens, Benjamin Golub, Volker Lesch, \n")
    f.write(comment+" Christian Schroeder, Marcello Sega, Veronika Zeindlhofer\n")
    f.write(comment+"."*70+"\n")
    return

#===============================================================================================
# input_output.conk
#===============================================================================================
def conk(f,comment):
    """
    comment is the symbol to indicate a comment in a particular md force field parameter file. 
    This comment symbol varies within the MD packages.
    """
    f.write(comment+"-::://////+ossoos+++///::::::+oo`                  F O R c e   f i e l d    C O N v e r s i o n\n")                  
    f.write(comment+"-:::://:::/+shsssoo++/::----:::/y-                           o n   X m l   b a s i s \n")                
    f.write(comment+":::::::::://ydyyso+o+/::-:----::/s+.\n")              
    f.write(comment+":::::-.::///dhyysoooo++/::::--::::/y/              Please cite:        \n")             
    f.write(comment+"--:--..:::/ydhysssossso/////:///::::so.            `A force field Conversion tool based on XML`\n")           
    f.write(comment+".--:-`-:::ohhyyoooossssshyoo+++////:/+ys.          \n")         
    f.write(comment+"-.-:-`-//sddhsossssysssyyy+:://:--:::://o+:        V. Lesch, D. Diddens,\n")       
    f.write(comment+"--::--:/+yhhyssssyssooossyo://::--:::://+/o+       C. E. S. Bernardes, B. Golub, J. N. Canongia Lopes,\n")      
    f.write(comment+"-:/:-:://+hyssoshhsooooosso/o//////////////++      A. Dequidt,\n")     
    f.write(comment+"-:::--:.-osooyyshhyyyyysso+//++++++/++//////+y`    V. Zeindlhofer, M. Sega and C. Schroeder\n")   
    f.write(comment+"-.-:-:/+osysossyyddmdhyso+++++/+/::::///:/://oy    J. Comput. Chem. 38 (2017), 629\n")   
    f.write(comment+"...-:+yhhysoo+++sdddddhsoooo+/:/::-:--::/:////o/\n")  
    f.write(comment+"`../hdhyyys+++////sdddhy+ss+//:++-/.:-.-:-/+///o\n")  
    f.write(comment+"..+ddddhyyy/s+/://:+sddhsss+/+/+///:---:::/+///o\n")  
    f.write(comment+"-/hddddhhhs+s+++/--//sdhsyo+/-:/+/://::::///::/s\n")  
    f.write(comment+":omdddddyhhsy+so/://:+shhyyo+++///////:/:::::+s.\n")  
    f.write(comment+"-:ymmddddyddyyssoo/++//yhhysooso++++////::::/y:\n")   
    f.write(comment+"..:ohdddmmmmmmddhyssysoydhhyysssooo+///:::::o.\n")    
    f.write(comment+"----::/+osyhhdddmNNNNNNNmdhyyyyssoooo++///+o.\n")     
    f.write(comment+":/:::/::::-----/+ydmNNNNNNNmddhsssososyo/:`\n")
    return

#===============================================================================================
# input_output.warranty
#===============================================================================================
def warranty(f,comment):
    """
    comment is the symbol to indicate a comment in a particular md force field parameter file. 
    This comment symbol varies within the MD packages.
    """
    f.write(comment+" WARNING: FORCE FIELD CONVERSION MAY NOT BE COMPLETE \n")
    f.write(comment+" NOR FULLY CORRECT. PLEASE CHECK THESE AUTOMATICALLY \n")
    f.write(comment+" GENERATED FILES CAREFULLY BEFORE USING THEM.\n")
    f.write(comment+" \n")
    f.write(comment+" THE PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL \n")
    f.write(comment+" BE USEFUL, BUT WITHOUT ANY WARRANTY.\n")
    f.write(comment+" IT IS PROVIDED 'AS IS' WITHOUT WARRANTY OF ANY KIND,\n")
    f.write(comment+" EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n")
    f.write(comment+" THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A \n")
    f.write(comment+" PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND \n")
    f.write(comment+" PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE \n")
    f.write(comment+" PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY\n")
    f.write(comment+" SERVICING, REPAIR OR CORRECTION.\n")
    f.write(comment+" \n")
    f.write(comment+" IN NO EVENT THE AUTHORS WILL BE LIABLE TO YOU FOR DAMAGES, \n") 
    f.write(comment+" INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES\n") 
    f.write(comment+" (INCLUDING BUT NOT LIMITED ARISING OUT OF THE USE OR \n")
    f.write(comment+" INABILITY TO USE THE PROGRAM TO LOSS OF DATA OR DATA BEING \n")
    f.write(comment+" RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES\n")
    f.write(comment+" OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),\n")
    f.write(comment+" EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.\n")
    f.write(comment+" \n")
    return
