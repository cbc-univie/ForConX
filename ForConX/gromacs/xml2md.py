#!/usr/bin/python
import xml.etree.ElementTree as ET
import re,sys
import os.path
from operator import itemgetter, attrgetter

#==================================================================================================
# Gromacs output
#==================================================================================================
def write(root):
    from ..md_xml import input_output

    current_output = input_output.mdElement(root,'output')
    energy_conversion   = current_output.convert_energy()
    distance_conversion = current_output.convert_distance()
    polarizability_conversion = current_output.convert_polarizability()
    
    import parameter
    parameter.create_forcefield_directory('forconx.ff')
    parameter.write_forcefielditp(root)
    parameter.write_ffnonbonded(root,energy_conversion,distance_conversion)
    parameter.write_atomtypes(root)
    parameter.write_ffbonded(root,energy_conversion,distance_conversion)
    parameter.write_rtp(root,energy_conversion,distance_conversion)
    parameter.write_forcefielddoc(root)
    
                
