def main():
    """
    ForConX main program
    """
    import os
    import sys
    import xml.etree.ElementTree as ET
    from md_xml import helpers     
    from md_xml import input_output 
    
    print '='*120
    print '                      FORce field CONversion on Xml basis'
    print '='*120
    print 'Authors in alphabetical order:'
    print 'Alain Dequidt, Diddo Diddens, Benjamin Golub, Volker Lesch, Christian Schroeder, Marcello Sega, Veronika Zeindlhofer\n'

    #.........................................................................................
    # Stage 1 - Reading the force field file or the XML structure 
    #.........................................................................................
    print '-'*120
    print '1. Reading MD force field files and convert to XML tree'
    print '-'*120
    # checking if input file is given in command line
    if len(sys.argv) < 2:
        print 'usage: forconx input.xml'
        sys.exit()
        
    # try to open input file and read xml structure    
    xmlfile = sys.argv[1]
    try:
        root =  ET.parse(xmlfile).getroot()
    except:
        helpers.error('XML file %s is not readable'%xmlfile)
        
    # Reading input section
    current_input = input_output.mdElement(root,'input')
    current_input.read_md()
    helpers.write_xml(root,'forconx_1.xml')

    #.........................................................................................
    # Stage 2 - Checking the sanity of the XML structure
    #.........................................................................................
    print "-"*120
    print "2. Checking XML tree"
    print "-"*120
    helpers.check_xml(root)
    helpers.write_xml(root,'forconx_2.xml')
    
    #.........................................................................................
    # Stage 3 - Based on the XML structure the new force field files are written 
    #.........................................................................................
    print "-"*120
    print "3. Writing force field files based on current XML structure"
    print "-"*120
    current_output = input_output.mdElement(root,'output')
    current_output.write_md()
    print
    helpers.write_xml(root,'forconx_3.xml')

