#!/usr/bin/env python3
# 
# Copyright (C) 2009 INRA
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie - Plateforme genomique Toulouse'
__copyright__ = 'Copyright (C) 2009 INRA'
__license__ = 'GNU General Public License'
__version__ = '2.4.0'
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'prod'

import argparse,re,sys


def cigar2statesRead (sam_line_obj):
    """
      Returns the status of each position on read from the CIGAR field.
      @param sam_line_obj : an hashage from a line in SAM (See 
                            parseSamLine)
      @return             : [S,S,S,M,M,M,I,M,M,...]
    """
    cigar = sam_line_obj['CIGAR']
      
    cigar=re.sub('([MIDNSHPX=])', r"\1 ", cigar)
    cigar_sub = cigar.split()
    
    status_on_read = ''
    for section in cigar_sub :
        result = re.search("(\d+)([MIDNSHPX=])", section)
        (nb_pos, state)=result.group(1,2)
        
        if( state != "D" and state != "P" and state != "N" ) :
            status_on_read += (state * int(nb_pos)) 
    
    return status_on_read


def getPositionSatus(sam_line_obj):
    """
      Returns the status of each position on read (status : clipping, 
      insertion, match, mismatch and undefined). 
      @param sam_line_obj : an hashage from a line in SAM (See 
                            parseSamLine)
      @return             : ['clipping','match','match','mismatch',...]
    """
    pos   = 0 
    prec  = ""
    isDel = False
    status_by_pos = []
    
    #Retrieve information from the CIGAR field
    cigar_states = cigar2statesRead(sam_line_obj)

    #Process the start clipping
    while(cigar_states[pos] == 'S' or cigar_states[pos] == 'H'):
        status_by_pos.append('clipping')
        pos+=1
    #For each character in the MD field
    for md_char in sam_line_obj['MD'] :
        #If the character is a number
        if( re.search("\d",md_char)):
            #If we were in a deleted area
            if( isDel == True) :
                isDel = False  #End of the deleted area
            prec += str(md_char)
        #If the character is a number in an undeleted area
        elif( re.search("[ATGCNatgcn]", md_char) and  not isDel ):
            if( prec != "" ):
                for i in range(0,int(prec)):
                    #While the position corresponds to an insertion
                    while( cigar_states[pos] == 'I' ):
                        status_by_pos.append('insertion')
                        pos+=1 
                    #Count the match  
                    status_by_pos.append('match')
                    pos+=1 
                prec ="" 
            #While the position corresponds to an insertion
            while( cigar_states[pos] == 'I' ):
                status_by_pos.append('insertion')
                pos+=1 
            #Count the mismatch
            status_by_pos.append('mismatch')
            pos+=1 
       #If the character is the start marker of deletion
        elif(re.search("\^",md_char)):
          #If you had one match or more  before
            if( prec != "" ):
                for i in range(0, int(prec)):
                    #While the position corresponds to an insertion
                    while( cigar_states[pos] == 'I' ):                  
                        status_by_pos.append('insertion')
                        pos+=1 
                    #Count the match
                    status_by_pos.append('match')
                    pos+=1
                prec =""
            #Indicate the start of an deleted area
            isDel = True ;
    #If the MD field ends with matches
    if( prec != "" ):
        for i in range(0, int(prec)):
            #While the position corresponds to an insertion
            while( cigar_states[pos] == 'I'):
                status_by_pos.append('insertion')
                pos+=1 
            #Count the match
            status_by_pos.append('match')
            pos+=1
        prec =""
    #Finish the positions of the read with the information present in the CIGAR field
    while( pos < len(cigar_states) ):
       #The position corresponds to one clipping
        if( cigar_states[pos] == 'S' or cigar_states[pos] == 'H' ):
            status_by_pos.append('clipping')
       #The position corresponds to one insertion
        elif( cigar_states[pos] == 'I' ):
            status_by_pos.append('insertion')
       #The position exists in the CIGAR field but not in the MD when it should be
        else:
            status_by_pos.append('undef')
            sys.stderr.write("[WARNING] " + sam_line_obj['QUERY'] + " the position " + str(pos) + " in the read isn't in MD field.\n") 
        pos+=1 
    #If the query is forward
    if( isForward(sam_line_obj) is False ):
        #Reverse status
        status_by_pos.reverse()

    return status_by_pos


def isAligned(sam_line_obj):
    """
    Returns True if the read in sam_line_obj is aligned.
      @param sam_line_obj : an hashage from a line in SAM (See 
                            parseSamLine)
      @return             : boolean
    """
    return ((sam_line_obj['REFNAME'] != '*') and (sam_line_obj['CIGAR'] != '*'))


def isForward(sam_line_obj):
    """
    Returns True if the strand of the query is forward.
      @param sam_line_obj : an hashage from a line in SAM (See 
                            parseSamLine)
      @return             : boolean
    """
    isForward = False 
    
    if( (int(sam_line_obj['FLAG']) & 16) == 0 ):
        isForward = True
            
    return isForward


def isPaired(sam_line_obj):
    """
    Returns True if the query is paired.
      @param sam_line_obj : an hashage from a line in SAM (See 
                            parseSamLine)
      @return             : boolean
    """
    isPaired = False
    
    if (int(sam_line_obj['FLAG']) & 1) != 0 :
        isPaired = True

    return isPaired


def getReadNumber(sam_line_obj):
    """
    Returns 1 if the query is the first in pair or if it is not paired.
      @param sam_line_obj : an hashage from a line in SAM (See 
                            parseSamLine)
      @return             : int
    """
    readNumber = 1

    if (int(sam_line_obj['FLAG']) & 128) != 0 :
        readNumber = 2

    return readNumber


def parseSamLine(line):
    """
    Returns an hash that contains each field of read line indexed by the
    name of field.
      @param line   : a line of one read extract from a SAM file
      @return       : {QUERY : HXV1, MD : 101, ...}
    """
    line_info = {}
    
    line.rstrip("\n") 
    
    line_subdivisions = line.split( "\t") 
    nb_fields = len(line_subdivisions)
    
    #Mandatory fields    
    line_info = {          
                'QUERY'   : line_subdivisions[0],
                'FLAG'    : line_subdivisions[1],
                'REFNAME' : line_subdivisions[2],
                'CIGAR'   : line_subdivisions[5],
                'SEQ'     : line_subdivisions[9] }
 
    #Optional fields
    regex=re.compile("([^\:]+)\:[^\:]+\:([^\s]*)")
    for i in range (11, nb_fields) : 
        result=re.search(regex,line_subdivisions[i])
        if(result is not None) :
            line_info[result.group(1)] = result.group(2)
        
    return line_info 


def sumStatus(input, nb_by_pos):
    """
      Sum by position the status of reads.
      @param input     : the SAM path.
      @param nb_by_pos : the hash use to write count by position. {read1
                         : [], read2 : []}
    """
    #For each line on SAM file
    for current_line in input :
        line_info = parseSamLine( current_line )

        #If the read is aligned
        if ( isAligned(line_info) ):
            read_number = 1
            status_by_pos = []
            
            #Unperfect alignment
            if( (re.search("^\d+M$", line_info['CIGAR']) is None) or (re.search("^\d+$", line_info['MD']) is None) ):  
                read_number   = getReadNumber(line_info)
                status_by_pos = getPositionSatus(line_info)
            #Perfect alignment
            else:
                read_number   = getReadNumber(line_info)
                seq_length    = len(line_info['SEQ']) 
                status_by_pos = ['match' for e in range( seq_length )] 
            
            #If the read is longer than nb_by_pos
            if len(status_by_pos) > len(nb_by_pos["read"+str(read_number)]) :
                start_pos = len(nb_by_pos["read"+str(read_number)])
                #Create the missing positions
                for pos in range(start_pos, len(status_by_pos) ):
                    nb_by_pos["read"+str(read_number)].append({'mismatch' : 0, 'match' : 0, 'clipping' : 0, 'insertion' : 0, 'undef' : 0})
            
            #Add alignment information to global count
            for pos in range( 0, len(status_by_pos) ):
                nb_by_pos["read"+str(read_number)][pos][status_by_pos[pos]] += 1

    input.close()


def writeOutput (output, separator, *reads):
    """
      Writes one output file from a list of count. 
      @param output    : the output file path.
      @param separator : the separator for each field of one line.
      @param reads     : a list of array. Example of array : 
                         [{match : 2, clipping : 3,...},   #pos 1
                          {match : 4, clipping : 1,...},   #pos 2
                          {match : 4, clipping : 1,...}]   #pos 3
    """
    status_list = 'match', 'mismatch', 'clipping', 'insertion'

    #Find the length of the longest read
    nb_pos = 0
    for read in reads:
        nb_pos = max( nb_pos, len(read) )

    if nb_pos == 0 :
        sys.stderr.write("[WARNING] " + str(output) + " has no reads.\n") 

    #Write output header
    output.write("pos" + separator +  separator.join(status_list) + "\n")

    #For each position of the longest read
    for i in range(0, nb_pos):
        status_at_pos = str(i+1)
        
        #For each status
        for status in status_list:
            nb = 0
            for read in reads:
                if len(read) > i:
                    nb += read[i][status]
            status_at_pos += separator + str(nb)
        
        output.write(status_at_pos + "\n")

    output.close()


def drawOutput (output, precision, *reads):
    """
      Draw one output graph from a list of count. 
      @param output    : the output file path.
      @param precision : The percentage decimal precision.
      @param reads     : a list of array. Example of array : 
                         [{match : 2, clipping : 3,...},   #pos 1
                          {match : 4, clipping : 1,...},   #pos 2
                          {match : 4, clipping : 1,...}]   #pos 3
    """
    colors = {"match": "g", "mismatch": "b", "insertion": "m", "clipping": "r"}
    graph_values = {"match": [], "mismatch": [], "insertion": [], "clipping": []}

    #Find the length of the longest read
    nb_pos = 0
    for read in reads:
        nb_pos = max( nb_pos, len(read) )

    if nb_pos == 0 :
        sys.stderr.write("[WARNING] " + str(output) + " has no reads.\n") 

    #Convert values to percentages by pos
    for pos in range(0, nb_pos):
        read_sum = 0
        
        #Count for each status
        for status in graph_values:
            graph_values[status].append(0)
            
            #For each read
            for read in reads:
                if len(read) > pos:
                    read_sum += read[pos][status]
                    graph_values[status][pos] += read[pos][status]

        #Percentages for each status
        for status in graph_values:
            if read_sum != 0:
                graph_values[status][pos] = (Decimal(graph_values[status][pos])*Decimal(100)/Decimal(read_sum)).quantize(Decimal(10) ** -precision)
            else:
                graph_values[status][pos] = 0
        
    #Build Graph
    keys = list(graph_values.keys())
    keys.sort()
    plot.clf()
    xscale = [e for e in range( 1, nb_pos+1 )]
    fig, (zoom_values, all_values) = plot.subplots(nrows=2, ncols=1)
    fig.subplots_adjust(hspace=0.4, bottom=0.15)

    for status in keys:
        all_values.plot(xscale, graph_values[status], colors[status] + "-", label=status)
        if status != "match":
            zoom_values.plot(xscale, graph_values[status], colors[status] + "-", label=status)
    
    #Set legend
    setSubplotLegend( all_values, 'CigarlineGraph', 0, nb_pos+1 )
    setSubplotLegend( zoom_values, 'CigarlineGraph without Match', 0, nb_pos+1 )
    leg = all_values.legend( keys, loc="upper center", ncol=5, bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True)
    for leg_title in leg.get_texts(): leg_title.set_fontsize('small')
    
    #Create output image
    plot.savefig(output)


def setSubplotLegend (graph, title, x_min, x_max):
    """
      Set the plot legend (axes limits, title and labels).
      @param graph : the plot.
      @param title : the title for plot.
      @param x_min : the x minimum limit of the graph axis.
      @param x_max : the x maximum limit of the graph axis.
    """
    graph.grid()        
    graph.set_title( title )
    graph.set_xlabel( 'Position on read' )
    graph.set_ylabel( 'Read %' )
    graph.set_xlim( x_min, x_max )


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        usage='samtools view INPUT_BAM_FILE | %(prog)s -i - [-o OUT_IMG OUT_IMG_R2 -t OUT_TXT OUT_TXT_R2 --readsplit]',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=('''
Description :
Sum by position the status of reads. This provides for example the number of reads that have a mismatch at first base.
Create percentage graph and/or tabular file with count.

Example :
  - Create one graph for reads 1 and one graph for reads 2
    samtools view hg18.bam | %(prog)s -i - -o r1.png r2.png --readsplit
  - Create one graph and count for reads 1 and one graph and count for reads 2
    samtools view hg18.bam | %(prog)s -i - -o r1.png r2.png -t r1.csv r2.csv --readsplit
  - Create one graph and count for reads
    samtools view hg18.bam | %(prog)s -i - -o all.png -t all.csv'''))
    
    parser.add_argument("-i", "--input", type=argparse.FileType('r'), default=sys.stdin, help="The alignment file. If you use directly a SAM on input (ex : %(prog)s -i hg18.sam -o hg18_all.png) the file does must not contain header section")
    parser.add_argument("-t", "--txt_output", nargs='*', type=argparse.FileType('w'), default=[], help="The file(s) with count")
    parser.add_argument("-o", "--img_output", nargs='*', type=argparse.FileType('w'), default=[], help="The graph(s)")
    parser.add_argument("--readsplit", action="store_true", help="If specified, reads statstics are splitted in two outputs (reads 1 and reads 2)")
    parser.add_argument("--separator", dest="separator", default="\t", help="The fields separator in txt output (default : tabulation)")
    parser.add_argument("-p", "--precision", type=int, dest="precision", default=2, help="The percentage decimal precision in img output (default=2)")
    parser.add_argument("--version", action='version', version=__version__)

    #Getting arguments from the command line
    args = parser.parse_args()

    if len(args.txt_output) == 0 and len(args.img_output) == 0:
        raise ValueError("[ERROR] No output files.")
    if len(args.txt_output)!= 0 and len(args.txt_output)!= 1 and len(args.txt_output)!= 2:
        raise ValueError("[ERROR] Incorrect number of txt output files, expected 1 output file or 2 output files with --readsplit option.")
    if len(args.img_output)!= 0:
        if len(args.img_output)!= 1 and len(args.img_output)!= 2:
            raise ValueError("[ERROR] Incorrect number of img output files, expected 1 output file or 2 output files with --readsplit option.")
        #Import the libraries for graphical output
        import matplotlib
        matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
        import matplotlib.pyplot as plot
        from decimal import Decimal
    if args.readsplit and ( (args.txt_output and len(args.txt_output)!= 2) or (args.img_output and len(args.img_output)!= 2) ):
        raise ValueError("[ERROR] Incorrect number of output files : with --readsplit option you need 2 output files.")
    
    #Init count
    nb_by_pos = {'read1':[],'read2':[]}

    #Sum by position the status of reads
    sumStatus( args.input, nb_by_pos )

    #Write txt output
    if len(args.txt_output)==2 :
        writeOutput( args.txt_output[0], args.separator, nb_by_pos["read1"] )
        writeOutput( args.txt_output[1], args.separator, nb_by_pos["read2"] )
    else:
        if len(args.txt_output)==1:
            writeOutput( args.txt_output[0], args.separator, nb_by_pos["read1"], nb_by_pos["read2"] )
    
    #Write img output
    if len(args.img_output)==2 :
        drawOutput( args.img_output[0], args.precision, nb_by_pos["read1"] )
        drawOutput( args.img_output[1], args.precision, nb_by_pos["read2"] )
    else:
        if len(args.img_output)==1:
            drawOutput( args.img_output[0], args.precision, nb_by_pos["read1"], nb_by_pos["read2"] )
