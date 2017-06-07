class BEDRecord:
    def __init__(self, chrom=None, start=None, end=None, name=None, score=None, strand=None, thickStart=None, thickEnd=None, itemRgb=None, blockCount=None, blockSizes=None, blockStarts=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts 

class BED:
    def __init__( self, filepath, mode="r" ):
        """
        @param filepath : [str] The filepath.
        @param mode : [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        self.file_handle = open( filepath, mode )
        self.current_line_nb = 0
        self.current_line = None
    
    def __del__( self ):
        self.close()
    
    def __iter__( self ):
        for line in self.file_handle:
            self.current_line = line.rstrip()
            self.current_line_nb += 1
            bed_record = None
            try:
                bed_record = self._parse_line()
            except:
                raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               "Line content : " + self.current_line )
            else:
                yield bed_record

    def close( self ) :
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None
            self.current_line = None

    def _parse_line( self ):
        """
        @summary : Returns a structured record from the BED current line.
        @return : [BEDrecord]
        """
        fields = [elt.strip() for elt in self.current_line.split('\t')]
        record = BEDrecord(fields)
        return record
