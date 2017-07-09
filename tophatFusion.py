import gzip
from bioStructures import Reference, Region


def is_gzip(file):
    """
    @return: [bool] True if the file is gziped.
    @param file: [str] Path to processed file.
    """
    is_gzip = None
    FH_input = gzip.open( file )
    try:
        FH_input.readline()
        is_gzip = True
    except:
        is_gzip = False
    finally:
        FH_input.close()
    return is_gzip


class TophatFusionIO:
    def __init__( self, filepath, mode="r" ):
        """
        @param filepath: [str] The filepath.
        @param mode: [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open( filepath, mode + "t" )
        else:
            self.file_handle = open( filepath, mode )
        self.current_line_nb = 0
        self.current_line = None

    def __del__( self ):
        self.close()

    def __iter__( self ):
        for line in self.file_handle:
            self.current_line = line.rstrip()
            self.current_line_nb += 1
            if self.current_line.startswith('#'):
                continue
            try:
                fusion_record = self._parse_line()
            except:
                raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n" +
                               "Line content : " + self.current_line )
            else:
                yield fusion_record

    def close( self ):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None
            self.current_line = None

    def __enter__(self):
        return(self)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _parse_line( self ):
        """
        @summary: Returns a structured record from the TophatFusion current line.
        @return: [dict] The fusion described by the current line.
        """
        fusion, trash_1, contig_a, contig_b, depth_a, depth_b, mate_distances = [elt.strip() for elt in self.current_line.split('@')]
        chrom, break_a, break_b, orientation, nb_splitted_reads, nb_splitted_pairs, nb_pairs_splitted_reads, nb_contradict, base_cover_left, base_cover_right, trash_1 = [field.strip() for field in fusion.split("\t")]
        chrom_a, chrom_b = chrom.split("-")
        break_a = int(break_a)
        break_b = int(break_b)
        strand_a, strand_b = [("+" if elt == "f" else "-") for elt in orientation]

        return {
            "partner_a": Region( break_a, break_a, strand_a, Reference(chrom_a, chrom_a, None) ),
            "partner_b": Region( break_b, break_b, strand_b, Reference(chrom_b, chrom_b, None) ),
            "nb_splitted_reads": int(nb_splitted_reads),
            "nb_splitted_pairs": int(nb_splitted_pairs),
            "nb_pairs_splitted_reads": int(nb_pairs_splitted_reads),
            "nb_contradict": int(nb_contradict),
            "base_cover_left": int(base_cover_left),
            "base_cover_right": int(base_cover_right)            
        }
