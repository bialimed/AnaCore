#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT-O
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

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import time
import argparse
import subprocess
from subprocess import Popen, PIPE



########################################################################
#
# FUNCTIONS
#
########################################################################
class Cmd:
    """
    @summary: Command wrapper.
    @copyright: FROGS's team INRA.
    """
    def __init__(self, program, description, exec_parameters, version_parameters=None):
        """
        @param exec_parameters: [str] The parameters to execute the program. Two possibles syntaxes.
                                If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                Otherwise the parameters will be added after the program in command line.
        @param version_parameters: [str] The parameters to get the program version. Two possibles syntaxes.
                                   If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                   Otherwise the parameters will be added after the program in command line.
        """
        self.program = program
        self.description = description
        self.exec_parameters = exec_parameters
        self.version_parameters = version_parameters

    def get_cmd(self):
        """
        @summary : Returns the command line.
        @return : [str] The command line.
        """
        cmd = None
        if '##PROGRAM##' in self.exec_parameters:
            cmd = self.exec_parameters.replace('##PROGRAM##', self.program)
        else:
            cmd = self.program + ' ' + self.exec_parameters
        return cmd

    def get_version(self, location='stdout'):
        """
        @summary : Returns the program version number.
        @param location : [str] If the version command returns the version number on 'stdout' or on 'stderr'.
        @return : [str] version number if this is possible, otherwise this method return 'unknown'.
        """
        if self.version_parameters is None:
            return "unknown"
        else:
            try:
                cmd = self.program + ' ' + self.version_parameters
                if '##PROGRAM##' in self.exec_parameters:
                    cmd = self.version_parameters.replace('##PROGRAM##', self.program)
                p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                if location == 'stderr':
                    return stderr.decode('ascii').strip()
                else:
                    return stdout.decode('ascii').strip()
            except:
                raise Exception( "Version cannot be retrieve for the software '" + self.program + "'." )

    def parser(self, log_file):
        """
        @summary : Parse the command results to add information in log_file.
        @log_file : [str] Path to the sample process log file.
        """
        pass

    def submit(self, log_file=None):
        """
        @summary : Launch command, trace this action in log and parse results.
        @log_file : [str] Path to the sample process log file.
        """
        # Log
        if log_file is not None:
            FH_log = Logger( log_file )
            FH_log.write( '# ' + self.description + '\n' )
            FH_log.write( '\tSoftware:\n\t\t' + os.path.basename(self.program) + ' version: ' + self.get_version() + '\n' )
            FH_log.write( '\tCommand:\n\t\t' + self.get_cmd() + '\n' )
            FH_log.write( '\tExecution:\n\t\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
            FH_log.close()
        # Process
        subprocess.check_output( self.get_cmd(), shell=True )
        # Log
        if log_file is not None:
            FH_log = Logger( log_file )
            FH_log.write( '\t\tend:   ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
            FH_log.close()
            # Post-process results
            self.parser(log_file)


class Logger:
    """
    @summary: Log file handler.
    @copyright: FROGS's team INRA.
    """
    def __init__(self, filepath=None):
        """
        @param filepath: [str] The log filepath. [default : STDOUT]
        """
        self.filepath = filepath
        self.file_handle = None
        if self.filepath is not None and self.filepath is not sys.stdout:
            self.file_handle = open( self.filepath, "a" )
        else:
            self.file_handle = sys.stdout

    def __del__(self):
        """
        @summary: Closed file handler when the logger is detroyed.
        """
        self.close()

    def close(self):
        """
        @summary: Closed file handler.
        """
        if self.filepath is not None and self.filepath is not sys.stdout:
            if self.file_handle is not None:
                self.file_handle.close()
                self.file_handle = None

    def write(self, msg):
        """
        @summary: Writes msg on file.
        @param msg: [str] The message to write.
        """
        self.file_handle.write( msg )

    @staticmethod
    def static_write(filepath, msg):
        """
        @summary: Writes msg on file.
        @param filepath: [str] The log filepath. [default : STDOUT]
        @param msg: [str] The message to write.
        """
        if filepath is not None and filepath is not sys.stdout:
            FH_log = open( filepath, "a" )
            FH_log.write( msg )
            FH_log.close()
        else:
            sys.stdout.write( msg )

class TmpFiles:
    """
    @summary: Manager for temporary files.
    @copyright: FROGS's team INRA.
    @note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add( "toto.fasta" )
            ...
            tmp_log = tmpFiles.add( "log.txt" )
            ...
        finaly:
            tmpFiles.deleteAll()
    """
    def __init__(self, tmp_dir, prefix=None):
        """
        @param tmp_dir: [str] The temporary directory path.
        @param prefix: [str] The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        """
        if prefix is None:
            prefix = str(time.time()) + "_" + str(os.getpid())
        self.files = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        @summary: Add a temporary file.
        @param filename: The filename without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def delete(self, filepath):
        """
        @summary: Deletes the specified temporary file.
        @param filepath: [str] The file path to delete.
        """
        self.files.remove(filepath)
        if os.path.exists(filepath): os.remove(filepath)

    def deleteAll(self):
        """
        @summary: Deletes all temporary files.
        """
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)

class SamtoolsIndex(Cmd):
    """
    @summary: Index alignment.
    """
    def __init__(self, in_aln):
        """
        @param in_aln: [str] Path to the alignment file (format: BAM).
        """
        cmd_param = "" + \
            " index" + \
            " " + in_aln

        Cmd.__init__( self,
                      "samtools",
                      "Index alignment.",
                      cmd_param,
                      "--version" )

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        output = Cmd.get_version(self)
        return output.split("\n")[0].split(" ")[1].strip()

class SplitBAMByRG(Cmd):
    """
    @summary: Splits BAM by groups of non-overlapping amplicons.
    """
    def __init__(self, in_design, in_aln, out_pattern):
        """
        @param in_design: [str] Path to the amplicons description file (format: BED).
        @param in_aln: [str] Path to the alignment file (format: BAM).
        @param out_pattern: [str] The path pattern for the outputted alignments files (format: BAM). In this path the keyword "{GP}" is replace by the group name for each group.
        """
        cmd_param = "" + \
            " --input-design " + in_design + \
            " --input-aln " + in_aln + \
            " --output-pattern " + out_pattern

        Cmd.__init__( self,
                      "splitBAMByRG.py",
                      "Splits BAM by groups of non-overlapping amplicons.",
                      cmd_param,
                      "--version" )

class AddRGOnBAM(Cmd):
    """
    @summary: Adds tag on reads by origin (amplicon ID).
    """
    def __init__(self, in_aln, out_aln, platform, sample, library):
        """
        @param in_aln: [str] Path to the alignments file (format: BAM).
        @param out_aln: [str] Path to the outputted alignments file (format: BAM).
        @param platform: [str] Platform/technology used to produce the reads.
        @param sample: [str] Sample. Use pool name where a pool is being sequenced.
        @param library: [str] Library.
        """
        cmd_param = "" + \
            " --pl " + platform + \
            " --sm " + sample + \
            " --lb " + library + \
            " --input-aln " + in_aln + \
            " --output-aln " + out_aln

        Cmd.__init__( self,
                      "addRGOnBAM.py",
                      "Adds tag on reads by origin (amplicon ID).",
                      cmd_param,
                      "--version" )

class GATKHaplotypeCaller(Cmd):
    def __init__(self, in_aln, out_variants, in_reference, in_dbsnp=None, in_intervals=None, min_confidence_calling=30, min_confidence_emitting=30):
        cmd_param = "" + \
            " --analysis_type HaplotypeCaller" + \
            " --standard_min_confidence_threshold_for_calling " + str(min_confidence_calling) + \
            " --standard_min_confidence_threshold_for_emitting " + str(min_confidence_emitting) + \
            " --annotateNDA" + \
            " --dontUseSoftClippedBases" + \
            " --reference_sequence " + in_reference + \
            ("" if in_dbsnp is None else " --dbsnp " + in_dbsnp) + \
            ("" if in_intervals is None else " --intervals " + in_intervals) + \
            " --input_file " + in_aln + \
            " --out " + out_variants + \
            " 2> /dev/null"

        Cmd.__init__( self,
                      "java -Xmx30g -jar /softs/tools/gatk/3.6/GenomeAnalysisTK.jar",
                      "Variant calling.",
                      cmd_param,
                      "--version" )

class FreeBayes(Cmd):
    def __init__(self, in_reference, in_aln, out_variants):
        cmd_param = "" + \
            " --fasta-reference " + in_reference + \
            " " + in_aln + \
            " > " + out_variants + \
            " 2> /dev/null"

        Cmd.__init__( self,
                      "freebayes",
                      "Variant calling.",
                      cmd_param,
                      "--version" )

class VarDictStep1(Cmd):
    """
    @summary: Dicovers variants.
    """
    def __init__(self, in_reference, in_regions, in_aln, out_file, min_AF=0.02):
        """
        @param in_reference: [str] Path to the reference sequences file (format: fasta).
        @param in_regions: [str] Path to the amplicons design (format: BED). Start and end of the amplicons must be with primers.
        @param in_aln: [str] Path to the alignments file (format: BAM).
        @param out_file: [str] Path to the outputted file.
        @param min_AF: [float] The threshold for allele frequency.
        """
        cmd_param = "" + \
            " -f " + str(min_AF) + \
            " -F 0" + \
            " -c 1 -S 2 -E 3 -g 4" + \
            " -b " + in_aln + \
            " -G " + in_reference + \
            " " + in_regions + \
            " > " + out_file

        Cmd.__init__( self,
                      "VarDict",
                      "Dicovers variants.",
                      cmd_param,
                      None )

class VarDictStep2(Cmd):
    """
    @summary: Filters variant on strand bias.
    """
    def __init__(self, in_file, out_file):
        """
        @param in_file: [str] Path to the input file.
        @param out_file: [str] Path to the outputted file.
        """
        cmd_param = "" + \
            " cat " + in_file + " | " + \
            " ##PROGRAM##" + \
            " > " + out_file

        Cmd.__init__( self,
                      "teststrandbias.R",
                      "Filters variants on strand bias.",
                      cmd_param,
                      None )

class VarDictStep3(Cmd):
    """
    @summary: Filters variants and converts to VCF.
    """
    def __init__(self, in_file, out_variants, min_AF=0.02):
        """
        @param in_file: [str] Path to the input file.
        @param out_variants: [str] Path to the outputted file (format: VCF).
        @param min_AF: [float] The threshold for allele frequency.
        """
        cmd_param = "" + \
            " cat " + in_file + " | " + \
            " ##PROGRAM##" + \
            " -a" + \
            " -E" + \
            " -f " + str(min_AF) + \
            " > " + out_variants

        Cmd.__init__( self,
                      "var2vcf_valid.pl",
                      "Filters variants and converts to VCF.",
                      cmd_param,
                      None )

class GatherOverlappingRegions(Cmd):
    """
    @summary: Gathers variants from non-overlapping groups.
    """
    def __init__(self, in_regions, in_variants, in_aln, out_variants):
        """
        @param in_regions: [str] Path to the amplicons design. Start and end of the amplicons must be without primers (format: BED).
        @param in_variants: [str] Path to the variants files (format: VCF).
        @param in_aln: [str] Path to the alignments files (format: BAM). Each alignment file correspond to a variants file.
        @param out_variants: [str] Path to the outputted file (format: VCF).
        """
        cmd_param = "" + \
            " --input-designs " + " ".join(in_regions) + \
            " --input-variants " + " ".join(in_variants) + \
            " --input-aln " + " ".join(in_aln) + \
            " --output-variants " + out_variants

        Cmd.__init__( self,
                      "mergeVCFAmpli.py",
                      "Gathers variants from non-overlapping groups.",
                      cmd_param,
                      "--version" )

class MeltOverlappingRegions(Cmd):
    """
    @summary: Melts all the samples contained in variant file in one sample.
    """
    def __init__(self, spl_name, in_variants, out_variants):
        """
        @param spl_name: [str] Name of the final sample.
        @param in_variants: [str] Path to the variants file with several samples (format: VCF).
        @param out_variants: [str] Path to the outputted file (format: VCF).
        """
        cmd_param = "" + \
            " --new-spl-name '" + spl_name + "'" + \
            " --input-variants " + in_variants + \
            " --output-variants " + out_variants

        Cmd.__init__( self,
                      "meltVCFSamples.py",
                      "Melts all the samples contained in variant file in one sample.",
                      cmd_param,
                      "--version" )

class FilterVCFPrimers(Cmd):
    """
    @summary: Removes variants located on amplicons primers.
    """
    def __init__(self, in_sequences, in_regions, in_variants, out_variants):
        """
        @param in_sequences: [str] Path to the reference sequences file (format: fasta). The reference used to discover variants.
        @param in_regions: [str] Path to the amplicons design with their primers (format: BED). The zone of interest is defined by thickStart and thickEnd. The amplicons must not have any overlap between them.
        @param in_variants: [str] Path to the variants file (format: VCF). This file should be sorted by coordinates otherwise the execution time will be dramatically increased.
        @param out_variants: [str] Path to the outputted variants file (format: VCF).
        """
        cmd_param = "" + \
            " --input-sequences " + in_sequences + \
            " --input-regions " + in_regions + \
            " --input-variants " + in_variants + \
            " --output-variants " + out_variants

        Cmd.__init__( self,
                      "filterVCFPrimers.py",
                      "Removes variants located on amplicons primers.",
                      cmd_param,
                      "--version" )

def filterBED(in_bed, in_names, out_bed, nb_col=None):
    """
    @summary: Filters a BED file with the list of names of regions to keep.
    @param in_bed: [str] Path to the initial file (format: BED).
    @param in_names: [str] Path to the file containing the list of names of the retained regions.
    @param out_bed: [str] Path to the filtered file (format: BED).
    @param nb_col: [int] Number of columns in output.
    """
    # Retrieve retained regions names
    retained_regions = dict()
    with open(in_names) as FH_names:
        for line in FH_names:
            retained_regions[line.strip()] = 1
    # Filter BED
    with open(in_bed) as FH_in:
        with open(out_bed, "w") as FH_out:
            for line in FH_in:
                if line.startswith("browser ") or line.startswith("track ") or line.startswith("#"):
                    FH_out.write( line )
                else:
                    fields = [field.strip() for field in line.split("\t")]
                    if fields[3] in retained_regions:
                        if nb_col is not None:
                            line = "\t".join(fields[:nb_col]) + "\n"
                        FH_out.write( line )

def VarDictFct(in_reference, in_regions, in_aln, out_variants, logger, tmp_file, min_AF=0.02):
    """
    @summary: Dicovers amplicons variants with VarDict.
    @param in_reference: [str] Path to the reference sequences file (format: fasta).
    @param in_regions: [str] Path to the amplicons design (format: BED). Start and end of the amplicons must be with primers.
    @param in_aln: [str] Path to the alignments file (format: BAM).
    @param out_variants: [str] Path to the outputted file (format: VCF).
    @param logger: [Logger] Logger used to trace sub-commands.
    @param tmp_file: [TmpFiles] Temporaries files manager.
    @param min_AF: [float] The threshold for allele frequency.
    """
    out_vardict = tmp.add("vardict.txt")
    out_strand_bias = tmp.add("strdBias.txt")
    VarDictStep1( in_reference, in_regions, in_aln, out_vardict, min_AF ).submit( logger )
    VarDictStep2( out_vardict, out_strand_bias ).submit( logger )
    VarDictStep3( out_strand_bias, out_variants, min_AF ).submit( logger )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Varaint calling on Illumina amplicon sequencing. It use VarDictJava (see: https://github.com/AstraZeneca-NGS/VarDictJava).' )
    parser.add_argument( '-m', '--min-AF', default=0.02, help='Variants with an allele frequency under this value are not emitted. [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_reference = parser.add_argument_group( 'Reference' ) # Reference
    group_reference.add_argument( '-g', '--input-genome', required=True, help='The path to the reference genome (format: fasta).' )
    group_panel = parser.add_argument_group( 'Design' ) # Design
    group_panel.add_argument( '-pi', '--input-design-with-primers', required=True, help='The path to the amplicons design with their primers (format: BED).' )
    group_panel.add_argument( '-po', '--input-design-wout-primers', required=True, help='The path to the amplicons design without their primers (format: BED).' )
    group_panel.add_argument( '-pg', '--input-non-overlapping-design', required=True, help='The path to the list of amplicons (format: TSV). The first column is the ID of the amplicon and the second is the name of the group where the amplicon has no overlap with other amplicons of this group.' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-l', '--library-name', help='The library name (example: "patient1_libA").' )
    group_input.add_argument( '-a', '--input-aln', required=True, help='The path to the alignment file (format: BAM).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-ov', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).' )
    group_output.add_argument( '-ol', '--output-log', default=sys.stdout, help='The path to the outputted log file (format: txt). [Default: STDOUT]' )
    args = parser.parse_args()

    Logger.static_write(args.output_log, "## Application\n\tSoftware:\n\t\t" + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\n\tCommand:\n\t\t" + " ".join(sys.argv) + "\n\n")
    tmp = TmpFiles( os.path.dirname(args.output_variants) )
    library_name = os.path.basename(args.input_aln).split(".")[0] if args.library_name is None else args.library_name

    # Get non-overlapping groups
    groups_names = set()
    with open(args.input_non_overlapping_design) as FH_gp:
        for line in FH_gp:
            if not line.startswith("#"):
                amplicon_id, group_name = [elt.strip() for elt in line.split("\t")]
                groups_names.add( group_name )
    groups_names = list(groups_names)

    # Split BAM in non-overlapping regions
    gp_alignment = [tmp.add(gp + ".bam") for gp in groups_names]
    out_bam_pattern = gp_alignment[-1][:-(len(groups_names[-1]) + 4)] + "{GP}.bam"
    SplitBAMByRG( args.input_non_overlapping_design, args.input_aln, out_bam_pattern ).submit( args.output_log )

    # Variant calling
    groups = list()
    for idx_gp, curr_gp in enumerate(groups_names):
        curr_gp_aln = gp_alignment[idx_gp]

        # Index BAM
        tmp.files.append( curr_gp_aln + ".bai" )
        SamtoolsIndex( curr_gp_aln ).submit( args.output_log )

        # Select regions of current group
        curr_gp_regions = tmp.add( curr_gp + "_amplicons.txt" )
        subprocess.check_output('grep "' + curr_gp + '$" ' + args.input_non_overlapping_design + ' | cut -f 1 > ' +  curr_gp_regions, shell=True )
        curr_gp_regions_with_prim = tmp.add( curr_gp + "_withPrimers.bed" )
        filterBED( args.input_design_with_primers, curr_gp_regions, curr_gp_regions_with_prim )
        curr_gp_regions_with_prim_4_col = tmp.add( curr_gp + "_withPrimers_4col.bed" )
        filterBED( args.input_design_with_primers, curr_gp_regions, curr_gp_regions_with_prim_4_col, 4 )
        curr_gp_regions_wout_prim = tmp.add( curr_gp + "_woutPrimers.bed" )
        filterBED( args.input_design_wout_primers, curr_gp_regions, curr_gp_regions_wout_prim )

        # Add RG on BAM
        curr_gp_aln_new_RG = tmp.add( curr_gp + "_RG.bam" )
        AddRGOnBAM( curr_gp_aln, curr_gp_aln_new_RG, "ILLUMINA", library_name, library_name ).submit( args.output_log )
        tmp.files.append( curr_gp_aln_new_RG + ".bai" )
        SamtoolsIndex( curr_gp_aln_new_RG ).submit( args.output_log )

        # Call variants
        curr_gp_vcf = tmp.add( curr_gp + ".vcf" )
        VarDictFct( args.input_genome, curr_gp_regions_with_prim_4_col, curr_gp_aln_new_RG, curr_gp_vcf, args.output_log, tmp, args.min_AF )

        # Filters variants located on primers
        curr_gp_clean_vcf = tmp.add( curr_gp + "_clean.vcf" )
        FilterVCFPrimers( args.input_genome, curr_gp_regions_with_prim, curr_gp_vcf, curr_gp_clean_vcf ).submit( args.output_log )

        # Store current group files
        groups.append({
            "name": curr_gp,
            "aln": curr_gp_aln,
            "vcf": curr_gp_clean_vcf,
            "design_wout_prim": curr_gp_regions_wout_prim,
            "design_with_prim": curr_gp_regions_with_prim
        })

    # Merge overlapping amplicons
    out_gather = tmp.add( "gatherOverlapping.vcf" )
    GatherOverlappingRegions(
        [curr_gp["design_wout_prim"] for curr_gp in groups],
        [curr_gp["vcf"] for curr_gp in groups],
        [curr_gp["aln"] for curr_gp in groups],
        out_gather
    ).submit( args.output_log )
    MeltOverlappingRegions( library_name, out_gather, args.output_variants ).submit( args.output_log )

    # Clean temporary files
    tmp.deleteAll()
