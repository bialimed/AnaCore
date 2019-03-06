#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
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
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import subprocess


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    if len(sys.argv) == 1:
        sys.exit(
            'usage: {} vep_path [vep_arg ...]\n\nWrap the execution of VEP to prevent bug when input file does not contain any variants.\n'.format(
                sys.argv[0]
            )
        )
    vep_cmd = sys.argv[1:]
    input_file = vep_cmd[vep_cmd.index("--input_file") + 1]
    output_file = vep_cmd[vep_cmd.index("--output_file") + 1]

    # Logger initialisation
    logging.basicConfig(level=logging.INFO, format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.info("Start with command: {}".format(" ".join(sys.argv)))

    # Check if VCF contains records
    is_empty = True
    with open(input_file) as FH_in:
        for line in FH_in:
            if not line.startswith("#"):
                is_empty = False

    # Process
    if is_empty:  # Create empty output if VCF does not contain any records
        log.info("Skip annotation")
        log.warn('[WARN] The input file "{}" does not contain any records.'.format(input_file))
        with open(input.vcf) as FH_in:
            with open(output_file, "w") as FH_out:
                for line in FH_in:
                    FH_out.write(line)
    else:  # Annotate
        log.info("Process annotation")
        subprocess.check_call(vep_cmd)
    log.info("End of process")
