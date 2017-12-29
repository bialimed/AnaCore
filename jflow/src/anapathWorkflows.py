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
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import time

from jflow.workflow import Workflow


class AnapathWorkflow(Workflow):
    def get_cmpt_by_nameid(self, selected_nameid):
        """
        @summary: Returns the component selected by his nameid.
        @param selected_nameid: [str] The nameid of the component (example: "BWAMem.default").
        @return: [Component] The selected component.
        """
        selected_cpmt = None
        for cmpt in self.components:
            if cmpt.get_nameid() == selected_nameid:
                selected_cpmt = cmpt
        return selected_cpmt


    def set_environment(self, workflow_path):
        """
        @summary: Add bin and lib directories from the workflow folder in PATH and PYTHONPATH.
        @param workflow_path: [str] The path to the workflow directory.
        """
        WORKFLOW_DIR = os.path.dirname(os.path.abspath(workflow_path))
        BIN_DIR = os.path.abspath(os.path.join(WORKFLOW_DIR, "bin"))
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + BIN_DIR
        LIB_DIR = os.path.abspath(os.path.join(WORKFLOW_DIR, "lib"))
        sys.path.append(LIB_DIR)
        if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
        else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR


    def write_log(self, log_path, version):
        """
        @summary: Writes a tiny log for user.
        @param log_path: [str] Path to the log file.
        @param version: [str] Version of the workflow.
        """
        with open(log_path, "w") as FH_log:
            FH_log.write(
                "Workflow={}\n".format(self.__class__.__name__) + \
                "Version={}\n".format(version) + \
                "Parameters={}\n".format(" ".join(sys.argv)) + \
                "Start_time={}\n".format(self.start_time) + \
                "End_time={}\n".format(time.time())
            )
