# -*- coding: utf-8 -*-

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import glob
import os


class Run:
    """Read and provide getters on Run from run folder."""

    def __init__(self, path):
        """
        Build and return an instance of Run.

        :param path: Path to the run folder.
        :type path: str
        :return: The new instance.
        :rtype: Run
        """
        self.path = path

    def isEnded(self):
        """
        Return True if all instrument steps are ended.

        :return: True if all instrument steps are ended.
        :rtype: boolean
        """
        return self.isSequenced()

    def isSequenced(self):
        """
        Return True if sequencing step is ended on instrument.

        :return: True if sequencing step is ended on instrument.
        :rtype: boolean
        """
        run_complete_file = os.path.join(self.path, "final_summary*.txt")
        return len(glob.glob(run_complete_file)) != 0
