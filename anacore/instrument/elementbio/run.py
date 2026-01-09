# -*- coding: utf-8 -*-

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

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

    def isCopied(self):
        """
        Return True if copy from temp to output folder is ended by instrument. Copy step is followed post processing step on instrument.

        :return: True if copy from temp to output folder is ended by instrument.
        :rtype: boolean
        """
        is_completed = False
        if self.isSequenced():
            completed_copy_file = os.path.join(self.path, "RunUploaded.json")  # https://docs.elembio.io/docs/elembio-cloud/runs/seq-run-output/
            is_completed = os.path.exists(completed_copy_file)
        return is_completed

    def isEnded(self):
        """
        Return True if all instrument steps are ended.

        :return: True if all instrument steps are ended.
        :rtype: boolean
        """
        return self.isCopied()

    def isSequenced(self):
        """
        Return True if sequencing step is ended on instrument.

        :return: True if sequencing step is ended on instrument.
        :rtype: boolean
        """
        complete_file = os.path.join(self.path, "AvitiRunStats.json")
        return os.path.exists(complete_file)
