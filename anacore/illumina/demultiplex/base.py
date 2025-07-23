# -*- coding: utf-8 -*-
"""Classes and functions required for demultiplex sub-packages."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'


class AbstractDemultStat:
    """
    Reader for demultiplex statistics.

    The attribute samples contains for each sample its ID, these barcodes
    sequences and their counts by lane. This count is a list for 0, 1 and 2
    mismatches.

    :Content example for instance.samples:

        .. highlight:: python
        .. code-block:: python

            [
                {
                    "id": "splA",
                    "barcodes": [
                        {
                            "seq": "ATGC",
                            "lanes": [
                                {
                                    "id": 1,
                                    "ct": [12, 14, 15]
                                }, ...
                            ]
                        }, ...
                    ]
                }, ...
            ]
    """

    def __init__(self, path):
        """
        Build and return an instance of DemultStat from demultiplexing statistics file.

        :param path: Path to the demultiplexing statistics file.
        :type path: str
        :return: The new instance.
        :rtype: DemultStat
        """
        self.samples = None
        self._path = path
        self._parse()

    def _parse(self):
        """Read self._path content and store information in self.samples."""
        raise NotImplementedError

    def expectedBarcodesCounts(self):
        """
        Return number of clusters by barcode of samples.

        :return: Number of clusters by barcode of samples.
        :rtype: dict
        """
        ct_by_bc = dict()
        for spl in self.samples:
            if spl["id"] != "Undetermined":
                for barcode in spl["barcodes"]:
                    nb_reads = sum(sum(lane["ct"]) for lane in barcode["lanes"])
                    ct_by_bc[barcode["seq"]] = nb_reads
        return ct_by_bc

    def samplesCounts(self):
        """
        Return number of clusters by samples ID.

        :return: Number of clusters by samples ID.
        :rtype: dict
        """
        ct_by_spl = dict()
        for spl in self.samples:
            nb_reads = 0
            for barcode in spl["barcodes"]:
                nb_reads += sum(sum(lane["ct"]) for lane in barcode["lanes"])
            ct_by_spl[spl["id"]] = nb_reads
        return ct_by_spl

    @property
    def undetermined(self):
        """
        Return undetermined barcodes.

        :return: Undetermined barcodes. Format: [{"seq": "ATGG+TTC", "lanes": [{"id": 1, "ct": 154782}, {"id": 2, "ct": 255567]}].
        :rtype: list
        """
        raise NotImplementedError

    def undeterminedCounts(self):
        """
        Return number of clusters by undetermined UDI.

        :return: Number of clusters by undetermined UDI.
        :rtype: dict
        """
        ct_by_barcode = dict()
        for barcode in self.undetermined:
            ct_by_barcode[barcode["seq"]] = sum(lane["ct"] for lane in barcode["lanes"])
        return ct_by_barcode

    def unexpectedBarcodes(self, parent_bc_max_dist=2, min_ct=1500, break_rate=0.80):
        """
        Return list of odd undetermined barcodes. Unexpected barcodes:
        - have clusters count greater than min_ct
        - have a sequence with at least parent_bc_max_dist differences with known barcodes. This filters out errors in barcode reading of samples with high number of clusters.
          or
          have clusters count greater than corresponding known barcode count. This identifies errors in barcode from sample sheet.
        - are in group of barcodes with count break against the others (smaller_count_in_group / greater_count_in_other < 1 / break_rate). This filters out noise (continuous values) against significant count (abnormal change value).

        :param parent_bc_max_dist: Maximum number off difference between an undetermined barcode and expected barcodes to consider that it is derived from an expected one and it is ignored in odd analysis.
        :type parent_bc_max_dist: int
        :param min_ct: Undetermined barcodes with fewer than this number of reads are ignored in odd analysis.
        :type min_ct: int
        :param break_rate: Rate between previous and next barcodes count of reads to consider break.
        :type break_rate: float
        :return: Odd undetermined index (list of {"spl": barcode, "ct": count, "parent": nearest_bc, "rate": rate}).
        :rtype: list
        """
        ct_by_bc = self.expectedBarcodesCounts()
        expected_bc = ct_by_bc.keys()
        unexpected_barcodes = list()
        prev_count = None
        exists_ct_break = False
        for barcode, count in sorted(self.undeterminedCounts().items(), key=lambda elt: elt[1]):  # List unknown barcode by ascending count
            if count < min_ct:
                prev_count = count
            else:
                if prev_count is None:
                    prev_count = count
                # Identified phiX
                without_udi = False
                for udi in barcode.split("+"):
                    udi = udi.replace("N", "")
                    if set(udi) == {"G"} or set(udi) == {"A"} or udi == "":
                        without_udi = True
                # For barcodes except phiX
                if not without_udi:
                    nearest_bc = _getNearestBarcode(barcode, expected_bc)
                    if nearest_bc["dist"] > parent_bc_max_dist or ct_by_bc[nearest_bc["bc"]] <= count:
                        # The barcode doesn't seem to be the result of reading errors on one of expected barcodes
                        # or
                        # There's an error on one of the expected barcodes because the derivative contains more read
                        rate = round(prev_count / count, 2)
                        if rate <= break_rate and not exists_ct_break:  # First break in count from lower to higher
                            exists_ct_break = True
                            unexpected_barcodes = list()
                        unexpected_barcodes.append({
                            "spl": barcode, "ct": count, "parent": nearest_bc,
                            "rate": rate
                        })
                        prev_count = count
        # Return
        if not exists_ct_break:
            unexpected_barcodes = list()  # 0 unexpected if it has no break in counts
        else:
            unexpected_barcodes = unexpected_barcodes[::-1]  # Sort by descending count
        return unexpected_barcodes


def _getNearestBarcode(query_bc, db_barcodes):
    """
    Return nearest barcode from list. Query must contain at least same number of index phases and same number of cycles than all barcordes in db_barcodes.

    :query_bc: Barcode to search.
    :query_bc: str
    :db_barcodes: List of barcodes.
    :db_barcodes: list
    :return: Nearest barcode as dict {"bc": str, "dist": int}.
    :rtype: dict
    """
    best_nb_diff = 999
    parts_query_bc = query_bc.split("+")
    len_parts_query_bc = [len(elt) for elt in parts_query_bc]
    match = None
    for curr_db_bc in db_barcodes:
        curr_query = query_bc
        # Check authorized configuration
        if "+" not in query_bc and "+" in curr_db_bc:
            raise Exception("Query barcode {} does not contain all index phases for subject {}.".format(query_bc, curr_db_bc))
        parts_db_bc = curr_db_bc.split("+")
        len_parts_db_bc = [len(elt) for elt in parts_db_bc]
        if len_parts_query_bc != len_parts_db_bc:
            for len_query_part, len_db_part in zip(len_parts_query_bc, len_parts_db_bc):
                if len_query_part < len_db_part:
                    raise Exception("Query barcode {} does not contain sufficient number of cycles for subject {}.".format(query_bc, curr_db_bc))
            # Adjust unknown barcode to match expected barcode lengths
            curr_query = "+".join(
                [bc_part[:bc_len] for bc_part, bc_len in zip(parts_query_bc, len_parts_db_bc)]
            )
        # Count diff
        nb_diff = 0
        for nt_a, nt_b in zip(curr_query, curr_db_bc):
            if nt_a != nt_b:
                nb_diff += 1
        # Select nearest
        if nb_diff < best_nb_diff:
            match = curr_db_bc
            best_nb_diff = nb_diff
    return {"bc": match, "dist": best_nb_diff}
