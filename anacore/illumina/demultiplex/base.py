# -*- coding: utf-8 -*-
"""Classes and functions required for demultiplex sub-packages."""


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

    def unexpectedBarcodes(self, skipped_spl=None):
        """
        Return UDI with a number of clusters greater than or equal to the smallest sample. Cluster without UDI are excluded (detected by repeat of only A or G in index).

        :param skipped_spl: List of samples excluded from the smallest sample finding. This can be used to exclude negative controls without reads.
        :type skipped_spl: set
        :return: UDI with a number of clusters greater than or equal to the smallest sample: {"udi": count, ...}.
        :rtype: dict
        """
        unexpected_barcodes = list()
        # Get smallest sample count
        if skipped_spl is None:
            skipped_spl = set()
        ct_by_identified = {spl: ct for spl, ct in self.samplesCounts().items() if spl not in skipped_spl}
        smallest_spl = min(ct_by_identified, key=ct_by_identified.get)
        identified_min_ct = ct_by_identified[smallest_spl]
        # Get undetermined barcodes with more read than smallest sample
        for barcode, count in self.undeterminedCounts().items():
            if count >= identified_min_ct:
                without_udi = False
                for udi in barcode.split("+"):
                    if set(udi) == {"G"} or set(udi) == {"A"}:
                        without_udi = True
                if not without_udi:
                    unexpected_barcodes.append({"spl": barcode, "ct": count})
        return sorted(unexpected_barcodes, key=lambda elt: elt["ct"], reverse=True)
