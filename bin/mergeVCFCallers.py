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
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import numpy
import logging
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)

from anacore.vcf import VCFIO, getAlleleRecord, VCFRecord, HeaderInfoAttr, HeaderFormatAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
def getNewHeaderAttr(args):
    """
    Return renamed and new VCFHeader elements for the merged VCF.

    :param args: The script's parameters.
    :type args: NameSpace
    :return: VCFHeader elements (filter, info, format, samples).
    :rtype: dict
    """
    final_filter = {}
    final_info = {
        "SRC": HeaderInfoAttr(
            "SRC", type="String", number=".", description="Variant callers where the variant is identified. Possible values: {}".format(
                {name: "s" + str(idx) for idx, name in enumerate(args.calling_sources)}
            )
        )
    }
    final_format = {
        "AD": HeaderFormatAttr("AD", type="Integer", number="A", description="Allele Depth"),
        "DP": HeaderFormatAttr("DP", type="Integer", number="1", description="Total Depth"),
        "ADSRC": HeaderFormatAttr("ADSRC", type="Integer", number=".", description="Allele Depth by source"),
        "DPSRC": HeaderFormatAttr("DPSRC", type="Integer", number=".", description="Total Depth by source")
    }
    final_samples = None
    for idx_in, curr_in in enumerate(args.inputs_variants):
        with VCFIO(curr_in) as FH_vcf:
            # Samples
            if final_samples is None:
                final_samples = FH_vcf.samples
            elif FH_vcf.samples != final_samples:
                raise Exception(
                    "The samples in VCF are not the same: {} in {} and {} in {}.".format(
                        final_samples,
                        args.inputs_variants[0],
                        FH_vcf.samples,
                        curr_in
                    )
                )
            # FILTER
            for tag, data in FH_vcf.filter.items():
                new_tag = tag
                if tag not in args.shared_filters:  # Rename filters not based on caller
                    new_tag = "s{}_{}".format(idx_in, tag)
                    data.id = new_tag
                    data.source = args.calling_sources[idx_in]
                final_filter[new_tag] = data
            # INFO
            for tag, data in FH_vcf.info.items():
                if tag == args.annotations_field:
                    final_info[tag] = data
                else:
                    new_tag = "s{}_{}".format(idx_in, tag)
                    data.id = new_tag
                    data.source = args.calling_sources[idx_in]
                    final_info[new_tag] = data
            qual_tag = "s{}_VCQUAL".format(idx_in)
            final_info[qual_tag] = HeaderInfoAttr(qual_tag, type="Float", number="1", description="The variant quality", source=args.calling_sources[idx_in])
            # FORMAT
            for tag, data in FH_vcf.format.items():
                new_tag = "s{}_{}".format(idx_in, tag)
                data.id = new_tag
                data.source = args.calling_sources[idx_in]
                final_format[new_tag] = data
    return {
        "filter": final_filter,
        "info": final_info,
        "format": final_format,
        "samples": final_samples
    }


def getMergedRecords(inputs_variants, calling_sources, annotations_field, shared_filters):
    """
    Merge VCFRecords coming from several variant callers.

    :param inputs_variants: Pathes to the variants files.
    :type inputs_variants: list
    :param calling_sources: Names of the variants callers (in same order as inputs_variants).
    :type calling_sources: list
    :param annotations_field: Field used to store annotations.
    :type annotations_field: str
    :param shared_filters: Filters tags applying to the variant and independent of caller like filters on annotations. These filters are not renamed to add caller ID as suffix.
    :type shared_filters: set
    :return: Merged VCF records.
    :rtype: list
    """
    variant_by_name = {}
    for idx_in, curr_in in enumerate(inputs_variants):
        curr_caller = calling_sources[idx_in]
        with VCFIO(curr_in) as FH_in:
            log.info("Process {}".format(curr_caller))
            for record in FH_in:
                variant_name = record.getName()
                # Extract AD and DP
                support_by_spl = {}
                for spl in FH_in.samples:
                    support_by_spl[spl] = {
                        "AD": record.getAltAD(spl)[0],
                        "DP": record.getDP(spl)
                    }
                # Rename filters
                if record.filter is not None:
                    new_filter = []
                    for tag in record.filter:
                        if tag != "PASS":
                            if tag in shared_filters:  # Rename filters not based on caller
                                new_filter.append(tag)
                            else:
                                new_filter.append("s{}_{}".format(idx_in, tag))
                    record.filter = new_filter
                # Rename INFO
                new_info = {}
                for key, val in record.info.items():
                    if key == annotations_field:
                        new_info[key] = val
                    else:
                        new_info["s{}_{}".format(idx_in, key)] = val
                record.info = new_info
                # Backup quality
                if record.qual is not None:
                    record.info["s{}_VCQUAL".format(idx_in)] = record.qual
                # Rename FORMAT
                record.format = ["s{}_{}".format(idx_in, curr_filter) for curr_filter in record.format]
                for spl_name, spl_info in record.samples.items():
                    renamed_info = {}
                    for key, val in spl_info.items():
                        renamed_info["s{}_{}".format(idx_in, key)] = val
                    record.samples[spl_name] = renamed_info
                # Add to storage
                if variant_name not in variant_by_name:
                    variant_by_name[variant_name] = record
                    # Data source
                    record.info["SRC"] = [curr_caller]
                    # Quality
                    if idx_in != 0:
                        record.qual = None  # For consistency, the quality of the variant comes only from the first caller of the variant
                    # AD and DP by sample (from the first caller finding the variant: callers are in user order)
                    record.format.insert(0, "ADSRC")
                    record.format.insert(0, "DPSRC")
                    record.format.insert(0, "AD")
                    record.format.insert(0, "DP")
                    for spl_name, spl_data in record.samples.items():
                        spl_data["AD"] = [support_by_spl[spl_name]["AD"]]
                        spl_data["DP"] = support_by_spl[spl_name]["DP"]
                        spl_data["ADSRC"] = [support_by_spl[spl_name]["AD"]]
                        spl_data["DPSRC"] = [support_by_spl[spl_name]["DP"]]
                else:
                    prev_variant = variant_by_name[variant_name]
                    prev_variant.info["SRC"].append(curr_caller)
                    # IDs
                    if record.id is not None:
                        prev_ids = prev_variant.id.split(";")
                        prev_ids.extend(record.id.split(";"))
                        prev_ids = sorted(list(set(prev_ids)))
                        prev_variant.id = ";".join(prev_ids)
                    # FILTERS
                    if record.filter is not None:
                        if prev_variant.filter is None:
                            prev_variant.filter = record.filter
                        else:
                            prev_variant.filter = list(set(prev_variant.filter) or set(record.filter))
                    # FORMAT
                    prev_variant.format.extend(record.format)
                    # INFO
                    prev_variant.info.update(record.info)
                    for spl_name, spl_data in prev_variant.samples.items():
                        spl_data.update(record.samples[spl_name])
                        spl_data["ADSRC"].append(support_by_spl[spl_name]["AD"])
                        spl_data["DPSRC"].append(support_by_spl[spl_name]["DP"])
    return variant_by_name.values()


def logACVariance(variants, log):
    """
    Display in log the variance on allele counts (AD and AF) between callers.

    :param variants: Merged VCF records.
    :type variants: list
    :param log: Logger object.
    :type log: logging.Logger
    """
    diff_AD = []
    diff_AF = []
    nb_var = 0
    for record in variants:
        if len(record.info["SRC"]) > 1:
            nb_var += 1
            for spl_name, spl_data in record.samples.items():
                # AD
                retained_AD = spl_data["AD"][0]
                max_diff = 0
                for curr_AD in spl_data["ADSRC"][1:]:
                    max_diff = max(max_diff, abs(retained_AD - curr_AD))
                diff_AD.append(max_diff)
                # AF
                AF = [AD / DP for AD, DP in zip(spl_data["ADSRC"], spl_data["DPSRC"])]
                retained_AF = AF[0]
                max_diff = 0
                for curr_AF in AF[1:]:
                    max_diff = max(max_diff, abs(retained_AF - curr_AF))
                diff_AF.append(max_diff)
    log.info("Differences between retained AF and others callers (without missing): median={:.1%}, upper_quartile={:.1%}, 90_persentile={:.1%} and max={:.1%} on {} variants".format(
        numpy.percentile(diff_AF, 50, interpolation='midpoint'),
        numpy.percentile(diff_AF, 75, interpolation='midpoint'),
        numpy.percentile(diff_AF, 90, interpolation='midpoint'),
        max(diff_AF),
        nb_var
    ))
    log.info("Differences between retained AD and others callers (without missing): median={}, upper_quartile={}, 90_persentile={} and max={} on {} variants".format(
        int(numpy.percentile(diff_AD, 50, interpolation='midpoint')),
        int(numpy.percentile(diff_AD, 75, interpolation='midpoint')),
        int(numpy.percentile(diff_AD, 90, interpolation='midpoint')),
        max(diff_AD),
        nb_var
    ))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Merge VCF coming from different calling on same sample(s). It is strongly recommended to apply this script after standardization and before annotation and filtering/tagging.')
    parser.add_argument('-a', '--annotations-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-s', '--shared-filters', nargs='*', default=["lowAF", "OOT", "homoP", "popAF", "CSQ", "ANN.COLLOC", "ANN.RNA", "ANN.CSQ", "ANN.popAF"], help='Filters tags applying to the variant and independent of caller like filters on annotations. These filters are not renamed to add caller ID as suffix. [Default: %(default)s]')
    parser.add_argument('-c', '--calling-sources', required=True, nargs='+', help='Name of the source in same order of --inputs-variants.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-variants', required=True, nargs='+', help='Path to the variants files coming from different callers (format: VCF). The order determine the which AF and AD are retained: the first caller where it is found in this list.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to the merged variants file (format: VCF).')
    args = parser.parse_args()
    args.shared_filters = set(args.shared_filters)

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Get merged records
    variants = getMergedRecords(args.inputs_variants, args.calling_sources, args.annotations_field, args.shared_filters)

    # Log differences in AF and AD
    logACVariance(variants, log)

    # Write
    with VCFIO(args.output_variants, "w") as FH_out:
        # Header
        new_header = getNewHeaderAttr(args)
        FH_out.samples = new_header["samples"]
        FH_out.info = new_header["info"]
        FH_out.format = new_header["format"]
        FH_out.filter = new_header["filter"]
        FH_out.writeHeader()
        # Records
        for record in sorted(variants, key=lambda record: (record.chrom, record.refStart(), record.refEnd())):
            if record.filter is not None and len(record.filter) == 0:
                record.filter = ["PASS"]
            FH_out.write(record)

    log.info("End of job")
