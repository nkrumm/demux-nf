#!/usr/bin/env python3
"""
Parse Illumina Sample Sheets.
Nik Krumm, 2020.

Example usage:

python convert_samplesheet.py
    --input ${samplesheet}
    --is-umi ${params.is_umi}
    --fwd-adapter ${params.fwd_adapter}
    --rev-adapter ${params.rev_adapter}
    --project-name ${params.project_name}
    --library-type ${params.library_type}
    > config.json

"""
import sys
import json
import logging
import argparse
import pandas as pd
from collections import defaultdict
from sample_sheet import SampleSheet, Sample

log = logging.getLogger()

## Utility functions
def neseted_dd():
    return defaultdict(neseted_dd)

def defaultdict_to_regular(d):
    if isinstance(d, defaultdict):
        d = {k: defaultdict_to_regular(v) for k, v in d.items()}
    return d

def parse_samplesheet(args):
    ss = SampleSheet(args.input)
    df = pd.DataFrame([s.to_json() for s in ss.samples])
    df["is_umi"] = args.is_umi
    df["fwd_adapter"] = args.fwd_adapter
    df["rev_adapter"] = args.rev_adapter
    if ("Lane" not in df) or (args.merge_lanes):
        log.info("Merging samples across all lanes!")
        df["Lane"] = "all"
    df["Sample_Project"] = args.project_name
    df["library_type"] = args.library_type
    
    if "Sample_ID" in df:
        df["Sample_Name"] = df.Sample_ID
    elif "Sample_Name" in df:
        df["Sample_ID"] = df.Sample_Name
    else:
        log.error("Samplesheet must specify Sample_ID or Sample_Name!")
        sys.exit(1)

    df = df.drop_duplicates() # needed in case merge_lanes is true, which can result in duplicates
    return df

def build_config(df, args):
    config = neseted_dd()
    cols = ['Sample_Name', 'index', 'index2', 
            "is_umi", "fwd_adapter", "rev_adapter", 
            "library_type"]

    for ix, row in df.iterrows():
        config["lanes"][row.Lane][row.Sample_Project][row.Sample_ID] = row[cols].to_dict()
    
    if args.is_umi:
        config["basemask"] = args.basemask

    config["number_of_lanes"] = df.Lane.nunique()

    return defaultdict_to_regular(config)

def build_samplesheet(df, args):
    samplesheet = SampleSheet()
    for ix, row in df.iterrows():
        s = {
            'Sample_ID': row.Sample_ID,
            'Sample_Name': row.Sample_Name,
            'Sample_Project': row.Sample_Project,
            'index': row["index"],
            'index2': row.index2
        }
        if not args.merge_lanes:
            s["Lane"] = row.Lane
        samplesheet.add_sample(Sample(s))
    return samplesheet

def main(args):
    df = parse_samplesheet(args)
    
    config = build_config(df, args)
    samplesheet = build_samplesheet(df, args)
    
    with open(f"{args.output}.config.json", "w") as config_fp:
        json.dump(config, config_fp, indent=2)
    with open(f"{args.output}.samplesheet.csv", "w") as samplesheet_fp:
        samplesheet.write(samplesheet_fp)
    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', help='Input sample sheet to parse.')
    parser.add_argument('--output', help='Output prefix for files')
    parser.add_argument('--is-umi', action='store_true', default=False)
    parser.add_argument('--merge-lanes', action='store_true', default=False)
    parser.add_argument('--basemask', help='Basemask string to use for UMI demultiplex', default='Y146,I8Y9,I8,Y146')
    parser.add_argument('--fwd-adapter', help='Forward adapter to trim.')
    parser.add_argument('--rev-adapter', help='Reverse adapter to trim.')
    parser.add_argument('--project-name', help='Project name', default='unknown')
    parser.add_argument('--library-type', help='Library type', default='unknown')
    parser.set_defaults(func=main)
    args = parser.parse_args()

    

    #try:
    args.func(args)
    # except AttributeError:
    #     parser.print_help(sys.stderr)
    #     sys.exit(1)        

