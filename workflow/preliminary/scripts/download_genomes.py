#!/usr/bin/env python

import sys
import os
import re
import logging
import argparse
import pandas as pd
import gzip
import urllib.request as ul
from Bio import Entrez
from Bio.SeqIO.FastaIO import SimpleFastaParser

__author__ = "Fong Yoke"
__license__ = "MIT"
__email__ = "bldmyoke@gmail.com"

def search_assembly_db(user_email, search_term, max_retrieval="1000"):
        """Download assembly db from  NCBI
        """

        Entrez.email = user_email
        handle = Entrez.esearch(db = "assembly", term = search_term, retmax = max_retrieval)
        record = Entrez.read(handle)
        return record

def get_assembly_summary(user_email, search_term, logger):
        """Download assembly summary for the search_term taxa
        Return dict
        """
        logger.info(f"Getting metadata from NCBI using account: {user_email} ")

        uid_list = search_assembly_db(user_email, search_term)["IdList"]
        str_uid_list = ",".join(uid_list)
        handle = Entrez.esummary(db = "assembly", id = str_uid_list)
        record = Entrez.read(handle)
        summary = record["DocumentSummarySet"]["DocumentSummary"]
        
        return summary

def download_summary(summary, outdir):
    """ Summary to df format, prior to download
    """
    final_df = []
    for uid, item in enumerate(summary):
        assembly_entry = []
        
        name = item["AssemblyName"]
        name = re.sub(" ", "_", name)
        assembly_entry.append(name)

        url = item["FtpPath_GenBank"]
        label = os.path.basename(url)
        link = os.path.join(url, label + "_genomic.fna.gz")
        assembly_entry.append(link)

        final_df.append(assembly_entry)

    final_df = pd.DataFrame(final_df)
    final_df.columns = ["assembly name", "link"]

    # check assembly names are unique
    dupls = final_df["assembly name"].value_counts() > 1
    dupls = dupls[dupls].index.tolist()

    for dupl in dupls:
        for idx, ind  in enumerate(final_df.loc[final_df["assembly name"] == dupl, :].index.tolist()):
            suffix = str.lower(chr(65 + idx))
            final_df.loc[ind, "assembly name"]  = final_df.loc[ind, "assembly name"] + suffix
        
    return final_df

def download_assemblies(processed_df, outdir, logger):
    """Download assemblies
    """
    failed_files = []
    for uid, item in enumerate(processed_df.to_dict(orient = "records")):
        # sys.stderr.write(item)
        name = item["assembly name"]
        name = re.sub(" ", "_", name)

        url = item["link"]
        # label = os.path.basename(url)
        # link = os.path.join(url, label + "_genomic.fna.gz")

        filename = os.path.join(outdir, f"{name}.fasta")

        try:
            if os.path.exists(filename) and os.path.getsize(filename) > 0:
                logger.info(f"Assembly file {outdir}/{name}.fasta already exists")    
            else:
                logger.info(f"Downloading assembly {name} to {outdir}/{name}.fasta from {url}")
                zip_file = ul.urlopen(url)
               
                with open(filename, "w") as f:
               	    f.write(gzip.decompress(zip_file.read()).decode("UTF8"))
        except:
            failed_files.append(url)
            os.remove(filename)

    # Trace failed files
    if len(failed_files) > 0:
        print(f"Number of failed downloads: {len(failed_files)}")
        for failed in failed_files:
            print(failed)

def get_assembly_features(assembly_file):
    """Get assembly size and number of contigs
    """
    assembly_size = 0
    number_of_contigs = 0

    with open(assembly_file, "r") as f:
        for title, seq in SimpleFastaParser(f):
            number_of_contigs += 1
            assembly_size += len(seq)

    return assembly_size, number_of_contigs
        
def give_summary(summary, outdir):
    """ Create assembly metadata
    The metadata are:
    - assembly accession
    - assembly name
    - species name
    - isolate ID
    - assembly release date
    - submitter
    - scaffold N50
    """
    final_df = []
    for uid, item in enumerate(summary):
        assembly_entry = []
        assembly_entry.append(item["AssemblyAccession"])
        
        name = item["AssemblyName"]
        name = re.sub(" ", "_", name)
        assembly_entry.append(name)
        
        final_df.append(assembly_entry)

    final_df = pd.DataFrame(final_df)
    final_df.columns = ["assembly accession", "assembly name"]

    # check assembly names are unique
    dupls = final_df["assembly name"].value_counts() > 1
    dupls = dupls[dupls].index.tolist()

    for dupl in dupls:
        for idx, ind  in enumerate(final_df.loc[final_df["assembly name"] == dupl, :].index.tolist()):
            suffix = str.lower(chr(65 + idx))
            final_df.loc[ind, "assembly name"] = final_df.loc[ind, "assembly name"] + suffix


    final_df["species name"] = [item["SpeciesName"] for item in summary]
    final_df["isolate ID "] = [item["Biosource"]["Isolate"] for item in summary]
    final_df["assembly release date"] = [item["AsmReleaseDate_GenBank"] for item in summary]
    final_df["submitter"] = [item["SubmitterOrganization"] for item in summary]
    final_df["scaffold N50"] = [item["ScaffoldN50"] for item in summary]

    col_assem_size = []
    col_n_contigs = []

    for name in final_df["assembly name"].tolist():
        assembly_file = os.path.join(outdir, f"{name}.fasta")
        print(assembly_file)
        assembly_size, number_of_contigs = get_assembly_features(assembly_file)
        col_assem_size.append(assembly_size)
        col_n_contigs.append(number_of_contigs)
 
    final_df["assembly_size"] = col_assem_size
    final_df["number_of_contigs"] = col_n_contigs
        
    return final_df

def parse_args(args):
    """Argument parsers
    """
    parser = argparse.ArgumentParser(description = "Download genomes using NCBI entrez and get assembly summary")
    parser.add_argument("-e", "--email_address",
        required = True,
        dest = "email_address",
        help = "Email address (required)")

    parser.add_argument("-t", "--taxa_name",
        required = True,
        dest = "taxa_name",
        help = "Taxa name, E.g. \"Rokubacteria\"")

    parser.add_argument("-o", "--outdir",
        required = True,
        dest = "outdir",
        help = "Output directory")
    if not args:
        parser.print_help(sys.stderr)
        sys.exit()
    return parser.parse_args(args)

def main():
    """Main function
    """
    args = sys.argv[1:]
    args = parse_args(args)

    logger = logging.getLogger("download genomes")
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(asctime)s - %(message)s"))
    logger.addHandler(sh)

    os.makedirs(args.outdir, exist_ok = True)
    
    summary = get_assembly_summary(args.email_address, args.taxa_name, logger)
    processed_summary = download_summary(summary, args.outdir)
    download_assemblies(processed_summary, args.outdir, logger)
    
    df = give_summary(summary, args.outdir)
    df.to_csv(os.path.join(args.outdir, "assembly_summary.csv"), sep = "\t", index = False)

if __name__ == "__main__":
    main()