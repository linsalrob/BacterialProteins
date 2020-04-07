# -*- coding: utf-8 -*-
"""
Extract clusters based on the number of members in the set

This takes the output from cd-hit and extracts protein
sequences into a fasta file. You can specify how many members
in the cluster are required for a sequence to be extracted.

The more members required, the fewer sequences you will get!
"""

import os
import sys
import re
import argparse

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, BacterialProteins'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


ENDC = '\033[0m'
GREEN = '\033[92m'
RED = '\033[91m'


def stream_fasta(fastafile):
    """
    Stream a fasta file, one read at a time. Saves memory!
    Originally part of the roblib libary.

    :param fastafile: The fasta file to stream
    :type fastafile: str
    :param whole_id: Whether to return the whole id (default) or just up to the first white space
    :type whole_id:bool
    :return:A single read
    :rtype:str, str
    """

    try:
        if fastafile.endswith('.gz'):
            f = gzip.open(fastafile, 'rt')
        elif fastafile.endswith('.lrz'):
            f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', fastafile], stdout=subprocess.PIPE).stdout
        else:
            f = open(fastafile, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + fastafile)

    posn = 0
    while f:
        # first line should start with >
        idline = f.readline()
        if not idline:
            break
        if not idline.startswith('>'):
            sys.exit("Do not have a fasta file at: {}".format(idline))
        idline = idline.strip().replace('>', '', 1)
        posn = f.tell()
        line = f.readline()
        seq = ""
        while not line.startswith('>'):
            seq += line.strip()
            posn = f.tell()
            line = f.readline()
            if not line:
                break
        f.seek(posn)
        yield idline, seq


def get_clusters(clstrfile, verbose=False):
    """
    count the members of each cluster and get the exemplar sequence
    :param clstrfile: the cluster input file
    :param verbose: more output
    :return: a dict of tples. The key is cluster name, value is (num members, exemplar)
    """

    if verbose:
        sys.stderr.write(f"{GREEN}Reading clusters from {clstrfile}{ENDC}\n")
    
    clstrs = {}
    lc = None
    c = 0
    exemplar = None
    srch = re.compile('>\w+')
    with open(clstrfile, 'r') as f:
        for l in f:
            l = l.strip()
            if l.startswith('>'):
                if lc:
                    if not exemplar:
                        sys.stderr.write(f"{RED}No exemplar found for {lc}{ENDC}\n")
                        exemplar = ""
                    clstrs[lc] = (c, exemplar)
                c = 0
                exemplar = None
                lc = l.replace('>', '', 1)
            else:
                c += 1
                if l.endswith('*'):
                    m = srch.search(l)
                    exemplar = m.group(0).replace('>', '', 1)

    clstrs[lc] = (c, exemplar)
    return clstrs

def print_sequences(fastafile, clstrs, minn, verbose):
    """
    Print out the fasta sequences
    """

    # filter the clusters for those we want
    if verbose:
        sys.stderr.write(f"{GREEN}Filtering clusters for matches{ENDC}\n")
    
    wanted = set(clstrs[y][1] for y in (filter(lambda x: clstrs[x][0] > minn, clstrs)))

    if verbose:
        sys.stderr.write(f"{GREEN}Streaming fasta for {len(wanted)} sequences{ENDC}\n")

    for seqid, seq in stream_fasta(fastafile):
        if seqid in wanted:
            print(f">{seqid}\n{seq}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fasta file to extract from', required=True)
    parser.add_argument('-c', help='cd-hit clstr file to extract from', required=True)
    parser.add_argument('-m', required=True, type=int,
                       help="minimum number of members in a cluster before extracting (will be > than this)")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    clst = get_clusters(args.c, args.v)
    print_sequences(args.f, clst, args.m, args.v)
