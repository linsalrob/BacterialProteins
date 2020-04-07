# -*- coding: utf-8 -*-
"""
Read one or more protein fasta files and calculate md5sums.

We accept a single fasta file, multiple fasta files, or a directory
of fasta files, read each of those files, generate the md5sum
from the protein sequences, and store the sequences based on 
their md5sum.

We also write a (user-definable) id file that contains the
md5sum and the fasta definition line (separated by a tab).

This code is released under the MIT Licence and is
provided as-is



"""

import os
import sys
import argparse
import hashlib

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, BacterialProteins'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


# colors to make the terminal pretty
GREEN = '\033[92m'
ENDC = '\033[0m'


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



def read_fasta(fafile, idmapfile, minlen=100, verbose=False):
    """
    Read a fasta file and return a dict of proteins and their md5 sums
    :param fafile: fasta file
    :param minlen: minimum protein sequence length (in amino acids) to be included
    :param verbose: more output
    :return:
    """

    if verbose:
        sys.stderr.write(f"{GREEN}Reading {fafile}{ENDC}\n")

    seqs = {}
    with open(idmapfile, 'a') as out:
        for seqid, seq in stream_fasta(fafile):
            # ignore a sequence with a stop codon
            if '*' in seq:
                continue
            if len(seq) < minlen:
                continue
            md5 = hashlib.md5(seq.upper().encode('utf-8')).hexdigest()
            seqs[md5] = seq
            out.write(f"{md5}\t{seqid}\n")

    return seqs





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='file', nargs='*')
    parser.add_argument('-d', help='directory of fasta files')
    parser.add_argument('-i', help='id map file to write')
    parser.add_argument('-m', help='minimum length of protein sequence to include (in amino acids). Default = 300', type=int, default=100)
    parser.add_argument('-o', help='output file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    seqs = {}
    for f in args.f:
        seqs.update(read_fasta(f, args.i, args.m, args.v))

    if args.d:
        for f in os.listdir(args.d):
            seqs.update(read_fasta(os.path.join(args.d, f), args.i, args.m, args.v))

    with open(args.o, 'w') as out:
        for m in seqs:
            out.write(f">{m}\n{seqs[m]}\n")
