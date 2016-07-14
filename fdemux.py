#!/usr/bin/env python3
#
# fdemux is a FASTA/FASTQ demultiplexer with support for (arbitrarily) fuzzy barcode matching,
# asymmetric barcoding, and parallel processing.
#
# Author: Ted Verhey (verheytb@gmail.com)

import os, sys, argparse, multiprocessing
from datetime import datetime
from Bio import SeqIO


class Counter(object):
    # a multiprocessing-safe counter
    def __init__(self, ):
        self.val = multiprocessing.Value('i', 0)

    def increment(self, n=1):
        with self.val.get_lock():
            self.val.value += n

    @property
    def value(self):
        return self.val.value


def printMessage(contents, ontop=False):
    if ontop:
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "   " + contents, end="\r", flush=True)
    else:
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "   " + contents)


def fuzzySubstring(needle, haystack):
    # Calculates the fuzzy match of needle in haystack,
    # using a modified version of the Levenshtein distance
    # algorithm.
    # The function is modified from the levenshtein function
    # in the bktree module by Adam Hupp
    m, n = len(needle), len(haystack)

    # base cases
    if m == 1:
        return not needle in haystack
    if not n:
        return m

    row1 = [0] * (n + 1)
    for i in range(0, m):
        row2 = [i + 1]
        for j in range(0, n):
            cost = (needle[i] != haystack[j])

            row2.append(min(row1[j + 1] + 1,  # deletion
                            row2[j] + 1,  # insertion
                            row1[j] + cost)  # substitution
                        )
        row1 = row2
    return min(row1)


def worker(sequence):
    # checks sense strand
    results = []
    for fcode, rcode in barcodes:
        if fuzzySubstring(fcode.seq, sequence.seq) <= args.fuzziness and\
                        fuzzySubstring(rcode.seq, sequence.seq) <= args.fuzziness:
            results.append([fcode, rcode])

    # checks antisense strand
    rcresults = []
    rcsequence = sequence.reverse_complement(id=sequence.id + "_rc", name="", description="")
    for fcode, rcode in barcodes:
        if fuzzySubstring(fcode.seq, rcsequence.seq) <= args.fuzziness and\
                        fuzzySubstring(rcode.seq, rcsequence.seq) <= args.fuzziness:
            rcresults.append([fcode, rcode])

    counter.increment()
    # returns something only if there is one barcode pair detected on both strands.
    if len(results) == 1 and len(rcresults) == 1:  # check for permissible rotational symmetry
        if results[0][0].seq == rcresults[0][0].seq and results[0][1].seq == rcresults[0][1].seq:
            return (sequence, results[0][0].id, results[0][1].id)
    elif len(results) == 1 and len(rcresults) == 0:
        return (sequence, results[0][0].id, results[0][1].id)
    elif len(results) == 0 and len(rcresults) == 1:
        return (rcsequence, rcresults[0][0].id, rcresults[0][1].id)
    else:  # reads are skipped if there is anything other than one barcode pair identified from both strands.
        skipCounter.increment()
        return None


# parses the arguments and displays usage information
parser = argparse.ArgumentParser(
    description="DEMUXP is a parallelized demultiplexer for sequences in a FASTA/FASTQ file that sorts them by barcode. DEMUXP supports assymmetric and symmetric barcoding schemes.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input_file", type=str, required=True, help="FASTA/FASTQ file containing sequences to demultiplex")
parser.add_argument("-b", "--barcodes_file", type=str, required=True,
                    help="FASTA file containing the barcodes. Must be in the PacBio order, ie. F1, R1, F2, R2, F3, R3, etc.")
parser.add_argument("-o", "--output_directory", type=str,
                    help="A directory for demultiplexed file output")
parser.add_argument("-z", "--fuzziness", type=int, default=0,
                    help="The total number of differences (including insertions, deletions, and substitutions) allowed when identifying a barcode in a sequence. This is important for identifying barcodes in sequences with high error rates (ie. PacBio, Oxford Nanopore). If set too high, multiple barcode pairs will be identified for each read and demultiplexing errors, skipped reads will increase.")
parser.add_argument("-c", "--cpus", default=multiprocessing.cpu_count(), type=int,
                    help="Number of processes to use")
args = parser.parse_args()

if any(args.input_file.lower().endswith(x) for x in {".fastq", ".fnq", ".fq"}):
    filetype = "fastq"
elif any(args.input_file.lower().endswith(x) for x in {".fasta", ".fna", ".fa"}):
    filetype = "fasta"
else:
    sys.exit("Error: Input file must be in a FASTA or FASTQ format.")

if not args.output_directory:
    outputdir = os.path.join(args.input_file[:args.input_file.rfind("/")],
                             "demuxp_" + datetime.now().strftime("%Y%m%d-%H%M%S"))
else:
    outputdir = args.output_directory

if not os.path.exists(outputdir):
    os.makedirs(outputdir)

# Reads the barcodes into memory
with open(args.barcodes_file, "r") as b:
    barcodeseqlist = list(SeqIO.parse(b, "fasta"))
    if len(barcodeseqlist) % 2 != 0:
        sys.exit("Barcode FASTA file contains an odd number of sequences; check the file.")
    barcodes = list(zip(barcodeseqlist[::2], barcodeseqlist[1::2]))

filesDict = {}
counter = Counter()
skipCounter = Counter()
p = multiprocessing.Pool()
with open(args.input_file, "r") as inHandle:
    for output in p.imap(worker, SeqIO.parse(inHandle, filetype)):
        printMessage("%d sequences completed, %d skipped" % (counter.value, skipCounter.value), ontop=True)
        if output:
            outSequence, FB, RB = output
            if not FB + "&" + RB in filesDict:
                filesDict[FB + "&" + RB] = open(os.path.join(outputdir, FB + "--" + RB + "." + filetype), "a")
            SeqIO.write(outSequence, filesDict[FB + "&" + RB], filetype)
for filehandle in filesDict.values():
    filehandle.close()
printMessage("Finished processing all %d sequences (%d skipped)." % (counter.value, skipCounter.value), ontop=False)
