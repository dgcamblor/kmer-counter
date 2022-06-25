#!/usr/bin/env python3

#------------------------------------------------------------------------------
#                               kmer_counter
#------------------------------------------------------------------------------

# Author: @dgcamblor
# Version: 1.0

import argparse
import os
import sys
import itertools
import threading
import time
import re 
from collections import defaultdict  
from multiprocessing import cpu_count, Pool 


parser = argparse.ArgumentParser(description="Counts k-mers (overlapping \
nucleotide words of length k) in a FASTA/multi-FASTA file storing a DNA \
    sequence")

parser.add_argument("filename", help="Name or path (path/to/file.fa) to the \
    FASTA/multi-FASTA file")
parser.add_argument("klength", help="Length of the k-mers (how much \
nucleotides are the words comprised of) (e.g., KLENGTH = 3 means \
    counting trinucleotides)", type=int)
parser.add_argument("-m", "--minimum", help="Filter output to only words \
    surpassing a minimum frequency", type=int)
parser.add_argument("-i", "--include", help="Include non-nucleotide symbols \
    in the output (N, R/purine, Y/pyrimidine, etc.), and also soft masked \
        nucleotides (acgt)", action="store_true")
parser.add_argument("-o", "--output", help="Store the results in an output \
    file (OUTPUT = name of the file)", type=str)
args = parser.parse_args()

FILENAME = args.filename
KLENGTH = args.klength
MINIMUM = args.minimum if args.minimum else 1

# Declare regexes
HEADER = re.compile(">.*\n")
NEWLINE = re.compile("\n")  
NONNUCLEOTIDE = re.compile("[^ACGT]")

ended = False  # Requirement for the loading bar


def check_input():
    """Check for valid inputs"""
    if KLENGTH <=0 :
        print("ERROR: K-mer length cannot be less than 1")
        sys.exit()
    if not os.path.exists(FILENAME):
        print(f"ERROR: {FILENAME} file does not exist")
        sys.exit()


def loading():
    """Displays a loading bar so that the user doesn't think the program is
    not working (hopefully it will be working)."""
    for c in itertools.cycle(['.  ', '.. ', '...']):
        if ended == True:
            break
        sys.stdout.write('\rProcessing ' + c)
        sys.stdout.flush()
        time.sleep(0.5)


def get_seq_frequencies(sequence, KLENGTH=KLENGTH):
    """Obtain a dictionary of kmer frequencies from a nucleotide sequence. Its
    main purpose is to be called for each fragment of the parallelization 
    process.

    Args:
        sequence (str): nucleotide sequence
        KLENGTH (int): length of kmers. Already specified so that the function
        is called seamlessly inside the parallelization process. 
            Defaults to KLENGTH (as introduced by the user).

    Returns:
        dict: Dictionary of pairs kmer:frequency for the sequence. It may 
        contain other symbols apart from nucleotides.
    """
    frequencies = defaultdict(int)  # Initialize dict with 0s
    
    # (- KLENGTH + 1) avoids endings with length(kmer) < KLENGTH
    for i in range(len(sequence) - KLENGTH + 1): 
        kmer = sequence[i:i+KLENGTH]
        frequencies[kmer] += 1
    
    return frequencies


def get_file_frequencies(filename, ncpu):
    """Get the global count of kmer frequencies in a fasta file. Divides the
    file in chunks to be processed in parallel, computing the frequencies in 
    each of them with get_seq_frequencies() and then uniting them for the
    global frequencies to be returned.

    Args:
        filename (str): the file containing the sequences which frequencies
        will be computed.
        ncpu (int): number of cpus used to perform parallelization

    Returns:
        dict: Dictionary of pairs kmer:frequency for the file. It may 
        contain other symbols apart from nucleotides.
    """
    frequencies = defaultdict(int)  # Initialize dict with 0s

    with open(filename, "r") as F:
        # We aim to split a multiFASTA into its constitutive FASTAs
        fastas = re.split(HEADER, F.read()) 

        # And then split each FASTA into chunks to be parallelly processed
        fasta_chunks = []

        for fasta in fastas:
            if fasta == "":  # Avoid "" items after splitting by header
                continue

            fasta = NEWLINE.sub("", fasta)
            n = len(fasta)//ncpu

            chunks = [fasta[i:i+n+KLENGTH-1] for i in range(0, len(fasta), n)]

            fasta_chunks.extend(chunks)

        # Parallel kmer counting (number of threads = number of cores)
        with Pool(ncpu) as pool:
            counts = pool.map(get_seq_frequencies, fasta_chunks)

        # Unify the counts
        for count in counts:
            for kmer, frequency in count.items():
                frequencies[kmer] += frequency

    return frequencies


def show_results(frequencies):
    """Display a summary of the kmer counting results (or save it to a file
    if the user has specified so)

    Args:
        frequencies (dict): a dictionary with the global frequencies of the 
        file.
    """

    # If the users wants to save the results, change the standard output to 
    # the new file 
    if args.output:
        original_stdout = sys.stdout  # Save the original standard output
        f = open(args.output, "w")
        sys.stdout = f

    print("---------------------------")
    print("       K-mer counter       ")
    print("       By: @dgcamblor      ")
    print("---------------------------")
    print(f"File: {FILENAME}          ")
    print("---------------------------")
    print("    k-mer  /  frequency    ")
    print("---------------------------")

    total = 0  # Store the total number of kmers
    different = 0 # Store the number of different kmers 


    for kmer, frequency in sorted(frequencies.items()):
        # Remove non-nucleotide letters (if --include was not specified and do 
        # not display kmers that had a letter removed.
        if not args.include:
            kmer = NONNUCLEOTIDE.sub("", kmer)  
        if len(kmer) == KLENGTH:
            if frequency >= MINIMUM:  # Filtering output by minimum
                print(f"*  {kmer}  {frequency}")  # Provide * for grep
                total += frequency
                different += 1

    print("----------------------------")
    print(f"{total} k-mers with k = {KLENGTH}")
    print(f"{different} different k-mers")
    print("----------------------------")

    if args.output:
        sys.stdout = original_stdout  # Restore standard output


def main():
    check_input()

    #-------Loading bar requirements---------
    global ended
    t = threading.Thread(target=loading)
    t.daemon = True
    t.start()
    #----------------------------------------

    ncpu = cpu_count()

    file_freqs = get_file_frequencies(FILENAME, ncpu)

    #-------Stop loading bar-----------------      
    ended = True
    sys.stdout.write('\rDone!         \n\n')
    sys.stdout.flush()
    #----------------------------------------

    show_results(file_freqs)


if __name__ == "__main__":
    main()