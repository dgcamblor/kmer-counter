#!/usr/bin/env python3

#------------------------------------------------------------------------------
# kmer_counter.py
#------------------------------------------------------------------------------

__author__ = "@dgcamblor"
__version__ = "1.1"
__license__ = "MIT"

import argparse
import os
import sys
import itertools
import threading
import time
import re 
from collections import defaultdict  
from multiprocessing import cpu_count, Process, Pipe

parser = argparse.ArgumentParser(description="Counts k-mers (overlapping \
    nucleotide words of length k) in a FASTA/multi-FASTA file storing a DNA \
    sequence")

parser.add_argument("filename", help="Name or path (path/to/file.fa) to the \
    FASTA/multi-FASTA file")
parser.add_argument("klength", help="Length of the k-mers (how much \
    nucleotides are the words comprised of) (e.g., KLENGTH = 3 means \
    counting trinucleotides)", type=int)

parser.add_argument("-v", "--version", action="version", version="%(prog)s " + 
    __version__)
parser.add_argument("-m", "--minimum", help="Filter output to only words \
    surpassing a minimum frequency", type=int)
parser.add_argument("-i", "--include", help="Include non-nucleotide symbols \
    in the output (N, R/purine, Y/pyrimidine, etc.), and also soft masked \
    nucleotides (acgt)", action="store_true")
parser.add_argument("-o", "--output", help="Store the results in an output \
    file (OUTPUT = name of the file)", type=str)
parser.add_argument("-x", "--processes", help="Number of processes to use \
    (default: number of cores)", type=int)

args = parser.parse_args()

# Global variables
FILENAME = args.filename
KLENGTH = args.klength
MINIMUM = args.minimum if args.minimum else 1
NCPU = args.processes if args.processes else cpu_count()

# Declare regexes
HEADER = re.compile(">.*\n")
NEWLINE = re.compile("\n")  
NONNUCLEOTIDE = re.compile("[^ACGT]")

# Requirement for the loading bar, tracks if the program has finished
ended = False  


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
        sys.stdout.write('\rCounting k-mers' + c)
        sys.stdout.flush()
        time.sleep(0.5)


def loading_start():
    """Starts the loading bar"""
    t = threading.Thread(target=loading)
    t.daemon = True
    t.start()


def loading_stop():
    """Stops the loading bar to display the results"""
    global ended
    ended = True
    sys.stdout.write('\r' + ' ' * 18 + '\n')  # Clear the loading bar
    sys.stdout.flush()


def get_seq_frequencies(sequences, conn, KLENGTH=KLENGTH):
    """Obtain a dictionary of kmer frequencies from a nucleotide sequence. Its
    main purpose is to be called for each fragment of the parallelization 
    process.

    Args:
        sequence s(str): nucleotide sequences assigned to the current process.
        conn (multiprocessing.Pipe): pipe to send the results to the main
            process.
        KLENGTH (int): length of kmers. Already specified so that the function
        is called seamlessly inside the parallelization process. 
            Defaults to KLENGTH (as introduced by the user).

    Returns:
        dict: Dictionary of pairs kmer:frequency for the sequence. It may 
        contain other symbols apart from nucleotides.
    """
    frequencies = defaultdict(int)  # Initialize dict with 0s

    for sequence in sequences:
        # (- KLENGTH + 1) avoids endings with length(kmer) < KLENGTH
        for i in range(len(sequence) - KLENGTH + 1): 
            kmer = sequence[i:i+KLENGTH]
            frequencies[kmer] += 1
    
    conn.send(frequencies)


def get_file_frequencies(filename, ncpu):
    """Get the global count of kmer frequencies in a fasta file. Divides the
    file in chunks to be processed in parallel, computing the frequencies in 
    each of them with get_seq_frequencies() and then uniting them for the
    global frequencies to be returned.

    Args:
        filename (str): the file containing the sequences which frequencies
        will be computed.
        ncpu (int): number of cpus/processes used to perform parallelization

    Returns:
        dict: Dictionary of pairs kmer:frequency for the file. It may 
        contain other symbols apart from nucleotides.
    """
    with open(filename, "r") as F:
        # We aim to split a multiFASTA into its constitutive FASTAs
        fastas = re.split(HEADER, F.read()) 

        # And then split each FASTA into chunks to be parallelly processed
        fasta_chunks = []

        for fasta in fastas:
            if fasta == "":  # Avoid "" items after splitting by header
                continue

            fasta = NEWLINE.sub("", fasta)
            fasta_len = len(fasta)

            n_chunk = fasta_len//ncpu  # Number of nucleotides per chunk 
            remainder = fasta_len % ncpu  # Remainder of nucleotides

            if n_chunk == 0: 
                loading_stop()
                print("ERROR: sequence in FASTA file is too short to be "  \
                      " processed with the current number of processes, " \
                        "try lowering it with -x")
                sys.exit()

            # Split the FASTA into chunks of n_chunk nucleotides, taking 
            # into account the remainder
            chunks = []
            start = 0
            overlap = KLENGTH - 1  # Overlap between chunks

            for _ in range(ncpu):
                end = start + n_chunk + overlap
                if remainder > 0:
                    end += 1
                    remainder -= 1
                chunks.append(fasta[start:end])
                start = end - overlap

            fasta_chunks.extend(chunks)

    # Now we can parallelize the computation of frequencies in each chunk
    conn_list = []
    proc_list = []

    for i_proc in range(ncpu):
        task = [fasta_chunks[i_chunk] for i_chunk in range(i_proc, len(fasta_chunks), ncpu)]
        
        parent_conn, child_conn = Pipe()  # Sending results to main process
        conn_list.append(parent_conn)

        p = Process(target=get_seq_frequencies, args=(task, child_conn))
        proc_list.append(p)
        
        proc_list[i_proc].start()

    # Wait for all processes to finish (when they send their results)
    while True:
        if all([conn.poll() for conn in conn_list]):
            break

    # Unify the counts
    file_frequencies = defaultdict(int)  # Initialize dict with 0s

    for conn in conn_list:
        proc_frequncies = conn.recv()
        for kmer, count in proc_frequncies.items():
            file_frequencies[kmer] += count

    return file_frequencies


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
    print(f"       By: {__author__}      ")
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
                print(f"*\t{kmer}\t{frequency}")  # Provide * for grep
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

    loading_start()
    file_freqs = get_file_frequencies(FILENAME, NCPU)
    loading_stop()

    show_results(file_freqs)


if __name__ == "__main__":
    main()