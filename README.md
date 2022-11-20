# kmer-counter

![Version](https://img.shields.io/github/v/tag/dgcamblor/kmer-counter?label=Version)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description

A simple program that harnesses the power of parallel processing to compute the absolute frequencies of overlapping nucleotide words of length k (k-mers) in a FASTA or multi-FASTA file[^1] containing a DNA sequence.

In essence, the whole sequence is split into pieces to be processed in parallel. Each of the processes created by the program is assigned to count the k-mers in their corresponding fragments, and the results are joined together for the final output. This considerably speeds up the counting, specially for larger k values.

## Installation

You only need to have installed `python3`. All the modules that were used belong to the Python standard library, so no further installation is required!

## Usage

This program can easily be executed with:

```{bash}
python3 kmer_counter.py FILENAME KLENGTH
```

Or

```{bash}
./kmer_counter.py FILENAME KLENGTH
```

Where `FILENAME` is the name or path to the FASTA/multi-FASTA file and `KLENGTH` is the length of the k-mers that will be counted. If `KLENGTH` is 3, for example, trinucleotide frequencies will be computed. If k-mers contain soft masked nucleotides (i.e., lowercase ones), they are not taken into account by default.

Additional options of interest are:

- Filter output by a minimum frequency: `-m MINIMUM`
- Include non-nucleotide symbols (and soft masked nucleotides): `-i`
- Store the output in a file: `-o OUTPUT`

The absolute frequencies can be extracted from the output with `grep "*"` for parsing purposes (with `awk`, `$2` contains the k-mers, `$3` contains the counts).

## Background

This project was initially developed as part of an optional competition among students at my MSc in Bioinformatics, in which our task was to design the fastest posible Python program capable of counting k-mers in a DNA sequence. This was the winning program! Some further work has been done to it in order to improve the experience for the user.

This is the first project I'm uploading to a GitHub repository, and one of my main aims with it is to use it as an introduction to working with git.

## Utility

k-mer frequencies have proven to be species-specific, and also tend to differ internally among sequences of different nature (e.g., coding/non-coding). Therefore, they can be used for the purpose of sequence identification; for example as part of a Markov chain approach.
This program is a simple Python approximation to computing them relatively fast (although not as fast as a C/C++ implementation would).

## License

Distributed under the MIT License (for more information, see the `LICENSE` file).

[^1]: Note that frequencies in multi-FASTA files are computed for each FASTA individually and then added up together. As of now, no results are given per individual sequence.
