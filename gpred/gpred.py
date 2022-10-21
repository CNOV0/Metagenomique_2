#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform prediction of gene."""

import argparse
import sys
import os
import csv
import re


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes (shine box not included).")
    parser.add_argument('-p', dest='predicted_genes_file', type=str,
                        default=os.curdir + os.sep +"predict_genes.csv",
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=str,
                        default=os.curdir + os.sep + "genes.fna",
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file):
    """Extract the complete genome sequence as a single string
    Parameters
    ----------
    fasta_file:
    string
        Name of the fasta file
    Returns
    -------
    String
        Sequence in capital letters
    """
    isfile(fasta_file)
    with open(fasta_file, "rt") as file_in:
        for line in file_in:
            if line.startswith('>'):
                seq = ""
            else:
                seq+=line.strip().upper()
    return seq


def find_start(start_regex, sequence, start, stop):
    """
    Find the start codon
    Parameters
    ----------
    start_regex:
    regex object
        identifier un codon start
    sequence:
    String
        Séquence du génome
    start :
    Int
        position de début de la recherche
    stop :
    Int
        position de fin de la recherche
    Returns
    -------
        Position of the first occurrence of a start codon in the search zone or None
    """
    first_occur =start_regex.search(sequence, start, stop)
    if first_occur is None:
        return first_occur
    return first_occur.start(0)


def find_stop(stop_regex, sequence, start):
    """
    Find the stop codon
    Parameters
    ----------
    stop_regex:
    regex object
        identifier un codon start
    sequence:
    String
        Sequence du genome
    start :
    Int
        position de début de la recherche
    Returns
    -------
        Position of the first occurrence of a stop codon in the search zone or None
    """
    multi_occur = stop_regex.finditer(sequence, start, len(sequence))
    for match in multi_occur:
        if match is not None:
            if (match.start(0) - start )% 3 == 0:
                return match.start(0)
    return None


def has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
    """Find a shine dalgarno motif before the start codon
    Parameters
    ----------
    shine_regex:
    regex object
        permettant d’identifier une séquence de Shine-Dalagarno
    sequence:
    String
        Sequence du genome
    start:
    Int
        Position of the start of the search
    max_shine_dalgarno_distance :
    Int
        Relative position of the Shine-Dalgarno distance from the start codon
    Returns
    -------
    Boolean
        True if a Shine_Dalgarno motif is found or False
    """
    first_occur = shine_regex.search(sequence, start-max_shine_dalgarno_distance, start-6-2)
    if first_occur is None:
        return False
    return True


def predict_genes(sequence, start_regex, stop_regex, shine_regex,
                  min_gene_len, max_shine_dalgarno_distance, min_gap):
    """Predict most probable genes
    Parameters
    ----------
    sequence :
    String
        Sequence du genome
    start_regex:
    regex object
        to identify a start codon
    stop_regex:
    regex object
        to identify a stop codon
    shine_regex:
    regex object
        to identify a Shine-Dalgarno sequence
    min_gene_len:
    Int
        minimal length of a gene
    max_shine_dalgarno_distance:
    Int
        Relative position of the Shine_Dalgarno sequence
    min_gap:
    Int
        Relative minimal distance between 2 genes
    Returns
    -------
    List of lists
        List of the predicted genes
    """
    list_genes = []
    current_pos = 0
    while len(sequence) - current_pos >= min_gap:
        current_pos = find_start(start_regex, sequence, current_pos, len(sequence))
        if current_pos is not None:
            stop = find_stop(stop_regex, sequence, current_pos)
            if stop is not None:
                if (stop - current_pos) >= min_gene_len:
                    if has_shine_dalgarno(shine_regex, sequence, current_pos,
                                          max_shine_dalgarno_distance):
                        list_genes.append([current_pos+1, stop+3])
                        current_pos = (stop+3) + min_gap
                    else:
                        current_pos = current_pos + 1
                else:
                    current_pos = current_pos + 1
            else:
                current_pos = current_pos + 1
    return list_genes


def write_genes_pos(predicted_genes_file, probable_genes):
    """Write list of gene positions
    """
    try:
        with open(predicted_genes_file, "wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_genes(fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp):
    """Write gene sequence in fasta format
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i, gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep,
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(kmer):
    """Get the reverse complement"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in kmer[::-1]])


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()
    # Let us do magic in 5' to 3'
    sequence = read_fasta(args.genome_file)
    sequence = sequence.replace('U', 'T')
    list_genes = predict_genes(sequence, start_regex, stop_regex, shine_regex,
                  args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)

    # We reverse and complement
    sequence_rc = reverse_complement(sequence)
    list_genes_rc = predict_genes(sequence_rc, start_regex, stop_regex, shine_regex,
                  args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)
    for gene in list_genes_rc:
        gene.reverse()
        gene[0] = len(sequence_rc) - gene[0] +1
        gene[1] = len(sequence_rc) - gene[1] +1
    # Don't forget to uncomment !!!
    # Call these function in the order that you want

    # Call to output functions
    write_genes_pos(args.predicted_genes_file, list_genes)
    write_genes(args.fasta_file, sequence, list_genes, sequence, list_genes_rc)


if __name__ == '__main__':
    main()
