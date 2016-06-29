/*
 *  alignment.h
 *
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License.
 *
 */

#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <vector>
#include <string>
#include <ostream>
#include "tree.h"


// This struct encodes the data of an alignment. For each species, the data is stored
// as a vector of letters (seq). The letters are stored as an integer between
// 0 and nalpha-1. In the vector alpha there is a translation to actual
// characters (A,C,G,T for example) used for printing.
struct Alignment {
  long nspecies;                       // number of species
  long nalpha;                         // number of letters
  long len;                            // length of the sequence
  std::vector<char> alpha;             // the actual letters
  std::vector<std::vector<int> > seq;  // a sequence for a given species
  std::vector<std::string> name;       // species names
};

// This struct stores the counts of patterns in an alignment in the vector c.
// Every pattern (or state) is encoded by a long integer, as described in state.h.
struct Counts {
  std::vector<std::string> species;    // Name of the species.
  std::vector<double> c;               // Vector of for every pattern. Has length nstates
  double N;                            // Sum of all the counts
  long nstates;                        // nstates = nalpha^nspecies
  long nspecies;                       // Number of species
  long nalpha;                         // Number of DNA bases
};


void print_alignment(Alignment &al, std::ostream &s=std::cout);
void save_alignment(Alignment &al, std::string const &fname);

void add_pseudocounts(double epsilon, Counts &data);
void add_pseudocounts(double epsilon, std::vector<double> &data, long N);
void get_counts(Alignment &align, Counts &data);
Counts read_counts(Tree &T, std::string fname, long nalpha=4);
void print_counts(Counts &data, std::ostream &s=std::cout);

#endif
