/*
 *  alignment.cpp
 *
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License.
 *
 */

#include <iostream>
#include <fstream>

#include "alignment.h"
#include "state.h"
#include "seqUtil.h"
#include "read_fasta.h"

#include <stdexcept>


// Cleans a species name. Removes initial '>' and the '\n'.
std::string clean_species_name(std::string str) {
  unsigned int i, i0, i1;
  i0 = 0;
  i1 = str.size()-1;
  for (i=0; i < str.size(); i++) {
    if (str[i] == '>')          i0 = i+1;
    if (str[i] == '\n')         i1 = i-1;
  }
  if (i0 > i1) {
    throw std::invalid_argument("Unmatched species name");

  }
  return str.substr(i0, i1-i0 + 1);
}


// Matches the species in the tree with the species in saved in s (gives an error). This is to ensure that the respective sequences in the fasta file
// match the order in the tree: produces a mapping d translating species in s to
// species in the tree: s[i] = T.name[d[i]].
// returns true if successful, false otherwise

bool match_species(Tree &T, std::string *s, int ns, std::vector<long> &d) {
  long i, j;
  bool found;
  std::string sname;
  if (ns != T.nleaves) {
    std::cout << "Number of leaves on tree does not match number of species on data." << std::endl;
    return false;
  }

  // build the map.
  d.resize(T.nleaves);
  for (i=0; i < T.nleaves; i++){
    sname = clean_species_name(s[i]);
    found = false;
    for (j=0; j < T.nleaves; j++) {
      if (sname == T.names[j]) {
        d[i] = j;
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << "There's a species in the fasta which is not a leaf in the tree." << std::endl;
      return false;
    }
  }

  // check (surjective?)
  for(i=0; i < T.nleaves; i++) {
    found = false;
    for (j=0; j < T.nleaves; j++) {
      if (d[j] == i) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << "There's a species in the tree which is not present in the fasta file." << std::endl;
      return false;
    }
  }

  return true;
}


// Reads the counts from a fasta file using the routines in read_fasta.cpp.
// Then formats the data as we want and stores it into the Counts data structure.
// Uses the order of the leaves on Tree for encoding the states.

Counts read_counts(Tree &T, std::string fname, long nalpha) {
  int i;
  unsigned int j;
  long id;
  State sta, sta2;
  std::string str;
  std::vector<long> d;
  Counts data;

  // The number of letters is fixed to 4 in read_fasta.
  if (nalpha != 4) {
    std::cout << "Reading counts only implemented for 4 letters." << std::endl;
  }

  // reads in the alignment. The data is in _orbits and _couts.
  read_alignment(fname);

  // Now we fill in the Counts structure.
  data.nalpha = 4;
  data.nspecies = g_numSpecies;
  data.nstates = 1;

  // calculate the power:
  for (i=0; i < data.nspecies; i++) {
    data.nstates = data.nstates*data.nalpha;
  }

  // Creates a state of a given dimension
  create_state(sta, data.nspecies, data.nalpha);
  create_state(sta2, data.nspecies, data.nalpha);

  // Matches the species in the fasta with the ones in the tree.
  if (!match_species(T, g_nameSpecies, g_numSpecies, d)) {
    throw std::length_error( "Could not match species in tree with species in fasta." );
  }

  // Stores the counts in the Counts structure.
  data.c.resize(data.nstates);
  data.N = 0;
  for (j=0; j < _orbitals.size(); j++) {
    str = transform_adn_chain_val_to_string(_orbitals[j]);
    string2state(str, sta);
    permute_state(d, sta, sta2);
    id = state2index(sta2);
    data.c[id] = (double) _counts[j];
    data.N = data.N + data.c[id];
  }

  // Stores the species names
  data.species.resize(data.nspecies);
  for (i=0; i < T.nleaves; i++) {
    data.species[i] = T.names[i];
  }

  add_pseudocounts(0.01, data);
}


// Prints out the data in the Counts
void print_counts(Counts &data, std::ostream &s) {
  unsigned int i;
  State sta;
  std::string str;
  create_state(sta, data.nspecies, data.nalpha);

  s << "species:  ";
  for (i=0; i < data.species.size(); i++) {
    s << data.species[i] << "  ";
  }
  s << std::endl;
  s << "nalpha:   " << data.nalpha << std::endl;
  s << "nspecies: " << data.nspecies << std::endl;
  s << "nstates:  " << data.nstates << std::endl;
  s << "total:    " << data.N << std::endl;
  s << "counts:   " << std::endl;
  for (i=0; i < data.c.size(); i++) {
    if (data.c[i] > 0) {
      index2state(i, sta);
      state2string(sta, str);
      s << "  " << str << ": " << data.c[i] << std::endl;
    }
  }
  s << std::endl;
}



// Adds epsilon to all data counts and rescales to keep the same sum.
void add_pseudocounts(double epsilon, Counts &data) {
  long i;

  //  double epsilon = 0.01;
  double factor = data.N/(data.N + data.nstates*epsilon);
  for(i=0; i < data.nstates; i++) {
    data.c[i] = (data.c[i] + epsilon)*factor;
  }
}

void add_pseudocounts(double epsilon, std::vector<double> &data, long N) {
	unsigned long i;
	//  double epsilon = 0.01;
	double factor = (double)N/((double)N + (double)data.size()*epsilon);
	for(i=0; i < data.size(); i++) {
		data[i] = (data[i] + epsilon)*factor;
	}
}

// Gets the counts from the input alignment.
void get_counts(Alignment &align, Counts &data) {
  unsigned int len;
  long i, j, a;
  State st;

  len = align.seq[0].size();
  for(i=0; i < align.nspecies; i++) {
    if (len > align.seq[i].size()) {
      len = align.seq[i].size();
    }
  }

  data.nspecies = align.nspecies;
  data.nalpha = align.nalpha;
  data.N = (double) len;

  data.nstates = 1;
  data.species.resize(data.nspecies);
  for (i=0; i < align.nspecies; i++) {
    data.species[i] = align.name[i];
    data.nstates = data.nstates * data.nalpha;
  }

  // assigns memory for the counts
  data.c.resize(data.nstates);
  for(a=0; a < data.nstates; a++) {
    data.c[a] = 0;
  }

  //ounting ...
  create_state(st, align.nspecies, align.nalpha);
  for(i=0; i < len; i++) {
    for(j=0; j < align.nspecies; j++) {
      st.s[j] = align.seq[j][i];
    }
    a = state2index(st);
    data.c[a] = data.c[a] + 1;
  }
}



// Prints an alignment in the fasta format.
void print_alignment(Alignment &al, std::ostream &s) {
  long i, j;
  int a;
  int width = 60;

  for (i=0; i < al.nspecies; i++) {
    s << ">" << al.name[i];
    for (j=0; j < al.len; j++) {
      a = al.seq[i][j];
      if (j % width == 0) {
	s << std::endl;
      }
      if (a >= 0 && a < al.nalpha) {
        s << al.alpha[a];
      } else {
	s << "X";
      }
    }
    s << std::endl;
  }
}

// Saves an alignment into a file in fasta format.
void save_alignment(Alignment &al, std::string const &fname) {
  std::ofstream ffile;
  ffile.open(fname.c_str(), std::ios::out);
  if (!ffile.is_open()) {
    std::cout << "Could not open file: " << fname << std::endl;
  }
  print_alignment(al, ffile);
  ffile.close();
}
