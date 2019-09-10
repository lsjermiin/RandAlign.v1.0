/*--------------------------------------------------------------------------
 Program name     : RandAlign.cpp
 
 Version          : 1.0
 
 Author           : Lars S Jermiin
 
 Institution      : Australian National University
                    Research School of Biology
                    Acton, ACT 2601, Australia
 
                    University College Dublin
                    School of Biology and Environmental Science
                    Belfield, Dublin 4, Ireland
 
 Date begun       : 12 November, 2018
 
 Date modified    : 10 September, 2019
 
 Copyright        : Copyright © 2019 Lars Sommer Jermiin. All rights reserved.
 
 Responsibility   : The copyright holder takes no legal responsibility for
                    the correctness of results obtained using this program.
 
 Summary          : RandAlign reads a multiple sequence alignment and shuffles
                    the content at the sites in three ways:
                        1. across sites within each sequence [retains entropy
                           of each sequence]
                        2. across sequences within each site [retains entropy
                           of each site]
                        3. across sites as well as sequences
 
                    All sites or variant sites are considered.
 
                    The resulting two alignments are printed to two files.
 
 Input format     : Sequences must be stored in the FASTA format.
 
 Input types      : Sequences can be read a strings of one, two, or three
                    nucleotides, strings of 10- or 14-state genotypes, or
                    as strings of amino acids.
 
 Output           : Sequences will be saved in the FASTA or sequential PHYLIP
                    formats.
 
 Reference        : Jermiin et al. (2019) ...
 
 --------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <iomanip>
#include <fstream>
#include <iostream>

// These variables are declared externally because they are needed in different functions
const unsigned FOUR(4);         // for 4-state alphabet (single nucleotides)
const unsigned SIXTEEN(16);     // for 16-state alphabet (duplet nucleotides)
const unsigned TEN(10);         // for 10-state alphabet (genotype data)
const unsigned FOURTEEN(14);    // for 14-state alphabet (genotype data)
const unsigned TWENTY(20);      // for 20-state alphabet (amino acids)
const unsigned SIXTYFOUR(64);   // for 64-state alphabet (triplet nucleotides)
std::string sites("");                  // string controlling whether a site is used or not
std::vector<std::string> taxon;           // 2D container for sequence names
std::vector<std::vector<unsigned> > alignment; // 2D container for sequence data


// This function translates a string of characters into a vector of integers
std::vector<unsigned> Translator(unsigned datatype, std::string seq) {
    int unit; // integer for singlet, duplet or triplet (codon)
    std::string duplet(""), triplet(""); // strings for dinucleotides and codons
    std::vector<unsigned> seq_data;
    
    switch (datatype) {
        case 1: // Single nucleotides (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 2: // Duplet nucleotides (AA|AC|..|TG|TT)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 2) {
                for (std::string::size_type j = i; j != i + 2; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': duplet.push_back('0'); break;
                        case 'C': duplet.push_back('1'); break;
                        case 'G': duplet.push_back('2'); break;
                        case 'T': duplet.push_back('3'); break;
                        case 'U': duplet.push_back('3'); break;
                        default : duplet.push_back('4'); break;
                    }
                }
                unit = stoi(duplet);
                duplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 3: // Triplet nucleotides (AAA|AAC|...|TTG|TTT)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 000: seq_data.push_back(0);  break; // AAA
                    case 001: seq_data.push_back(1);  break; // AAC
                    case 002: seq_data.push_back(2);  break; // AAG
                    case 003: seq_data.push_back(3);  break; // AAT
                    case 010: seq_data.push_back(4);  break; // ACA
                    case 011: seq_data.push_back(5);  break; // ACC
                    case 012: seq_data.push_back(6);  break; // ACG
                    case 013: seq_data.push_back(7);  break; // ACT
                    case 020: seq_data.push_back(8);  break; // AGA
                    case 021: seq_data.push_back(9);  break; // AGC
                    case 022: seq_data.push_back(10); break; // AGG
                    case 023: seq_data.push_back(11); break; // AGT
                    case 030: seq_data.push_back(12); break; // ATA
                    case 031: seq_data.push_back(13); break; // ATC
                    case 032: seq_data.push_back(14); break; // ATG
                    case 033: seq_data.push_back(15); break; // ATT
                    case 100: seq_data.push_back(16); break; // CAA
                    case 101: seq_data.push_back(17); break; // CAC
                    case 102: seq_data.push_back(18); break; // CAG
                    case 103: seq_data.push_back(19); break; // CAT
                    case 110: seq_data.push_back(20); break; // CCA
                    case 111: seq_data.push_back(21); break; // CCC
                    case 112: seq_data.push_back(22); break; // CCG
                    case 113: seq_data.push_back(23); break; // CCT
                    case 120: seq_data.push_back(24); break; // CGA
                    case 121: seq_data.push_back(25); break; // CGC
                    case 122: seq_data.push_back(26); break; // CGG
                    case 123: seq_data.push_back(27); break; // CGT
                    case 130: seq_data.push_back(28); break; // CTA
                    case 131: seq_data.push_back(29); break; // CTC
                    case 132: seq_data.push_back(30); break; // CTG
                    case 133: seq_data.push_back(31); break; // CTT
                    case 200: seq_data.push_back(32); break; // GAA
                    case 201: seq_data.push_back(33); break; // GAC
                    case 202: seq_data.push_back(34); break; // GAG
                    case 203: seq_data.push_back(35); break; // GAT
                    case 210: seq_data.push_back(36); break; // GCA
                    case 211: seq_data.push_back(37); break; // GCC
                    case 212: seq_data.push_back(38); break; // GCG
                    case 213: seq_data.push_back(39); break; // GCT
                    case 220: seq_data.push_back(40); break; // GGA
                    case 221: seq_data.push_back(41); break; // GGC
                    case 222: seq_data.push_back(42); break; // GGG
                    case 223: seq_data.push_back(43); break; // GGT
                    case 230: seq_data.push_back(44); break; // GTA
                    case 231: seq_data.push_back(45); break; // GTC
                    case 232: seq_data.push_back(46); break; // GTG
                    case 233: seq_data.push_back(47); break; // GTT
                    case 300: seq_data.push_back(48); break; // TAA
                    case 301: seq_data.push_back(49); break; // TAC
                    case 302: seq_data.push_back(50); break; // TAG
                    case 303: seq_data.push_back(51); break; // TAT
                    case 310: seq_data.push_back(52); break; // TCA
                    case 311: seq_data.push_back(53); break; // TCC
                    case 312: seq_data.push_back(54); break; // TCG
                    case 313: seq_data.push_back(55); break; // TCT
                    case 320: seq_data.push_back(56); break; // TGA
                    case 321: seq_data.push_back(57); break; // TGC
                    case 322: seq_data.push_back(58); break; // TGG
                    case 323: seq_data.push_back(59); break; // TGT
                    case 330: seq_data.push_back(60); break; // TTA
                    case 331: seq_data.push_back(61); break; // TTC
                    case 332: seq_data.push_back(62); break; // TTG
                    case 333: seq_data.push_back(63); break; // TTT
                    default:  seq_data.push_back(64); break; // In case of other characters
                }
            }
            break;
        case 4: // 10-state genotype data
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    default : seq_data.push_back(10);break; // In case of other characters
                }
            }
            break;
        case 5: // 14-state genotype data
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    case 'B': seq_data.push_back(10);break;
                    case 'D': seq_data.push_back(11);break;
                    case 'H': seq_data.push_back(12);break;
                    case 'V': seq_data.push_back(13);break;
                    default : seq_data.push_back(14);break; // In case of other characters
                }
            }
            break;
        default: // amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(2); break;
                    case 'E': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'G': seq_data.push_back(5); break;
                    case 'H': seq_data.push_back(6); break;
                    case 'I': seq_data.push_back(7); break;
                    case 'K': seq_data.push_back(8); break;
                    case 'L': seq_data.push_back(9); break;
                    case 'M': seq_data.push_back(10);break;
                    case 'N': seq_data.push_back(11);break;
                    case 'P': seq_data.push_back(12);break;
                    case 'Q': seq_data.push_back(13);break;
                    case 'R': seq_data.push_back(14);break;
                    case 'S': seq_data.push_back(15);break;
                    case 'T': seq_data.push_back(16);break;
                    case 'V': seq_data.push_back(17);break;
                    case 'W': seq_data.push_back(18);break;
                    case 'Y': seq_data.push_back(19);break;
                    default : seq_data.push_back(20);break; // In case of other characters
                }
            }
            break;
    }
    return(seq_data);
}


// This function translates a vector of intergers into a string of characters
std::string Back_translator(unsigned datatype, std::vector<unsigned> seq_data) {
    unsigned number(0);
    std::string str("");
    
    switch (datatype) {
        case 1:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                number = seq_data[i];
                switch (number) {
                    case 0:  str.push_back('A'); break;
                    case 1:  str.push_back('C'); break;
                    case 2:  str.push_back('G'); break;
                    case 3:  str.push_back('T'); break;
                    default: str.push_back('-'); break;
                }
            }
            break;
        case 2:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                number = seq_data[i];
                switch (number) {
                    case 00: str = str + "AA"; break;
                    case 01: str = str + "AC"; break;
                    case 02: str = str + "AG"; break;
                    case 03: str = str + "AT"; break;
                    case 10: str = str + "CA"; break;
                    case 11: str = str + "CC"; break;
                    case 12: str = str + "CG"; break;
                    case 13: str = str + "CT"; break;
                    case 20: str = str + "GA"; break;
                    case 21: str = str + "GC"; break;
                    case 22: str = str + "GG"; break;
                    case 23: str = str + "GT"; break;
                    case 30: str = str + "TA"; break;
                    case 31: str = str + "TC"; break;
                    case 32: str = str + "TG"; break;
                    case 33: str = str + "TT"; break;
                    default: str = str + "--"; break;
                }
            }
            break;
        case 3:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                number = seq_data[i];
                switch (number) {
                    case 000: str = str + "AAA"; break;
                    case 001: str = str + "AAC"; break;
                    case 002: str = str + "AAG"; break;
                    case 003: str = str + "AAT"; break;
                    case 010: str = str + "ACA"; break;
                    case 011: str = str + "ACC"; break;
                    case 012: str = str + "ACG"; break;
                    case 013: str = str + "ACT"; break;
                    case 020: str = str + "AGA"; break;
                    case 021: str = str + "AGC"; break;
                    case 022: str = str + "AGG"; break;
                    case 023: str = str + "AGT"; break;
                    case 030: str = str + "ATA"; break;
                    case 031: str = str + "ATC"; break;
                    case 032: str = str + "ATG"; break;
                    case 033: str = str + "ATT"; break;
                    case 100: str = str + "CAA"; break;
                    case 101: str = str + "CAC"; break;
                    case 102: str = str + "CAG"; break;
                    case 103: str = str + "CAT"; break;
                    case 110: str = str + "CCA"; break;
                    case 111: str = str + "CCC"; break;
                    case 112: str = str + "CCG"; break;
                    case 113: str = str + "CCT"; break;
                    case 120: str = str + "CGA"; break;
                    case 121: str = str + "CGC"; break;
                    case 122: str = str + "CGG"; break;
                    case 123: str = str + "CGT"; break;
                    case 130: str = str + "CTA"; break;
                    case 131: str = str + "CTC"; break;
                    case 132: str = str + "CTG"; break;
                    case 133: str = str + "CTT"; break;
                    case 200: str = str + "GAA"; break;
                    case 201: str = str + "GAC"; break;
                    case 202: str = str + "GAG"; break;
                    case 203: str = str + "GAT"; break;
                    case 210: str = str + "GCA"; break;
                    case 211: str = str + "GCC"; break;
                    case 212: str = str + "GCG"; break;
                    case 213: str = str + "GCT"; break;
                    case 220: str = str + "GGA"; break;
                    case 221: str = str + "GGC"; break;
                    case 222: str = str + "GGG"; break;
                    case 223: str = str + "GGT"; break;
                    case 230: str = str + "GTA"; break;
                    case 231: str = str + "GTC"; break;
                    case 232: str = str + "GTG"; break;
                    case 233: str = str + "GTT"; break;
                    case 300: str = str + "TAA"; break;
                    case 301: str = str + "TAC"; break;
                    case 302: str = str + "TAG"; break;
                    case 303: str = str + "TAT"; break;
                    case 310: str = str + "TCA"; break;
                    case 311: str = str + "TCC"; break;
                    case 312: str = str + "TCG"; break;
                    case 313: str = str + "TCT"; break;
                    case 320: str = str + "TGA"; break;
                    case 321: str = str + "TGC"; break;
                    case 322: str = str + "TGG"; break;
                    case 323: str = str + "TGT"; break;
                    case 330: str = str + "TTA"; break;
                    case 331: str = str + "TTC"; break;
                    case 332: str = str + "TTG"; break;
                    case 333: str = str + "TTT"; break;
                    default:  str = str + "---"; break;
                }
            }
            break;
        case 4:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                number = seq_data[i];
                switch (number) {
                    case 0:  str.push_back('A'); break;
                    case 1:  str.push_back('C'); break;
                    case 2:  str.push_back('G'); break;
                    case 3:  str.push_back('T'); break;
                    case 4:  str.push_back('K'); break;
                    case 5:  str.push_back('M'); break;
                    case 6:  str.push_back('R'); break;
                    case 7:  str.push_back('Y'); break;
                    case 8:  str.push_back('S'); break;
                    case 9:  str.push_back('W'); break;
                    default: str.push_back('-'); break;
                }
            }
            break;
        case 5:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                number = seq_data[i];
                switch (number) {
                    case  0: str.push_back('A'); break;
                    case  1: str.push_back('C'); break;
                    case  2: str.push_back('G'); break;
                    case  3: str.push_back('T'); break;
                    case  4: str.push_back('K'); break;
                    case  5: str.push_back('M'); break;
                    case  6: str.push_back('R'); break;
                    case  7: str.push_back('Y'); break;
                    case  8: str.push_back('S'); break;
                    case  9: str.push_back('W'); break;
                    case 10: str.push_back('B'); break;
                    case 11: str.push_back('D'); break;
                    case 12: str.push_back('H'); break;
                    case 13: str.push_back('V'); break;
                    default: str.push_back('-'); break;
                }
            }
            break;
        default:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                number = seq_data[i];
                switch (number) {
                    case  0: str.push_back('A'); break;
                    case  1: str.push_back('C'); break;
                    case  2: str.push_back('D'); break;
                    case  3: str.push_back('E'); break;
                    case  4: str.push_back('F'); break;
                    case  5: str.push_back('G'); break;
                    case  6: str.push_back('H'); break;
                    case  7: str.push_back('I'); break;
                    case  8: str.push_back('K'); break;
                    case  9: str.push_back('L'); break;
                    case 10: str.push_back('M'); break;
                    case 11: str.push_back('N'); break;
                    case 12: str.push_back('P'); break;
                    case 13: str.push_back('Q'); break;
                    case 14: str.push_back('R'); break;
                    case 15: str.push_back('S'); break;
                    case 16: str.push_back('T'); break;
                    case 17: str.push_back('V'); break;
                    case 18: str.push_back('W'); break;
                    case 19: str.push_back('Y'); break;
                    default: str.push_back('-'); break;
                }
            }
            break;
    }
    return(str);
}

// Function that reads input file and stores data in two 2D containers
unsigned long Read_Input(std::string inname, unsigned datatype){
    unsigned long alignment_length(0);
    unsigned long counter(0);
    std::string seq(""), str(""); // temporary string used to store input
    std::vector<unsigned> sequence;    // temporary vector used to store input
    std::ifstream infile;
    
    infile.open(inname.c_str());
    if (!infile) {
        std::cerr << "\nERROR: input file not found" << std::endl;
        exit(1);
    }
    while (getline(infile, str)) {
        if (!str.empty()) {
            if (str[0] == '>') {
                if (seq.size() > 0) {
                    if (datatype == 2) {
                        if (seq.size() % 2 != 0) {
                            std::cerr << "\nERROR: expected sequence of di-nucleotides" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    if (datatype == 3) {
                        if (seq.size() % 3 != 0) {
                            std::cerr << "\nERROR: expected sequence of codons" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    sequence = Translator(datatype, seq);
                    alignment.push_back(sequence); // stores sequence in vector
                    if (alignment_length == 0)
                        alignment_length = sequence.size();
                    sequence.clear();
                    seq.clear();
                }
                str.erase(str.begin()); // removes first character from name
                taxon.push_back(str); // stores sequence name in vector
                str.clear();
            } else {
                seq += str;
            }
        }
    }
    // Store last sequence in vector
    if (seq.size() > 0) {
        sequence = Translator(datatype, seq);
        alignment.push_back(sequence);
    } else {
        std::cerr << "\nERROR: last sequence empty" << "\n" << std::endl;
        exit(1);
    }
    //Check whether the sequence names are unique
    for (std::vector<std::string>::const_iterator iter1 = taxon.begin(); iter1 != taxon.end(); ++iter1) {
        for (std::vector<std::string>::const_iterator iter2 = iter1 + 1; iter2 != taxon.end(); ++iter2) {
            if (*iter1 == *iter2) {
                std::cerr << "\nERROR: sequence names not unique -- look for " << *iter1 << "\n" << std::endl;
                exit(1);
            }
        }
    }
    // Check whether the sequences have the same length
    for (std::vector<std::vector<unsigned> >::const_iterator iter = alignment.begin()+1; iter != alignment.end(); ++iter) {
        ++counter;
        sequence = *iter;
        if (sequence.size() != alignment_length) {
            std::cerr << "\nERROR: sequences 1 and " << counter << " differ in length!\n" << std::endl;
            exit(1);
        }
    }
    return(alignment_length);
}


// Function used to identify variant sites
void Identify_Variant_Sites(std::string choice_of_sites, unsigned states, unsigned long alignment_length) {
    std::vector<int> column;
    
    for (std::string::size_type j = 0; j != sites.size(); ++j) {
        for (std::string::size_type i = 0; i != taxon.size(); ++i) {
            column.push_back(alignment[i][j]);
        }
        sort(column.begin(),column.end());
        std::vector<int>::iterator end_unique = unique(column.begin(),column.end());
        column.erase(end_unique,column.end());
        if (column.size() == 1 || (column.size() == 2 && column[1] == states)) {
            sites[j] = '0'; // exclude constant sites from the analysis
        }
        column.clear();
    }
}


int main(int argc, char** argv) {
    unsigned long k(0);
    unsigned alphabet(0);
    unsigned dataType(0);
    unsigned long alignment_length(0);
    unsigned long sequence_number(0);
    std::string inName(""), outName(""); // File names
    std::string choice_of_sites, method, nature_of_data, forMat;
    std::vector<unsigned> data, variant;
    std::mt19937_64 generator;
    std::ifstream infile;
    std::ofstream outfile;

    if(argc != 6) {
        std::cerr << "\n\nERROR -- use command: randalign <infile> <a|v> <r|c|b> <1|...|6> <FASTA|PHYLIP>\n" << std::endl;
        std::cerr << "  infile   Fasta-formatted alignment" << std::endl;
        std::cerr << "       a   All sites" << std::endl;
        std::cerr << "       v   Variant sites" << std::endl;
        std::cerr << "       r   Randomize across sites within each sequence" << std::endl;
        std::cerr << "       c   Randomize across sequences within each site" << std::endl;
        std::cerr << "       b   Randomize across sites as well as sequences" << std::endl;
        std::cerr << "       1   Single nucleotides; 4 states (A|C|G|T)" << std::endl;
        std::cerr << "       2   Duplet nucleotides; 16 states (AA|AC|...|TG|TT)" << std::endl;
        std::cerr << "       3   Triplet nucleotides; 64 states (AAA|AAC|...|TTG|TTT)" << std::endl;
        std::cerr << "       4   Genotypes; 10 states (A|C|G|T|K|M|R|Y|S|W)" << std::endl;
        std::cerr << "       5   Genotypes; 14 states (A|C|G|T|K|M|R|Y|S|W|B|D|H|V)" << std::endl;
        std::cerr << "       6   Amino acids; 20 states (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)" << std::endl;
        std::cerr << "   FASTA   File format" << std::endl;
        std::cerr << "  PHYLIP   File format (sequential Phylip)\n" << std::endl;
        exit(1);
    }
    inName = argv[1];
    choice_of_sites = argv[2];
    method = argv[3];
    nature_of_data = argv[4];
    forMat = argv[5];
    // check availability of input file
    infile.open(inName.c_str());
    if (!infile) {
        std::cerr << "\nERROR: file not found...\n" << std::endl;
        exit(1);
    }
    // check choice of sites
    if (toupper(choice_of_sites[0]) != 'A' && toupper(choice_of_sites[0]) != 'V') {
        std::cerr << "\nERROR: incorrect choice of sites: [a|v]\n" << std::endl;
        exit(1);
    }
    // check choice of method
    if (toupper(method[0]) != 'R' && toupper(method[0]) != 'C' && toupper(method[0]) != 'B') {
        std::cerr << "\nERROR: permute within rows (r), columns (c) or both (b)\n" << std::endl;
        exit(1);
    }
    // check choice of data and alphabet
    dataType = stoi(nature_of_data);
    if (dataType < 1 || dataType > 6) {
        std::cerr << "\nERROR: incorrect choice of data: [1|...|6]\n" << std::endl;
        exit(1);
    }
    for (std::string::size_type i = 0; i != forMat.size(); i++) {
        forMat[i] = toupper(forMat[i]);
    }
    // check choice of output file format
    if (forMat != "FASTA" && forMat != "PHYLIP") {
        std::cerr << "\nERROR: format of outfile must be FASTA or PHYLIP\n" << std::endl;
        exit(1);
    }
    for (std::string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
        outName += inName[i];
    }
    if (forMat == "FASTA") {
        outName = outName + "_R" + choice_of_sites[0] + method[0] + ".fst";  // Naming of output file
    } else {
        outName = outName + "_R" + choice_of_sites[0] + method[0] + ".phy";  // Naming of output file
    }
    switch (dataType) {
        case 1: alphabet = FOUR; break;
        case 2: alphabet = SIXTEEN; break;
        case 3: alphabet = SIXTYFOUR; break;
        case 4: alphabet = TEN; break;
        case 5: alphabet = FOURTEEN; break;
        default: alphabet = TWENTY; break;
    }
    alignment_length = Read_Input(inName, dataType);
    // populate guide string sites with 1s
    for (std::string::size_type i = 0; i != alignment_length; ++i) {
        sites.push_back('1');
    }
    // change 1s to 0s in guide string sites if they are constant
    if (toupper(choice_of_sites[0]) == 'V') {
        Identify_Variant_Sites(choice_of_sites, alphabet, alignment_length);
    }
//    std::cout << sites << std::endl;
    // obtain a seed for the random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // priming the Mersene Twister 19937 generator (64 bit)
    generator = std::mt19937_64(seed);
    outfile.open(outName.c_str());
    switch (method[0]) {
        case 'r':
            // Permutation across sites within each sequence
            for (std::vector<std::vector<unsigned> >::size_type i = 0; i != alignment.size(); ++i) {
                data = alignment[i];
                for (std::vector<unsigned>::size_type j = 0; j != data.size(); ++j) {
                    if (sites[j] == '1') {
                        variant.push_back(data[j]);
                    }
                }
                std::shuffle(variant.begin(), variant.end(), generator);
                k = 0;
                for (std::vector<unsigned>::size_type j = 0; j != data.size(); ++j) {
                    if (sites[j] == '1') {
                        data[j] = variant[k];
                        ++k;
                    } 
                }
                alignment[i] = data;
                variant.clear();
            }
            break;
        case 'c':
            // Permutation across sequences with each site
            for (unsigned i = 0; i != alignment_length; ++i) {
                for (std::vector<std::vector<unsigned> >::size_type j = 0; j != alignment.size(); ++j) {
                    data.push_back(alignment[j][i]);
                }
                std::shuffle(data.begin(), data.end(), generator);
                for (std::vector<std::vector<unsigned> >::size_type j = 0; j != alignment.size(); ++j) {
                    alignment[j][i] = data[j];
                }
                data.clear();
            }
            break;
        default:
            // Permutation across sites as well as sequences
            for (std::vector<std::vector<unsigned> >::size_type i = 0; i != alignment.size(); ++i) {
                data = alignment[i];
                for (std::vector<unsigned>::size_type j = 0; j != data.size(); ++j) {
                    if (sites[j] == '1') {
                        variant.push_back(data[j]);
                    }
                }
                data.clear();
            }
//            for (std::vector<unsigned>::size_type j = 0; j != variant.size(); ++j) {
//                std::cout << variant[j];
//            }
//            std::cout << std::endl;
            std::shuffle(variant.begin(), variant.end(), generator);
//            for (std::vector<unsigned>::size_type j = 0; j != variant.size(); ++j) {
//                std::cout << variant[j];
//            }
//            std::cout << std::endl;
            k = 0;
            for (std::vector<std::vector<unsigned> >::size_type i = 0; i != alignment.size(); ++i) {
                data = alignment[i];
                for (std::vector<unsigned>::size_type j = 0; j != data.size(); ++j) {
                    if (sites[j] == '1') {
                        data[j] = variant[i+k];
                        ++k;
                    }
                }
                --k;
                alignment[i] = data;
                data.clear();
            }
            break;
    }
    // Printing output to file
    if (forMat == "FASTA") {
        for (std::vector<std::vector<unsigned> >::size_type i = 0; i != alignment.size(); i++) {
            outfile << ">" << taxon[i] << std::endl;
            data = alignment[i];
            Back_translator(dataType, data);
            outfile << Back_translator(dataType, data) << std::endl;
            data.clear();
        }
    } else {
        outfile << sequence_number + 1 << "  " << alignment_length + 1 << std::endl;
        for (std::vector<std::string>::size_type i = 0; i != alignment.size(); i++) {
            outfile << std::left << std::setw(10) << taxon[i];
            Back_translator(dataType, data);
            outfile << Back_translator(dataType, data) << std::endl;
            data.clear();
        }
    }
    outfile << std::endl;
    outfile.close();
    std::cout << "Permutation completed" << std::endl;
    return 0;
}
