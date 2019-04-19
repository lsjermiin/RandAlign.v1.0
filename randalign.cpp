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
 
 Date modified    : 5 March, 2019
 
 Copyright        : Copyright © 2019 Lars Sommer Jermiin. All rights reserved.
 
 Responsibility   : The copyright holder takes no legal responsibility for
                    the correctness of results obtained using this program.
 
 Summary          : RandAlign reads a multiple sequence alignment, shuffles
                    the content in three ways:
                        1. across sites within each sequence [retains entropy
                           of each sequence]
                        2. across sequences within each site [retains entropy
                           of each site]
                        3. across sites as well as sequences
                    The resulting two alignments are printed to two files.
 
 Input            : Sequences must be stored in the FASTA format.
 
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

// The following are declared here because they are needed in different functions

std::ifstream infile;
unsigned long alignment_length(0);
unsigned long sequence_number(0);
std::vector<std::string> taxon;               // 2D container for sequence names
std::vector<std::string> alignment;            // 2D container for sequence data


// Function that reads input file and stores data in two 2D containers

int Read_Input(std::string inName){
    std::string str1("");        // temporary string used to store input
    std::string str2("");        // temporary string used to store input
    std::string temp;
    bool success(true);
    
    infile.open(inName.c_str());
    while (getline(infile, str1)) {
        if (!str1.empty()) {
            if (str1[0] == '>') {
                temp = str1;
                str1.clear();
                for (std::string::size_type i = 1; i != temp.size(); i++) {
                    str1.push_back(temp[i]);
                }
                // Store sequence name in vector
                if (str1.size() > 0) {
                    taxon.push_back(str1);
                    ++sequence_number; // Needed to control output
                    str1 = "";
                }
                // Store sequence in vector
                if (str2.size() > 0){
                    alignment.push_back(str2);
                    str2 = "";
                }
            } else {
                for (std::string::size_type i = 0; i != str1.size(); ++i) {
                    if (str1[i] == '-' || isalpha(str1[i])) {
                        str2 += toupper(str1[i]);
                    }
                }
            }
        }
    }
    infile.close();
    // Store last sequence in vector
    if (str2.size() > 0) {
        alignment.push_back(str2);
        alignment_length = str2.size();  // Needed to control output
    } else {
        std::cout << "Last sequence empty -- check input file!\n" << std::endl;
        success = false;
    }
    //Check whether the sequence names are unique
    std::cout << std::endl;
    for (std::vector<std::string>::const_iterator iter1 = taxon.begin(); iter1 != taxon.end(); ++iter1) {
        for (std::vector<std::string>::const_iterator iter2 = iter1 + 1; iter2 != taxon.end(); ++iter2) {
            if (*iter1 == *iter2) {
                std::cout << "Sequences have the same name -- look for " << *iter1 << "\n" << std::endl;
                success = false;
            }
        }
    }
    // Check whether the sequences have the same length
    unsigned i(0);
    for (std::vector<std::string>::const_iterator iter1 = alignment.begin(); iter1 != alignment.begin() + 1; ++iter1) {
        str1 = *iter1;
        ++i;
        for (std::vector<std::string>::const_iterator iter2 = iter1 + 1; iter2 != alignment.end(); ++iter2) {
            str2 = *iter2;
            ++i;
            if (str1.size() != str2.size()) {
                std::cout << "Sequences 1 and " << i << " differ in length!\n" << std::endl;
                success = false;
            }
        }
    }
    return success;
}

int main(int argc, char** argv) {
    bool  success(false);
    std::string data("");
    std::ofstream outfile;
    std::string inName, outName, forMat;
    std::string method;
    std::mt19937_64 generator;

    if(argc != 4) {
        printf("\n\nERROR -- use command: randalign <infile> <r|c|b> <FASTA|PHYLIP>\n\n");
        printf("  r       Randomize across sites within each sequence\n");
        printf("  c       Randomize across sequences within each site\n");
        printf("  b       Randomize across sites as well as sequences\n\n");
        printf("  FASTA   File format\n");
        printf("  PHYLIP  File format (sequential Phylip)\n\n");
        exit(1);
    }
    inName = argv[1];
    method = argv[2];
    forMat = argv[3];
    if (method != "r" && method != "c" && method != "b") {
        printf("\n\nERROR: permute within rows (r), columns (c) or both (b)\n\n");
        exit(1);
    }
    for (std::string::size_type i = 0; i != forMat.size(); i++) {
        forMat[i] = toupper(forMat[i]);
    }
    if (forMat != "FASTA" && forMat != "PHYLIP") {
        printf("\n\nERROR: format of outfile must be FASTA or PHYLIP\n\n");
        exit(1);
    }
    infile.open(inName.c_str());
    if (infile) {
        infile.close();
        success = true;
    }
    else std::cout << "\nWARNING -- file not found...\n\n\n";
    if (success == true) {
        outName = "";
        for (std::string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
            outName += inName[i];
        }
        if (forMat == "FASTA") {
            outName = outName + "_R" + method[0] + ".fst";  // Permutation across sites within each sequence
        } else {
            outName = outName + "_R" + method[0] + ".phy";  // Permutation across sites within each sequence
        }
        if (Read_Input(inName)) {
            outfile.open(outName.c_str());
            // obtain a seed for the random number generator
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            // priming the Mersene Twister 19937 generator (64 bit)
            generator = std::mt19937_64(seed);
            switch (method[0]) {
                case 'r':
                    // Permutation across sites within each sequence
                    for (std::vector<std::string>::size_type i = 0; i != alignment.size(); ++i) {
                        data = alignment[i];
                        std::shuffle(data.begin(), data.end(), generator);
                        alignment[i] = data;
                    }
                    break;
                case 'c':
                    // Permutation across sequences with each site
                    for (std::string::size_type i = 0; i != alignment_length; ++i) {
                        for (std::string::size_type j = 0; j != sequence_number; ++j) {
                            data.push_back(alignment[j][i]);
                        }
                        std::shuffle(data.begin(), data.end(), generator);
                        for (std::string::size_type j = 0; j != sequence_number; ++j) {
                            alignment[j][i] = data[j];
                        }
                        data.clear();
                    }
                    break;
                default:
                    // Permutation across sites as well as sequences
                    for (std::vector<std::string>::const_iterator iter = alignment.begin(); iter != alignment.end(); ++iter) {
                        data += *iter;
                    }
                    std::shuffle(data.begin(), data.end(), generator);
                    unsigned long k(0);
                    for (std::string::size_type i = 0; i != sequence_number; ++i) {
                        for (std::string::size_type j = 0; j != alignment_length; ++j) {
                            alignment[i][j] = data[k];
                            ++k;
                        }
                    }
                    break;
            }
            // Printing output to file
            if (forMat == "FASTA") {
                for (std::vector<std::string>::size_type i = 0; i != alignment.size(); i++) {
                    outfile << ">" << taxon[i] << std::endl;
                    outfile << alignment[i] << std::endl;

                }
            } else {
                outfile << sequence_number + 1 << "  " << alignment_length + 1 << std::endl;
                for (std::vector<std::string>::size_type i = 0; i != alignment.size(); i++) {
                    outfile << std::left << std::setw(10) << taxon[i];
                    outfile << alignment[i] << std::endl;
                }
            }
            outfile << std::endl;
            outfile.close();
        }
    }
    std::cout << "Permutation completed" << std::endl;
    return 0;
}
