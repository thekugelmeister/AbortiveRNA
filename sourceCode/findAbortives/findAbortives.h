/* findAbortives.h
   given a gene sequence, find abortive fragments of a given length
 */

#ifndef BABAT_H
#define BABAT_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include "time.h"
#include <sstream>
#include "jstats.h"
#include <vector>

using namespace std;

//GENOME SEQUENCE CLASS:
//provides retrieval, storage, and accessibility for the genome
class GenomeSequence {
public:
    string getTranscript(int begin, int end);
    void populate(char *fileName);
private:
    string wholeSequence;
};

//TRANSCRIPT UNIT CLASS:
//provides storage and operations for a given transcript unit
class TranscriptUnit {
public:
    /*** Public TU Functions ***/
    TranscriptUnit(int numPermutations);  //constructor initializes match counters
    
    //sets the length of the abortive, then uses that information to retrieve
    //the abortive sequence from the genome
    void setAbortiveLength(int len); 

    //searches for exact sequence matches in the original sequence
    void sequenceMatchSearch(int allowedMismatches);

    //counts matches found in random permutations of the original sequence
    //for statistical testing; populates the permuted count array
    void permutedSequenceSearch(int allowedMismatches, int numPermutations);

    //searches for wobble matches in the original sequence
    void wobbleMatchSearch(int allowedMismatches);

    //same as sequence counterpart but with wobble base pairing
    void permutedWobbleSearch(int allowedMismatches, int numPermutations);

    //calculates the number of standard deviations between the mean of the
    //permuted searches and the original search match counts
    void calculateSignificance(int numPermutations);

    /*** Public TU Attributes ***/
    //Information obtained from the TU source file
    string operon;
    string tu;
    string promoter;
    int tss;  //only use a transcript if this is not ""
    string strand;
    string fg;
    int fglCoord;
    int fgrCoord;
    string lg;
    int lglCoord;
    int lgrCoord;
    string trm;
    int trmlCoord;
    int trmrCoord;
    int utrlCoord;  //should suffice for boundaries of transcript
    int utrrCoord;  //^
    int flCoord;
    int frCoord;
    string fSeq;  //for assertion
    int tlCoord;
    int trCoord;
    string tSeq;  //for assertion

    //Information to be obtained or generated elsewhere

    //the exact sequence of the forward strand of this TU (from genome file)
    //5' to 3'
    string forwardSequence;
    //the exact sequence of the reverse strand of this TU (from forwardSequence)
    //5' to 3'
    string reverseSequence;
    //the exact sequence of the abortive (from the appropriate strand)
    //3' to 5'
    string abortiveSequence;

    //storage for permuted versions of the original sequence
    string permutedForward;
    string permutedReverse;
    string permutedAbortive;

    int abortiveLength;

    //the number of matches for the current transcript's orig. sequence
    int matchCount;
    vector <int> matchStartSites;
    vector <string> matchSequences;

    //stores the counts calculated for each permuted search
    int *permutedMatchCounts;

    //stores the number of standard deviations calculated following all
    //searches for this TU
    double sDev;
    
private:
    void permuteTranscriptSequence();

    //both functions tell if two given sequences match with or without wobble
    //matching respectively
    bool validSequenceMatch(string abortive, string target, int allowedMismatches);
    bool validWobbleMatch(string abortive, string target, int allowedMismatches);
};

//TRANSCRIPT LIST CLASS
//provides storage and functions associated with a list of transcription units,
//along with storing the entire genome sequence
class TranscriptList {
public:
    //constructor populates both the genome sequence and the transcript list
    TranscriptList(char *gFile, char *tFile, int nPerm);

    //takes the sequences for each transcript from the parsed genome file and
    //stores them for easy access later
    void setTranscriptSequences();

    //sets the length for the next search and stores the appropriate sequences
    //for each TU
    void setAbortiveLengths(int len);

    void setOutputFile(char *fName);

    //resets the TU list for the next search
    void reset();

    //sorts the TUs by their sDev score (descending)
    void sortBySDev();

    //performs the appropriate searches for all TUs with the given parameters
    //already set (e.g. abortive length)
    void searchForSequenceMatches(int allowedMismatches);
    void searchForWobbleMatches(int allowedMismatches);

    void printMatchStatistics();
    void printAllMatches();

    GenomeSequence geneSeq;

   
    
    
private:
    //vector storing the list of TUs
    vector <TranscriptUnit> tUnits;

    //gets one transcript from the file; returns false if there are no
    //transcripts left in the file, else returns true
    bool getNextTranscript();

    ifstream transcriptFile;
    ofstream outputFile;

    int totalMatchCount;
    int transcriptsWithMatches;
    int numberOfTranscripts;

    int abortiveLength;

    int numPermutations;
};

#endif
