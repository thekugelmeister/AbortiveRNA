#ifndef RBS_H
#define RBS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

using namespace std;

class RBS {
public:
    RBS(string line);
    bool containsMatch(int l, int r, int len, string strand);  //given the boundaries of a match, tells if
                                                               //that RBS overlaps that match
    bool sameStrand(string strand);  //equates the '+/-' direction identifier with the 'forward/reverse' identifier

    string name;
    int leftPos;
    int rightPos;
    string strand;
    int centerPos;
    string sequence;
};

class RBSList {
public:
    RBSList(char *fName);
    void findRBSMatches(int l, int r, int len, string strand, vector <RBS *> &list);  //populates a list of RBS's
                                                                                      //overlaping the given match position
    vector <RBS *> RBSs;
};

struct AbortiveMatch {
    string operon;
    string transcriptName;
    int leftBound;
    int rightBound;
    string strand;
    int startPosition;
    string bindSite; //5' - 3'
    string abortive; //3' - 5'
    double sDev;
    int matchLength;
    
    vector <RBS *> RBSMatches;
};


#endif
