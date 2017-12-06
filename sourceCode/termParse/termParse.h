#include <iostream>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;


/*******Abortive Match Struct********
 * stores information about an individual match generated from the abortive
 * match program
 */
struct AbortiveMatch
{
    string operon;
    string transcriptName;
    int leftBound;
    int rightBound;
    string strand;
    int startPosition;
    string bindSite; //5' - 3'
    string abortive; //3' - 5'
    double sDev;
};


/*******Terminator Class********
 * stores all information that can possibly be stored for an individual
 * terminator from the files Nik gave me
 */

class Terminators
{
public:
    void getAbortiveMatches(char *fName);
    void setTermFile(char *fName);
    void setOutputFile(char *fName);
    bool getNextTerminator(string currentStrand);
    void searchForAbortiveOverlaps();
    
    string strand;
    int fileIndex;
    int startPosition;
    string upstreamLeg;
    string downstreamLeg;
    string bulb;
    string tail;
    int upstreamLegLen;
    int downstreamLegLen;
    int stemLen;
    int bulbLen;
    int numMismatch;
    int numGap;
    double deltaG;
    string geneName;
    int geneStartPos;
    int geneStopPos;
    int stopToStartDist; //distance from gene stop to structure start
    int stopToStopDist;  //distance from gene stop to structure stop
    int stopToMidDist;   //distance from gene stop to structure middle
    int seqStartPos;        //Not sure what these are for
    string sequence;        //^^
    string description;     //Will eventually need to be broken down further?
private:
    ifstream terminatorFile;
    ofstream outputFile;

    AbortiveMatch * matchList;
    int abortiveLength;
    int numAbortiveMatches;
};
