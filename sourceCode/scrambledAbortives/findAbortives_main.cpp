#include "findAbortives.h"

int main(int argc, char *argv[])
{
    srand(time(NULL));
    if (argc != 8) {
        cerr << "ERROR: Usage is ./findAbortives [genome file] [transcript file] [output file] [sequence or wobble] [numPermutations] [minLen] [maxLen]" << endl;
        exit(1);
    }
    char *gFile = argv[1];
    char *tFile = argv[2];
    string oFile = argv[3];
    string matchType = argv[4];
    int numPermutations = atoi(argv[5]);
    
    int minLen = atoi(argv[6]);
    int maxLen = atoi(argv[7]);
    int maxMismatches = 0;

    if (matchType != "sequence" && matchType != "wobble") {
        cerr << "ERROR: unknown matching type" << endl;
        exit(1);
    }

    TranscriptList t(gFile, tFile, numPermutations);
    t.setTranscriptSequences();

    //perform the analysis on abortives of different lengths
    for (int len = minLen; len <= maxLen; len++) {
        //generate an appropriate file name variant based on the abortive length
        stringstream convert;
        convert << len;
        string currFile = oFile + convert.str();
        t.setOutputFile((char *) currFile.c_str());

        t.setAbortiveLengths(len);

        if (matchType == "sequence") {
            t.searchForSequenceMatches(maxMismatches);
        } else if (matchType == "wobble") {
            t.searchForWobbleMatches(maxMismatches);
        }

        // t.printMatchStatistics();
        t.printAllMatches();
        t.reset();
    }
    return 0;
}
