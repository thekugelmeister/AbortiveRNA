#include "findAbortives.h"

int main(int argc, char *argv[])
{
    if (argc != 5) {
        cerr << "ERROR: Usage is ./findAbortives [genome file] [transcript file] [output file] [maxLen]" << endl;
        exit(1);
    }
    char *gFile = argv[1];
    char *tFile = argv[2];
    string oFile = argv[3];

    int maxLen = atoi(argv[4]);
    int numPermutations = 0;
    int maxMismatches = 0;

    TranscriptList t(gFile, tFile, numPermutations);
    t.setTranscriptSequences();

    //generate an appropriate file name variant based on the abortive length
    stringstream convert;
    convert << maxLen;
    string currFile = oFile + convert.str();
    t.setOutputFile((char *) currFile.c_str());

    t.setAbortiveLengths(maxLen);

    t.printAllSequences();
    t.reset();
    return 0;
}
