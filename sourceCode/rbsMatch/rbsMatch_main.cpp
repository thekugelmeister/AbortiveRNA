#include "rbsMatch.h"

int main(int argc, char *argv[])
{
    if (argc != 3) {
        cerr << "ERROR: Usage is ./RBSMatch [RBS file] [abortive match file]" << endl;
        exit(1);
    }

    char *RBSFile = argv[1];
    char *abortiveFile = argv[2];
    string outputFile = argv[2];
    outputFile = outputFile + ".RBSMatches";

    //create the gene list from the file
    RBSList r(RBSFile);

    //create the list of abortive matches from the file
    ifstream aFile(abortiveFile);
    string buffer;
    int abortiveLength;
    int numMatches;
    //get the length of the matches
    aFile >> abortiveLength;
    //discard the rest of that line and two more lines
    for (int i = 0; i < 3; i++) {
        getline(aFile, buffer);
    }
    //get rid of the first three words on the next line
    for (int i = 0; i < 3; i++) {
        aFile >> buffer;
    }
    aFile >> numMatches;
    
    for (int i = 0; i < 7; i++) {
        getline(aFile, buffer);
    }
    //take in the relevant information and populate the geneMatches list
    AbortiveMatch matches[numMatches];
    for (int i = 0; i < numMatches; i++) {
        aFile >> matches[i].operon
              >> matches[i].transcriptName
              >> matches[i].leftBound
              >> matches[i].rightBound
              >> matches[i].strand
              >> matches[i].startPosition
              >> matches[i].bindSite
              >> matches[i].abortive
              >> matches[i].sDev;
        matches[i].matchLength = abortiveLength;

        r.findRBSMatches(matches[i].startPosition,
                          (matches[i].startPosition + abortiveLength - 1),
                          matches[i].matchLength, matches[i].strand,  matches[i].RBSMatches);
    }

    //open output file and print those with >0 rbs matches
    ofstream oFile(outputFile.c_str());
    int count = 0;
    for (int i = 0; i < numMatches; i++) {
        if (matches[i].RBSMatches.size() != 0) {
            count += 1;
            oFile << "Match:\t\t"
                  << matches[i].startPosition << "\t"
                  << matches[i].startPosition + matches[i].matchLength - 1 << "\t"
                  << matches[i].sDev << "\t"
                  << matches[i].bindSite << "\t" << matches[i].abortive << endl
                  << "Transcript:\t"
                  << matches[i].operon << "\t"
                  << matches[i].transcriptName << "\t"
                  << matches[i].leftBound << "\t"
                  << matches[i].rightBound << "\t"
                  << matches[i].strand << endl;
            oFile << "--------------------------------------------------" << endl;
            for (unsigned j = 0; j < matches[i].RBSMatches.size(); j++) {
                oFile << "RBS:\t\t"
                      << matches[i].RBSMatches[j]->name << "\t"
                      << matches[i].RBSMatches[j]->strand << "\t"
                      << matches[i].RBSMatches[j]->leftPos << "\t"
                      << matches[i].RBSMatches[j]->rightPos << endl
                      << matches[i].RBSMatches[j]->sequence << "\t"
                      << min((matches[i].startPosition + matches[i].matchLength - 1), (matches[i].RBSMatches[j]->rightPos)) - max((matches[i].startPosition), (matches[i].RBSMatches[j]->leftPos)) + 1 << endl;
            }
            oFile << endl;
        }
    }
    oFile << "Total:\t" << count << endl << endl;
    
    return 0;
}
