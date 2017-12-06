#include "termParse.h"

int main(int argc, char *argv[])
{
    if (argc != 4) {
        cerr << "ERROR: Usage is ./termParse [abortive match file identifier] "
             << "[terminator file identifier] [strand]" << endl;
        exit(1);
    }
    
    Terminators t;
    t.getAbortiveMatches(argv[1]);
    t.setTermFile(argv[2]);
    string outputFile = argv[1];
    string strand = argv[3];
    outputFile = outputFile + ".trmout." + strand;
    cerr << "Output in: " << outputFile << endl;
    t.setOutputFile((char *) outputFile.c_str());
    while(t.getNextTerminator(argv[3])) {
        t.searchForAbortiveOverlaps();
    }
    return 0;
}
