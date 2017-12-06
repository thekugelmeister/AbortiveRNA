#include "rbsMatch.h"

RBS::RBS(string line)
{
    unsigned pos1 = 0;
    unsigned pos2 = 0;

    string buffer;

    //get name;
    pos1 = line.find("\t", 0) + 1;
    pos2 = line.find("\t", pos1);
    name = line.substr(pos1, (pos2 - pos1));

    //get left position
    pos1 = line.find("\t", pos2) + 1;
    pos2 = line.find("\t", pos1);
    buffer = line.substr(pos1, (pos2 - pos1));
    leftPos = atoi(buffer.c_str());

    //get right position
    pos1 = line.find("\t", pos2) + 1;
    pos2 = line.find("\t", pos1);
    buffer = line.substr(pos1, (pos2 - pos1));
    rightPos = atoi(buffer.c_str());

    //get strand
    pos1 = line.find("\t", pos2) + 1;
    pos2 = line.find("\t", pos1);
    strand = line.substr(pos1, (pos2 - pos1));

    //get center position
    pos1 = line.find("\t", pos2) + 1;
    pos2 = line.find("\t", pos1);
    buffer = line.substr(pos1, (pos2 - pos1));
    centerPos = atoi(buffer.c_str());

    //get sequence
    pos1 = line.find("\t", pos2) + 1;
    pos2 = line.find("\t", pos1);
    sequence = line.substr(pos1, (pos2 - pos1));
}

RBSList::RBSList(char *fName)
{
    ifstream rFile(fName);
    string buffer;

    while (getline(rFile, buffer)) {
        RBS *r = new RBS(buffer);
        RBSs.push_back(r);
    }
}

void RBSList::findRBSMatches(int l, int r, int len, string strand, vector <RBS *> &list)
{
    for (unsigned i = 0; i < RBSs.size(); i++) {
        if (RBSs[i]->containsMatch(l, r, len, strand)) {
            list.push_back(RBSs[i]);
        }
    }
}

bool RBS::containsMatch(int l, int r, int len, string strand)
{
    if ((l >= (leftPos - (len - 1))) && (r <= (rightPos + (len - 1))) && sameStrand(strand)) {
        return true;
    } else {
        return false;
    }
}

bool RBS::sameStrand(string s)
{
    if ((strand == "forward") && (s == "forward")) {
        return true;
    } else if ((strand == "reverse") && (s == "reverse")) {
        return true;
    } else {
        return false;
    }
}
