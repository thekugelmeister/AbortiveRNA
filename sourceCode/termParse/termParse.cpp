#include "termParse.h"

string removeParens(string s)
{
    s.erase(remove(s.begin(), s.end(), '('), s.end());
    s.erase(remove(s.begin(), s.end(), ')'), s.end());
    return s;
}

void Terminators::getAbortiveMatches(char *fName)
{
    ifstream matchFile(fName);
    string buffer;
    matchFile >> abortiveLength; //the first thing in the file is the length of the abortives
    //throw away unnecessary information
    for (int i = 0; i < 14; i++) {
        matchFile >> buffer;
    }
    
    matchFile >> numAbortiveMatches; //get the number of matches

    matchList = new AbortiveMatch[numAbortiveMatches];
    for (int i = 0; i < numAbortiveMatches; i++) {
        matchFile >> matchList[i].operon;
        matchFile >> matchList[i].transcriptName;
        matchFile >> matchList[i].leftBound;
        matchFile >> matchList[i].rightBound;
        matchFile >> matchList[i].strand;
        matchFile >> matchList[i].startPosition;
        matchFile >> matchList[i].bindSite;
        matchFile >> matchList[i].abortive;
        matchFile >> matchList[i].sDev;
    }
}

void Terminators::setTermFile(char *fName)
{
    terminatorFile.open(fName);
    //discard first line of file
    string buffer;
    getline(terminatorFile, buffer);
}

void Terminators::setOutputFile(char *fName)
{
    outputFile.open(fName);
}

bool Terminators::getNextTerminator(string currentStrand)
{
    strand = currentStrand;
    string buffer;
    if (!getline(terminatorFile, buffer)) {
        cerr << "DONE" << endl;
        return false;
    }
    //1. Have it interpret the abbreviations
    //2. Have it locate things based on the '/' and the '='
    int separatorCount = count(buffer.begin(), buffer.end(), '/');
    int pos1 = 0;
    int pos2 = 0;
    string acronym; //holds the acronym for determining what each entry describes
    string temp;
    for (int i = 0; i < separatorCount; i++) {
        pos1 = buffer.find('/', pos2);
        pos2 = buffer.find('=', pos1);
        acronym = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));

        pos1 = pos2;
        pos2 = buffer.find(' ', pos1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        //essentially a switch statement on the acronym parsed out
        if (acronym == "No") {
            fileIndex = atoi(temp.c_str());
        } else if (acronym == "LP") {
            startPosition = atoi(temp.c_str());
        } else if (acronym == "US") {
            upstreamLeg = removeParens(temp);
        } else if (acronym == "DS") {
            downstreamLeg = removeParens(temp);
        } else if (acronym == "B") {
            bulb = temp;
        } else if (acronym == "T") {
            tail = temp;
        } else if (acronym == "USL") {
            upstreamLegLen = atoi(temp.c_str());
        } else if (acronym == "DSL") {
            downstreamLegLen = atoi(temp.c_str());
        } else if (acronym == "SL") {
            stemLen = atoi(temp.c_str());
        } else if (acronym == "BL") {
            bulbLen = atoi(temp.c_str());
        } else if (acronym == "Mm") {
            numMismatch = atoi(temp.c_str());
        } else if (acronym == "Gp") {
            numGap = atoi(temp.c_str());
        } else if (acronym == "DG") {
            deltaG = atof(temp.c_str());
        } else if (acronym == "G") {
            geneName = temp;
        } else if (acronym == "G>") {
            geneStartPos = atoi(temp.c_str());
        } else if (acronym == "G<") {
            geneStopPos = atoi(temp.c_str());
        } else if (acronym == "DS>") {
            stopToStartDist = atoi(temp.c_str());
        } else if (acronym == "DS<") {
            stopToStopDist = atoi(temp.c_str());
        } else if (acronym == "DM") {
            stopToMidDist = atoi(temp.c_str());
        } else if (acronym == "Seq>") {
            seqStartPos = atoi(temp.c_str());
        } else if (acronym == "Seq") {
            sequence = temp;
        } else if (acronym == "D") {
            pos2 = buffer.length();
            description = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        } else {
        }
    }
    return true;
}

void Terminators::searchForAbortiveOverlaps()
{
    int terminatorLength = (downstreamLegLen + upstreamLegLen + bulbLen);
    //for all abortive matches found previously
    for (int i = 0; i < numAbortiveMatches; i++) {
        //check if the abortive match starts anywhere between
        //terminatorStartSite - abortiveLength and terminatorEndSite
        //NOTE: FORWARD AND REVERSE MATCHING IS DIFFERENT; start position is 5'
        //on both strands, not just based on the forward strand
        if (matchList[i].strand == "forward"
            && matchList[i].startPosition >= startPosition - abortiveLength + 1
            && matchList[i].startPosition <= startPosition + terminatorLength - 1
            && matchList[i].strand == strand) {
            //print information about the sequence match
            outputFile << "Match:\t\t"
                       << matchList[i].startPosition << "\t"
                       << matchList[i].startPosition + abortiveLength - 1 << "\t"
                       << matchList[i].sDev << "\t"
                       << min((matchList[i].startPosition + abortiveLength - 1), (startPosition + terminatorLength - 1)) - max(matchList[i].startPosition, startPosition) + 1 << "\t"
                       << matchList[i].bindSite << "\t" << matchList[i].abortive << endl
                       << "Transcript:\t"
                       << matchList[i].operon << "\t"
                       << matchList[i].transcriptName << "\t"
                       << matchList[i].leftBound << "\t"
                       << matchList[i].rightBound << "\t"
                       << matchList[i].strand << endl;
            outputFile << "--------------------------------------------------" << endl;
            //print information about the terminator
            outputFile << "Terminator:\t"
                       << geneName << "\t"
                       << strand << "\t"
                       << startPosition << "\t"
                       << startPosition + upstreamLegLen + downstreamLegLen + bulbLen - 1 << "\t"
                       << upstreamLeg << bulb << downstreamLeg << "\t"
                       << deltaG << endl << endl;
        }else if (matchList[i].strand == "reverse"
            && matchList[i].startPosition >= startPosition - terminatorLength - abortiveLength + 2
            && matchList[i].startPosition <= startPosition
            && matchList[i].strand == strand) {
            //print information about the sequence match
            outputFile << "Match:\t\t"
                       << matchList[i].startPosition << "\t"
                       << matchList[i].startPosition + abortiveLength - 1 << "\t"
                       << matchList[i].sDev << "\t"
                       << min((matchList[i].startPosition + abortiveLength - 1), (startPosition)) - max(matchList[i].startPosition, (startPosition - terminatorLength + 1)) + 1 << "\t"
                       << matchList[i].bindSite << "\t" << matchList[i].abortive << endl
                       << "Transcript:\t"
                       << matchList[i].operon << "\t"
                       << matchList[i].transcriptName << "\t"
                       << matchList[i].leftBound << "\t"
                       << matchList[i].rightBound << "\t"
                       << matchList[i].strand << endl;
            outputFile << "--------------------------------------------------" << endl;
            //print information about the terminator
            outputFile << "Terminator:\t"
                       << geneName << "\t"
                       << strand << "\t"
                       << startPosition - terminatorLength + 1 << "\t"
                       << startPosition << "\t"
                       << upstreamLeg << bulb << downstreamLeg << "\t"
                       << deltaG << endl << endl;
        }
    }
}
