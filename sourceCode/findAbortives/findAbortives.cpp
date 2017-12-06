#include "findAbortives.h"

///////////////////////////////
// General Purpose Functions //
///////////////////////////////

//compare two TranscriptUnits based on their standard deviations (for sorting)
bool devCompare(const TranscriptUnit &lhs, const TranscriptUnit &rhs) {
    return (lhs.sDev > rhs.sDev);
}

//reverse complement a given string
string revComp(string s)
{
    string rVal = s;
    //reverse
    for (unsigned i = 0; i < s.length(); i++) {
        rVal[i] = s[s.length() - i - 1];
    }
    //complement
    for (unsigned i = 0; i < rVal.length(); i++) {
        if (rVal[i] == 'a') {
            rVal[i] = 't';
        } else if (rVal[i] == 't') {
            rVal[i] = 'a';
        } else if (rVal[i] == 'g') {
            rVal[i] = 'c';
        } else if (rVal[i] == 'c') {
            rVal[i] = 'g';
        }
    }
    return rVal;
}

//reverse a given string
string rev(string s)
{
    string rVal = s;
    //reverse
    for (unsigned i = 0; i < s.length(); i++) {
        rVal[i] = s[s.length() - i - 1];
    }
    return rVal;
}

//////////////////////////////
// GenomeSequence Functions //
//////////////////////////////

void GenomeSequence::populate(char *fileName)
{
    ifstream geneFile(fileName);
    string buffer;
    int i;
    wholeSequence = "";
    //for progress output
    int progCount = 1;
    cerr << "Loading genome sequence:" << endl;
    while (geneFile >> buffer) {
        for (i = 0; i < 6; i++) {
            geneFile >> buffer;
            wholeSequence = wholeSequence + buffer;
        }
        if (progCount % 500 == 0) {
            cerr << " . " << flush;
        }
        if (progCount % 4000 == 0) {
            cerr << "]\n\33[1A\33[2K" << flush;
        }
        progCount += 1;
    }
    cerr << "]\n\33[1A\33[2K" << flush;
    cerr << "Finished loading." << endl;

}

string GenomeSequence::getTranscript(int begin, int end)
{
    string transcript = wholeSequence.substr(begin, (end - begin) + 1);
    return transcript;
}

//////////////////////////////
// TranscriptList Functions //
//////////////////////////////

TranscriptList::TranscriptList(char *gFile, char *tFile, int nPerm)
{
    abortiveLength = 0;
    numPermutations = nPerm;
    
    //initialize counters to 0
    totalMatchCount = 0;
    transcriptsWithMatches = 0;
    numberOfTranscripts = 0;
    
    //take in the gene sequence
    geneSeq.populate(gFile);

    //take in all of the transcripts
    transcriptFile.open(tFile);
    while (getNextTranscript()) {
    }
}

void TranscriptList::setOutputFile(char *fName)
{
    outputFile.close();
    outputFile.open(fName);
}

bool TranscriptList::getNextTranscript()
{
    string buffer;
    string temp;
    int pos1 = 0;
    int pos2 = 0;

    //file is comma separated (FOR NOW; fix this to be universal?)
    //therefore, find locations of commas and read between them

    //only use a transcript if tss != "" (use a bool for easy checking)
    bool validTSS = false;

    //make a new TranscriptUnit
    TranscriptUnit newTU(numPermutations);
    while (!validTSS) {
        if (!getline(transcriptFile, buffer)) {
            return false;
        }
        
        //get operon name
        pos1 = buffer.find(',');
        pos2 = buffer.find(',', pos1 + 1);
        newTU.operon = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << operon << " ";

        //get transcript unit name
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.tu = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << tu << " ";

        //get promoter name
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.promoter = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << promoter << " ";

        //get transcription start site
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        if (temp != "\"\"") { //if the entry for the start site is empty, we can't use this transcript
            validTSS = true;
        }
        newTU.tss = atoi(temp.c_str());
        // cerr << tss << " ";
    
        //get strand
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.strand = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << strand << " ";

        //get name of first gene
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.fg = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << fg << " ";

        //get left coord of first gene
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.fglCoord = atoi(temp.c_str());
        // cerr << fglCoord << " ";

        //get right coord of first gene
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.fgrCoord = atoi(temp.c_str());
        // cerr << fgrCoord << " ";

        //get name of last gene
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.lg = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << lg << " ";

        //get left coord of last gene
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.lglCoord = atoi(temp.c_str());
        // cerr << lglCoord << " ";

        //get right coord of last gene
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.lgrCoord = atoi(temp.c_str());
        // cerr << lgrCoord << " ";

        //get terminator type
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.trm = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << trm << " ";

        //get left coord of terminator
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.trmlCoord = atoi(temp.c_str());
        // cerr << trmlCoord << " ";

        //get right coord of terminator
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.trmrCoord = atoi(temp.c_str());
        // cerr << trmrCoord << " ";

        //get left coord of untranslated region (transcript)
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.utrlCoord = atoi(temp.c_str());
        // cerr << utrlCoord << " ";

        //get right coord of untranslated region (transcript)
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.utrrCoord = atoi(temp.c_str());
        // cerr << utrrCoord << " ";

        //get left coord of 5' untranslated region
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.flCoord = atoi(temp.c_str());
        // cerr << flCoord << " ";

        //get right coord of 5' untranslated region
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.frCoord = atoi(temp.c_str());
        // cerr << frCoord << " ";

        //get 5' utr sequence
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.fSeq = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << fSeq << " ";

        //get left coord of 3' untranslated region
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.tlCoord = atoi(temp.c_str());
        // cerr << tlCoord << " ";

        //get right coord of 3' untranslated region
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        temp = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        newTU.trCoord = atoi(temp.c_str());
        // cerr << trCoord << " ";

        //get 3' utr sequence
        pos1 = pos2;
        pos2 = buffer.find(',', pos1 + 1);
        newTU.tSeq = buffer.substr((pos1 + 1), (pos2 - pos1 - 1));
        // cerr << tSeq << " ";

        // cerr << endl;
    }
    numberOfTranscripts += 1;
    tUnits.push_back(newTU);
    return true;
}

void TranscriptList::setTranscriptSequences()
{
    for (int i = 0; i < numberOfTranscripts; i++) {
        //get forward sequence directly from genome file sequence
        tUnits[i].forwardSequence = geneSeq.getTranscript(tUnits[i].utrlCoord - 1, tUnits[i].utrrCoord - 1);
        //store reverse sequence by reverse complementing the forward sequence
        tUnits[i].reverseSequence = revComp(tUnits[i].forwardSequence);
    }
}

void TranscriptList::setAbortiveLengths(int len)
{
    for (int i = 0; i < numberOfTranscripts; i++) {
        tUnits[i].setAbortiveLength(len);
    }

    abortiveLength = len;
}

void TranscriptList::searchForSequenceMatches(int allowedMismatches)
{
    for (int i = 0; i < numberOfTranscripts; i++) {
        cerr << "]\n\33[1A\33[2K" << flush;
        cerr << "Performing search " << i + 1;
        tUnits[i].sequenceMatchSearch(allowedMismatches);
        totalMatchCount += tUnits[i].matchCount;
        if (tUnits[i].matchCount > 0) {
            transcriptsWithMatches += 1;

            //perform permuted search if matches were found
            tUnits[i].permutedSequenceSearch(allowedMismatches, numPermutations);
            tUnits[i].calculateSignificance(numPermutations);
        }
    }
    cerr << "]\n\33[1A\33[2K" << flush;
    cerr << "Done" << endl;
}

void TranscriptList::searchForWobbleMatches(int allowedMismatches)
{
    for (int i = 0; i < numberOfTranscripts; i++) {
        cerr << "]\n\33[1A\33[2K" << flush;
        cerr << "Performing search " << i + 1;
        tUnits[i].wobbleMatchSearch(allowedMismatches);
        totalMatchCount += tUnits[i].matchCount;
        if (tUnits[i].matchCount > 0) {
            transcriptsWithMatches += 1;

            //perform permuted search if matches were found
            tUnits[i].permutedWobbleSearch(allowedMismatches, numPermutations);
            tUnits[i].calculateSignificance(numPermutations);
        }
    }
    cerr << "]\n\33[1A\33[2K" << flush;
    cerr << "Done" << endl;
}

void TranscriptList::printMatchStatistics()
{
    outputFile << abortiveLength << " BASE ABORTIVES" << endl
               << "Total Number of Transcripts: " << numberOfTranscripts << endl
               << "Transcripts With Matches:    " << transcriptsWithMatches << endl
               << "Total Match Count:           " << totalMatchCount << endl << endl;
}

void TranscriptList::printAllMatches()
{
    // for (int i = 0; i < numberOfTranscripts; i++) {
    //     if (tUnits[i].matchCount > 0) {
    //         tUnits[i].matchCount = 1;
    //     }
    // }
    
    for (int i = 0; i < numberOfTranscripts; i++) {
        for (int j = 0; j < tUnits[i].matchCount; j++) {
            outputFile << tUnits[i].operon
                       << " " << tUnits[i].tu
                       << " " << tUnits[i].utrlCoord
                       << " " << tUnits[i].utrrCoord
                       << " " << tUnits[i].strand
                       << " " << tUnits[i].matchStartSites[j]
                       << " " << tUnits[i].matchSequences[j]
                       << " " << tUnits[i].abortiveSequence
                       // << " " << tUnits[i].sDev
                       << " " << 0 // to prevent 'inf' from messing up further scripting... REMOVE LATER
                       << endl;
        }
    }
}

void TranscriptList::reset()
{
    //reset list counters
    totalMatchCount = 0;
    transcriptsWithMatches = 0;
    abortiveLength = 0;

    //reset counters for each TU
    for (int i = 0; i < numberOfTranscripts; i++) {
        tUnits[i].matchCount = 0;
        tUnits[i].matchStartSites.clear();
        tUnits[i].matchSequences.clear();
        for (int j = 0; j < numPermutations; j++) {
            tUnits[i].permutedMatchCounts[j] = 0;
        }
        tUnits[i].sDev = 0;
    }
}

void TranscriptList::sortBySDev()
{
    sort(tUnits.begin(), tUnits.end(), devCompare);
}


//////////////////////////////
// TranscriptUnit Functions //
//////////////////////////////

TranscriptUnit::TranscriptUnit(int numPermutations)
{
    //initialize all match counters to 0
    matchCount = 0;
    permutedMatchCounts = new int[numPermutations];
    for (int i = 0; i < numPermutations; i++) {
        permutedMatchCounts[i] = 0;
    }
}

void TranscriptUnit::setAbortiveLength(int len)
{
    abortiveLength = len;

    //store the abortive sequence of the new length
    if (strand == "forward") {
        abortiveSequence = forwardSequence.substr(0, abortiveLength);
        abortiveSequence = rev(abortiveSequence);
    } else if (strand == "reverse") {
        abortiveSequence = reverseSequence.substr(0, abortiveLength);
        abortiveSequence = rev(abortiveSequence);
    } else {
        cerr << "ERROR: unknown strand identifier" << endl;
        exit(1);
    }
}

void TranscriptUnit::sequenceMatchSearch(int allowedMismatches)
{
    //since both strands are stored 5' to 3', the logic for the searches is the
    //same; however, the strand to search through is different
    if (strand == "forward") {
        //self-complementary matches don't matter, as RNAP prevents binding
        //for lack of a better starting place, begin the search immediately
        //following the abortive sequence
        for (unsigned i = abortiveLength; i <= forwardSequence.length() - abortiveLength; i++) {
            if (validSequenceMatch(abortiveSequence,
                                   forwardSequence.substr(i, abortiveLength),
                                   allowedMismatches)) {
                matchStartSites.push_back(i + utrlCoord);
                matchCount += 1;
            }
        }
    } else if (strand == "reverse") {
        for (unsigned i = abortiveLength; i <= reverseSequence.length() - abortiveLength; i++) {
            if (validSequenceMatch(abortiveSequence,
                                 reverseSequence.substr(i, abortiveLength),
                                 allowedMismatches)) {
                matchStartSites.push_back(
                    (utrlCoord + reverseSequence.length() - 1) -
                    (i + abortiveLength - 1));
                matchCount += 1;
            }
        }
    }
}

bool TranscriptUnit::validSequenceMatch(string abortive, string target, int allowedMismatches)
{
    int mismatchCount = 0;
    
    for (unsigned i = 0; i < target.length(); i++) {
        switch(target[i]) {
        case 'a':
            if (abortive[i] != 't') {
                mismatchCount += 1;
            }
            break;
        case 't':
            if (abortive[i] != 'a') {
                mismatchCount += 1;
            }
            break;
        case 'c':
            if (abortive[i] != 'g') {
                mismatchCount += 1;
            }
            break;
        case 'g':
            if (abortive[i] != 'c') {
                mismatchCount += 1;
            }
            break;
        default:
            cerr << endl << "ERROR: Unknown base in sequence" << endl;
            exit(1);
            break;
        }
    }

    if (mismatchCount <= allowedMismatches) {
        matchSequences.push_back(target);
        return true;
    } else {
        return false;
    }
}

void TranscriptUnit::wobbleMatchSearch(int allowedMismatches)
{
    //since both strands are stored 5' to 3', the logic for the searches is the
    //same; however, the strand to search through is different
    if (strand == "forward") {
        //self-complementary matches don't matter, as RNAP prevents binding
        //for lack of a better starting place, begin the search immediately
        //following the abortive sequence
        for (unsigned i = abortiveLength; i <= forwardSequence.length() - abortiveLength; i++) {
            if (validWobbleMatch(abortiveSequence,
                                   forwardSequence.substr(i, abortiveLength),
                                   allowedMismatches)) {
                matchStartSites.push_back(i + utrlCoord);
                matchCount += 1;
            }
        }
    } else if (strand == "reverse") {
        for (unsigned i = abortiveLength; i <= reverseSequence.length() - abortiveLength; i++) {
            if (validWobbleMatch(abortiveSequence,
                                 reverseSequence.substr(i, abortiveLength),
                                 allowedMismatches)) {
                matchStartSites.push_back(
                    (utrlCoord + reverseSequence.length() - 1) -
                    (i + abortiveLength - 1));
                matchCount += 1;
            }
        }
    }
}

bool TranscriptUnit::validWobbleMatch(string abortive, string target, int allowedMismatches)
{
    int mismatchCount = 0;
    
    for (unsigned i = 0; i < target.length(); i++) {
        switch(target[i]) {
        case 'a':
            if (abortive[i] != 't') {
                mismatchCount += 1;
            }
            break;
        case 't':
            if (abortive[i] != 'a' && abortive[i] != 'g') {
                mismatchCount += 1;
            }
            break;
        case 'c':
            if (abortive[i] != 'g') {
                mismatchCount += 1;
            }
            break;
        case 'g':
            if (abortive[i] != 'c' && abortive[i] != 't') {
                mismatchCount += 1;
            }
            break;
        default:
            cerr << endl << "ERROR: Unknown base in sequence" << endl;
            exit(1);
            break;
        }
    }

    if (mismatchCount <= allowedMismatches) {
        matchSequences.push_back(target);
        return true;
    } else {
        return false;
    }
}

void TranscriptUnit::permutedSequenceSearch(int allowedMismatches, int numPermutations)
{
    if (strand == "forward") {
        permutedForward = forwardSequence;
        for (int k = 0; k < numPermutations; k++) {
            permuteTranscriptSequence();
            for (unsigned i = abortiveLength; i < permutedForward.length() - abortiveLength; i++) {
                if (validSequenceMatch(abortiveSequence,
                                       permutedForward.substr(i, abortiveLength),
                                       allowedMismatches)) {
                    permutedMatchCounts[k] += 1;
                }
            }
        }
    } else if (strand == "reverse") {
        permutedReverse = reverseSequence;
        for (int k = 0; k < numPermutations; k++) {
            permuteTranscriptSequence();
            for (unsigned i = abortiveLength; i < permutedReverse.length() - abortiveLength; i++) {
                if (validSequenceMatch(abortiveSequence,
                                       permutedReverse.substr(i, abortiveLength),
                                       allowedMismatches)) {
                    permutedMatchCounts[k] += 1;
                }
            }
        }
    }
}

void TranscriptUnit::permutedWobbleSearch(int allowedMismatches, int numPermutations)
{
    if (strand == "forward") {
        permutedForward = forwardSequence;
        for (int k = 0; k < numPermutations; k++) {
            permuteTranscriptSequence();
            for (unsigned i = abortiveLength; i < permutedForward.length() - abortiveLength; i++) {
                if (validWobbleMatch(abortiveSequence,
                                       permutedForward.substr(i, abortiveLength),
                                       allowedMismatches)) {
                    permutedMatchCounts[k] += 1;
                }
            }
        }
    } else if (strand == "reverse") {
        permutedReverse = reverseSequence;
        for (int k = 0; k < numPermutations; k++) {
            permuteTranscriptSequence();
            for (unsigned i = abortiveLength; i < permutedReverse.length() - abortiveLength; i++) {
                if (validWobbleMatch(abortiveSequence,
                                       permutedReverse.substr(i, abortiveLength),
                                       allowedMismatches)) {
                    permutedMatchCounts[k] += 1;
                }
            }
        }
    }
}

void TranscriptUnit::permuteTranscriptSequence()
{
    int j, tLen;
    srand(time(NULL));
    //Durnstenfeld's implementation of the Fisher-Yates shuffle algorithm
    //for all of the bases except for the abortive, perform swaps
    if (strand == "forward") {
        tLen = permutedForward.length();
        for (int i = tLen - 1; i > abortiveLength; i--) {
            j = rand() % (i - abortiveLength);
            swap(permutedForward[i], permutedForward[j + abortiveLength]);
        }
        //bug catching; ensure that the abortive is unaffected
        if (rev(abortiveSequence) != permutedForward.substr(0, abortiveLength)) {
            cerr << "Crashing in permuteTranscriptSequence: " << abortiveSequence << " " << permutedForward.substr(0, abortiveLength) << endl;
            exit(1);
        }
    } else if (strand == "reverse") {
        tLen = permutedReverse.length();
        for (int i = tLen - 1; i > abortiveLength; i--) {
            j = rand() % (i - abortiveLength);
            swap(permutedReverse[i], permutedReverse[j + abortiveLength]);
        }
        //bug catching; ensure that the abortive is unaffected
        if (rev(abortiveSequence) != permutedReverse.substr(0, abortiveLength)) {
            cerr << "Crashing in permuteTranscriptSequence: " << abortiveSequence << " " << permutedReverse.substr(0, abortiveLength) << endl;
            exit(1);
        }
    }
}

void TranscriptUnit::calculateSignificance(int numPermutations)
{
    sDev = jSTT(permutedMatchCounts, numPermutations, matchCount);
}
