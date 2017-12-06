import csv
import regex
from random import shuffle
from itertools import permutations, product
import statistics
from tqdm import tqdm
import gc

sequenceMatch = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

wobbleMatch = {
    'A': 'T',
    'T': '[AG]',
    'G': '[CT]',
    'C': 'G'
}


def bruteforce(charset, length):
    return (''.join(candidate) for candidate in product(charset, repeat=length))


class SequenceUnit(object):
    def __init__(self, operon, transunit, codingstrand, abortseq325, codingseq523, abortivelength):
        self.operon = operon
        self.transUnit = transunit
        self.codingStrand = codingstrand
        self.abortSeq325 = abortseq325[len(abortseq325) - abortivelength:]
        self.abortiveLength = abortivelength
        self.codingSeq523 = codingseq523
        self.baseDistribution = {base: self.codingSeq523.count(base) for base in ('A', 'T', 'G', 'C')}
        self.trueMatches = []
        self.scrambleMatches = []
        self.numScrambles = 0

    def findTrueMatches(self):
        self.trueMatches = self._findMatches(self.abortSeq325)

    def scrambleMatch(self):
        maxIterations = 1000
        permutationsList = list(permutations(self.abortSeq325))
        shuffle(permutationsList)
        permutationsList = permutationsList[:maxIterations]
        self.numScrambles = maxIterations
        self.scrambleMatches.clear()
        for i in range(maxIterations):
            try:
                scrambledAbortive = ''.join(permutationsList[i])
            except IndexError:
                self.numScrambles = i
                break
            else:
                scrambleMatches = self._findMatches(scrambledAbortive)
                self.scrambleMatches.append(len(scrambleMatches))

    def _findMatches(self, abortseq325):
        bindSiteSeq523 = ''.join([wobbleMatch[base] for base in abortseq325])
        return regex.findall(bindSiteSeq523, self.codingSeq523, pos=len(abortseq325), overlapped=True)

    def __str__(self):
        return " ".join([self.operon, self.transUnit, self.codingStrand, str(self.abortiveLength), self.abortSeq325,
                         str(self.numScrambles), str(statistics.mean(self.scrambleMatches)),
                         str(statistics.median(self.scrambleMatches)), str(statistics.stdev(self.scrambleMatches)),
                         str(len(self.trueMatches))])


abortiveLength = 4

with open('allSequences15', 'r') as seqCSV:
    seqReader = csv.reader(seqCSV, delimiter=' ')
    seqUnits = [SequenceUnit(*row, abortiveLength) for row in seqReader]
    for seqUnit in tqdm(seqUnits):
        seqUnit.findTrueMatches()
        seqUnit.scrambleMatch()
        gc.collect()

with open('scrambleResultsWobble{al}'.format(al=abortiveLength), 'w') as outFile:
    outFile.write('operon transUnit codingStrand abortiveLength abortSeq325 numScrambles scrambleMean scrambleMedian scrambleStdev numTrueMatches\n')
    for seqUnit in seqUnits:
        outFile.write(str(seqUnit))
        outFile.write('\n')
