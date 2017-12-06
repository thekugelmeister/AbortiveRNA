from fishers_exact_test import FishersExactTest
import csv
from collections import Counter


def is_int(c):
    try:
        int(c)
        return True
    except ValueError:
        return False

# First read in terminator files:
terminators = []
# skipped_terminators = []
with open('E_coli_MG1655_term_RawData\89318805\greatestdGreg', 'r') as fTermFile:
    fTermFile.__next__()
    fTermReader = csv.reader(fTermFile, delimiter='/')
    # Multiple genes called insL; label them for ease of comparison between files
    insLCount = 1
    for row in fTermReader:
        termData = {}
        for field in row[1:]:
            kv = field.split(sep="=", maxsplit=1)
            try:
                termData[kv[0]] = kv[1].rstrip()
            except IndexError:
                pass
        if termData["G"] == "insL":
            termData["G"] = "insL-{i}".format(i=insLCount)
            insLCount += 1
        termData["strand"] = "forward"
        termData["length"] = int(termData["USL"]) + int(termData["DSL"]) + int(termData["BL"])
        termData["lBound"] = int(termData["LP"])
        termData["rBound"] = int(termData["LP"]) + termData["length"] - 1
        # if any([is_int(c) for c in termData["G"]]):
        #     skipped_terminators.append(termData)
        # else:
        terminators.append(termData)

with open('E_coli_MG1655_term_RawData\89318805\greatestdGcmp', 'r') as rTermFile:
    rTermFile.__next__()
    rTermReader = csv.reader(rTermFile, delimiter='/')
    for row in rTermReader:
        termData = {}
        for field in row[1:]:
            kv = field.split(sep="=", maxsplit=1)
            try:
                termData[kv[0]] = kv[1]
            except IndexError:
                pass
        termData["strand"] = "reverse"
        termData["length"] = int(termData["USL"]) + int(termData["DSL"]) + int(termData["BL"])
        termData["rBound"] = int(termData["LP"])
        termData["lBound"] = int(termData["LP"]) - termData["length"] + 1
        # if any([is_int(c) for c in termData["G"]]):
        #     skipped_terminators.append(termData)
        # else:
        terminators.append(termData)

print("Num Terminators:", len(terminators))
# print("Num Skipped Terms:", len(skipped_terminators))

with open('RBSSet.txt', 'r') as RBSFile:
    fieldnames = ['id', 'gene', 'leftBound', 'rightBound', 'strand', 'centerPos', 'sequence', 'evidence']
    RBSReader = csv.DictReader(RBSFile, fieldnames=fieldnames, delimiter='\t')
    RBSs = [row for row in RBSReader]
for rbs in RBSs:
    rbs["leftBound"] = int(rbs["leftBound"])
    rbs["rightBound"] = int(rbs["rightBound"])

# 1: get len(transcript) - len(RBS) - len(terminator) for each transcript
# 2: get number of matches in each transcript
# 3: get number of matches in RBS + number of matches in terminator for each transcript
# 4: get (2) - (3)
# 5: run test
# 6: calculate significance
for aLen in range(4, 5):
    # get terminator matches per transcript; keyed by operon+tUnit+TSS
    with open("wob{i}_trmout.csv".format(i=aLen), "r") as termFile:
        termReader = csv.reader(termFile, delimiter=" ")
        termMatches = Counter([row[0] + row[1] + str(row[-3]) for row in termReader])
    print("Num TUs with Term Matches:", len(termMatches))

    # TODO: Need to figure out how many RBSs are in each transcriptional unit, then multiply by 6
    # get RBS matches per transcript; keyed by operon+tUnit+TSS
    with open("wob{i}_RBSMatches.csv".format(i=aLen), "r") as rbsFile:
        rbsReader = csv.reader(rbsFile, delimiter=" ")
        rbsMatches = Counter([row[0] + row[1] + str(row[-3]) for row in rbsReader])
    print("Num TUs with RBS Matches:", len(rbsMatches))

    TUs = []
    # skipped_TUs = []
    with open('wobbleProbabilities{i}'.format(i=aLen), 'r') as probFile:
        fieldnames = ["prob", "operon", "tUnit", "lBound", "rBound", "codingStrand", "totalMatches"]
        allMatchReader = csv.DictReader(probFile, fieldnames=fieldnames, delimiter=' ')
        for row in allMatchReader:
            # store tu length
            row["tuLength"] = int(row["rBound"]) - int(row["lBound"]) + 1
            row["rBound"] = int(row["rBound"])
            row["lBound"] = int(row["lBound"])
            row["totalMatches"] = int(row["totalMatches"])
            # # skipping any TU with numbers in its name
            # if any([is_int(c) for c in row["tUnit"]]):
            #     skipped_TUs.append(row)
            # else:
            # store signatures of all genes in tu (for finding terminators)
            # genes = []
            # for s in row["tUnit"].split("-"):
            #     if len(s) == 3:
            #         genes.append(s)
            #     else:
            #         tlCode = s[:3]
            #         rest = s[3:]
            #         for c in rest:
            #             genes.append(tlCode + c)
            # row["genes"] = genes
            row["terminators"] = []
            row["RBSs"] = []
            row["termMatches"] = 0
            row["rbsMatches"] = 0
            TUs.append(row)
    print("Num TUs:", len(TUs))
    # print("Num Skipped TUs:", len(skipped_TUs))

    for term in terminators:
        for tu in TUs:
            # if term["G"] in tu["genes"]:
            if (tu["lBound"] < term["lBound"] < tu["rBound"] or tu["lBound"] < term["rBound"] < tu["rBound"]) and tu["codingStrand"] == term["strand"]:
                tu["terminators"].append(term)
    for tu in TUs:
        tu["totalTermLength"] = sum([term["length"] for term in tu["terminators"]])
    print("Num TUs with 1+ Terms:", len([tu for tu in TUs if len(tu["terminators"]) != 0]))

    for rbs in RBSs:
        for tu in TUs:
            if (tu["lBound"] < rbs["leftBound"] < tu["rBound"] or tu["lBound"] < rbs["rightBound"] < tu["rBound"]) and tu["codingStrand"] == rbs["strand"]:
                tu["RBSs"].append(rbs)
    for tu in TUs:
        tu["totalRBSLength"] = len(tu["RBSs"]) * 6
    print("Num TUs with 1+ RBSs:", len([tu for tu in TUs if len(tu["RBSs"]) != 0]))

    i = 0
    temp = []
    for tu in TUs:
        if tu["codingStrand"] == "forward":
            if tu["operon"] + tu["tUnit"] + str(tu["lBound"]) in termMatches:
                i += 1
                temp.append(tu["operon"] + tu["tUnit"] + str(tu["lBound"]))
                tu["termMatches"] = termMatches[tu["operon"] + tu["tUnit"] + str(tu["lBound"])]
        elif tu["codingStrand"] == "reverse":
            if tu["operon"] + tu["tUnit"] + str(tu["rBound"]) in termMatches:
                i += 1
                temp.append(tu["operon"] + tu["tUnit"] + str(tu["rBound"]))
                tu["termMatches"] = termMatches[tu["operon"] + tu["tUnit"] + str(tu["rBound"])]
    print("Num termMatches retrieved:", i)
    # cTemp = Counter(temp)
    # print([(k, v) for k, v in cTemp.items() if v > 1])

    i = 0
    temp = []
    for tu in TUs:
        if tu["codingStrand"] == "forward":
            if tu["operon"] + tu["tUnit"] + str(tu["lBound"]) in rbsMatches:
                i += 1
                temp.append(tu["operon"] + tu["tUnit"] + str(tu["lBound"]))
                tu["rbsMatches"] = rbsMatches[tu["operon"] + tu["tUnit"] + str(tu["lBound"])]
        elif tu["codingStrand"] == "reverse":
            if tu["operon"] + tu["tUnit"] + str(tu["rBound"]) in rbsMatches:
                i += 1
                temp.append(tu["operon"] + tu["tUnit"] + str(tu["rBound"]))
                tu["rbsMatches"] = rbsMatches[tu["operon"] + tu["tUnit"] + str(tu["rBound"])]
    print("Num RBSMatches retrieved:", i)

    # # Fisher Test for individual TUs
    # i = 0
    # for tu in TUs:
    #     table = [[tu["totalMatches"] - tu["termMatches"],
    #               tu["termMatches"]],
    #               [tu["tuLength"] - tu["totalMatches"] - tu["termMatches"],
    #               tu["totalTermLength"] - tu["termMatches"]]]
    #     if table[1][1] < 0:
    #         print(i, "ERROR: no assigned terminator, but matches exist", tu)
    #     else:
    #         print(i, tu["tUnit"], table, tu["termMatches"], end=" ")
    #         fet = FishersExactTest(table)
    #         try:
    #             print(fet.left_tail_p(), fet.right_tail_p(), fet.two_tail_p())
    #         except OverflowError:
    #             print("OFE")
    #         except AssertionError:
    #             print("AE")
    #     i += 1

    # Fisher Test for individual RBSs
    i = 0
    for tu in TUs:
        table = [[tu["totalMatches"] - tu["rbsMatches"],
                  tu["rbsMatches"]],
                  [tu["tuLength"] - tu["totalMatches"] - tu["rbsMatches"],
                  tu["totalRBSLength"] - tu["rbsMatches"]]]
        if table[1][1] < 0:
            print(i, "ERROR: no assigned RBS, but matches exist", tu)
        else:
            print(i, tu["tUnit"], table, tu["rbsMatches"], end=" ")
            fet = FishersExactTest(table)
            try:
                print(fet.left_tail_p(), fet.right_tail_p(), fet.two_tail_p())
            except OverflowError:
                print("OFE")
            except AssertionError:
                print("AE")
        i += 1

    # # Fisher Test for entire genome
    # totalTUMatches = sum([tu["totalMatches"] for tu in TUs])
    # totalTermMatches = sum(tu["termMatches"] for tu in TUs)
    # totalTULength = sum([tu["tuLength"] for tu in TUs])
    # totalTermLength = sum([tu["totalTermLength"] for tu in TUs])
    # table = [[totalTUMatches - totalTermMatches,
    #           totalTermMatches],
    #          [totalTULength - totalTUMatches - totalTermMatches,
    #           totalTermLength - totalTermMatches]]
    # fet = FishersExactTest(table)
    # print(table)
    # try:
    #     print(fet.left_tail_p(), fet.right_tail_p(), fet.two_tail_p())
    # except OverflowError:
    #     print("OFE")
