##################################################################################################
# dgcalc.py
##################################################################################################
# NOTE: Requires installation of ViennaRNA package (specifically utilizes "RNAcofold")
##################################################################################################
#
# This script takes as an argument a base-filename (e.g. "seq10") and does the following:
# 1. parses "*.trmout.forward", "*.trmout.reverse", "*.RBSMatches", and "*"
# 2. retrieves the names, the 5'-3' binding site sequence, and the 3'-5' abortive sequence
#    of each entry in each file
# 3. reverses the abortive sequences
# 4. for each entry, put the binding site sequence in "bind.seq" and the abortive sequence
#    in "abort.seq"
# 5. run the following command: perl hybrid2.pl --NA=RNA --tmin=0 --tmax=100 --A0=0.00001
#    --B0=0.00001 bind.seq abort.seq
# 6. find the deltaG by subtracting the Tmax and Tmin free energy values found in
#    "bind-abort.ens.dG" and store it in the output file
# 7. remove all the extraneous files and move to the next feature
##################################################################################################
import os, sys
import re

if len(sys.argv) != 2:
    print "ERROR: Expected base-filename as argument"
    sys.exit(1)

# input files
baseName = str(sys.argv[1])
termForward = open(baseName + ".trmout.forward", 'r')
termReverse = open(baseName + ".trmout.reverse", 'r')
rbsMatch = open(baseName + ".RBSMatches", 'r')
baseMatch = open(baseName, 'r')

# output files
# tOFile = open(baseName + ".trmout.dG", 'w')
# rOFile = open(baseName + ".RBSMatches.dG", 'w')
bOFile = open(baseName + ".dG", 'w')

# tOFile.write("Operon TranscriptionUnit Terminator Strand BindingSiteSequence faRNASequence faRNALength DeltaGTerminator DeltaGBinding TSS matchLB matchRB\n")
# # terminator file analysis
# for f in [termForward, termReverse]:
#     lines = f.readlines()
#     numEntries = len(lines) / 5
#     for i in xrange(0,numEntries):
#         # get the required values for the given entry
#         bindSeq = lines[5 * i].split()[5] # 5'-3'
#         abortSeq = lines[5 * i].split()[6] # 3'-5'
#         # reverse abortSeq
#         abortSeq = abortSeq[::-1]
#         # get abortive length
#         abortLen = len(bindSeq)
#         operonName = lines[(5 * i) + 1].split()[1]
#         tranName = lines[(5 * i) + 1].split()[2]
#         strand = lines[(5 * i) + 1].split()[5]
#         termName = lines[(5 * i) + 3].split()[1]
#         termDG = lines[(5 * i) + 3].split()[6]
#         # get TSS based on gene orientation
#         tss = ""
#         if (strand == "forward"):
#             tss = lines[(5 * i) + 1].split()[3]
#         else:
#             tss = lines[(5 * i) + 1].split()[4]
#         matchLB = lines[5 * i].split()[1]
#         matchRB = lines[5 * i].split()[2]
#         # write those vals to the output file
#         tOFile.write(operonName + " " + tranName + " " + termName + " " + strand + " " + bindSeq + " " + abortSeq + " " + str(abortLen) + " " + termDG + " ")

#         with open("test.seq", "w") as seqFile:
#             seqFile.write("{bseq}&{aseq}".format(bseq=bindSeq, aseq=abortSeq))
#         # perform the dG calculation for the given entry
#         os.system("RNAcofold --noPS -T 25 < test.seq > result.txt")
#         # open the dG file and calculate dG as specified above, then write to output file
#         with open("result.txt", "r") as resultFile:
#             dgLine = resultFile.readlines()[1]
#             match = re.match("([\(\.]+)\&([\)\.]+)\s+\(\s*(-?[\d.]+)\)", dgLine)
#             try:
#                 deltaG = match.group(3)
#             except AttributeError:
#                 deltaG = "ERROR"
#                 print "ERROR: " + tranName + " " + termName + " " + strand + " " + bindSeq + " " + abortSeq
#         tOFile.write(deltaG + " " + tss +  " " + matchLB + " " + matchRB + "\n")
#         # clean up directory
#         os.system("rm test.seq result.txt")    


# # RBS file analysis
# lines = rbsMatch.readlines()
# numEntries = (len(lines) - 2) / 6
# rOFile.write("Operon TranscriptionUnit RBS Strand BindingSiteSequence faRNASequence faRNALength DeltaGBinding TSS matchLB matchRB\n")
# for i in xrange(0,numEntries):
#     # get the required values for the given entry
#     bindSeq = lines[6 * i].split()[4] # 5'-3'
#     abortSeq = lines[6 * i].split()[5] # 3'-5'
#     # reverse abortSeq
#     abortSeq = abortSeq[::-1]
#     # get abortive length
#     abortLen = len(bindSeq)
#     operonName = lines[(6 * i) + 1].split()[1]
#     tranName = lines[(6 * i) + 1].split()[2]
#     strand = lines[(6 * i) + 1].split()[5]
#     RBSName = lines[(6 * i) + 3].split()[1]
#     # get TSS based on gene orientation
#     tss = ""
#     if (strand == "forward"):
#         tss = lines[(6 * i) + 1].split()[3]
#     else:
#         tss = lines[(6 * i) + 1].split()[4]
#     matchLB = lines[6 * i].split()[1]
#     matchRB = lines[6 * i].split()[2]
#     rOFile.write(operonName + " " + tranName + " " + RBSName + " " + strand + " " + bindSeq + " " + abortSeq + " " + str(abortLen) + " ")

#     with open("test.seq", "w") as seqFile:
#             seqFile.write("{bseq}&{aseq}".format(bseq=bindSeq, aseq=abortSeq))
#     # perform the dG calculation for the given entry
#     os.system("RNAcofold --noPS -T 25 < test.seq > result.txt")
#     # open the dG file and calculate dG as specified above, then write to output file
#     with open("result.txt", "r") as resultFile:
#         dgLine = resultFile.readlines()[1]
#         match = re.match("([\(\.]+)\&([\)\.]+)\s+\(\s*(-?[\d.]+)\)", dgLine)
#         try:
#             deltaG = match.group(3)
#         except AttributeError:
#             deltaG = "ERROR"
#             print "ERROR: " + tranName + " " + termName + " " + strand + " " + bindSeq + " " + abortSeq
#     rOFile.write(deltaG + " " + tss +  " " + matchLB + " " + matchRB + "\n")
#     # clean up directory
#     os.system("rm test.seq result.txt")

# base file analysis
lines = baseMatch.readlines()
for i in range(0,len(lines)):
    #skip header lines
    if i in range(0,5):
        continue
    # get the required values for the given entry
    bindSeq = lines[i].split()[6] # 5'-3'
    abortSeq = lines[i].split()[7] # 3'-5'
    # reverse abortSeq
    abortSeq = abortSeq[::-1]
    tranName = lines[i].split()[1]
    strand = lines[i].split()[4]
    bOFile.write(tranName + " " + strand + " " + bindSeq + " " + abortSeq + " ")
    with open("test.seq", "w") as seqFile:
        seqFile.write("{bseq}&{aseq}".format(bseq=bindSeq, aseq=abortSeq))
    # perform the dG calculation for the given entry
    os.system("RNAcofold --noPS -T 25 < test.seq > result.txt")
    # open the dG file and calculate dG as specified above, then write to output file
    with open("result.txt", "r") as resultFile:
        dgLine = resultFile.readlines()[1]
        match = re.match("([\(\.]+)\&([\)\.]+)\s+\(\s*(-?[\d.]+)\)", dgLine)
        try:
            deltaG = match.group(3)
        except AttributeError:
            deltaG = "ERROR"
            print "ERROR: " + tranName + " " + termName + " " + strand + " " + bindSeq + " " + abortSeq
    bOFile.write(deltaG + "\n")
    # clean up directory
    os.system("rm test.seq result.txt")    
