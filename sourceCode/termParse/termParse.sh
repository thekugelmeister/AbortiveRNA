# script for running the termParse program on both strands for a range
# of abortive lengths

printf "Match File Identifier: "
read mFile
printf "Using greatestDg terminator files\n"
printf "Shortest Abortive Length: "
read sLen
printf "Longest Abortive Length: "
read lLen

for i in $(seq $sLen $lLen)
do
    ./termParse $mFile$i ../terminators/greatestdGreg forward
    ./termParse $mFile$i ../terminators/greatestdGcmp reverse
done

