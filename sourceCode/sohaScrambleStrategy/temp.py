import csv

for i in range(4, 16):
    with open('wob{i}.dg'.format(i=i), 'r') as sigCSV:
        sigMatches = []
        reader = csv.reader(sigCSV, delimiter=' ')
        for row in reader:
            if float(row[-1]) < -7:
                sigMatches.append(row[0] + row[1] + row[2] + row[3] + row[4])
        # print(sigMatches)
    with open('wobbleProbabilities{i}'.format(i=i), 'r') as probCSV:
        expected = 0
        measured = 0
        sigexpected = 0
        sigmeasured = 0
        reader = csv.reader(probCSV, delimiter=' ')
        for row in reader:
            # if row[5] == 'forward':
            #     tss = row[3]
            # elif row[5] == 'reverse':
            #     tss = row[4]
            # else:
            #     exit(1)
            matchID = row[1] + row[2] + row[3] + row[4] + row[5]
            if matchID in sigMatches:
                # print(matchID, row[0], row[-1])
                sigexpected += float(row[0])
                sigmeasured += float(row[-1])
            expected += float(row[0])
            measured += float(row[-1])
        print("length", i)
        print("sig:", sigexpected, sigmeasured)
        print("all:", expected, measured)
