cleavelist = set()
with open("netchop_numbered.txt", 'r') as filter_file:
    next(filter_file)
    counter = -1
    for line in filter_file:
        tokens = line.strip().split()
        if int(tokens[0]) == 1:
            counter += 1
        if float(tokens[3]) > 0.1:
            cleavelist.add(int(tokens[5])-1+counter)
with open("protein_no_space.fasta", "r") as f:
        covids = f.read()


#counter = 0
for i in range(len(covids)):
        #counter += 1
        for x in [8,9,10,11,12]:
                #if i>100000:
        #               break
                if (i+1) in cleavelist and (i+1+x) in cleavelist:
                        print(covids[i:i+x])
