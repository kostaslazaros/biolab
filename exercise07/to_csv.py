def parse_real(fname):
    """Parse file with uniprot/validated real data"""
    lns = []
    with open(fname) as fil:
        for lin in fil.readlines():
            lns.append(lin.strip())
    votes = ""
    for i in range(len(lns[0])):
        v0 = 0
        for j in range(2, 7):
            v0 += 1 if lns[j][i] == "M" else 0
        if v0 >= 3:
            votes += '1'
        else:
            votes += '0'
    lns.append(votes)
    for lin in lns:
        print(','.join(lin))


if __name__ == "__main__":
    parse_real("transpred_res.txt")
