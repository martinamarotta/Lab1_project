import sys

def get_ids(idlist):
    f = open(idlist)
    return f.read().strip().split("\n")

def get_seq(pidlist,seqfile):
    f = open(seqfile)
    s = 0
    for line in f:
        if line.startswith(">"):
            pid = line.split("|")[1]
            s = 0
            if pid in pidlist : s = 1
        if s == 1:  print(line.strip())


if __name__ == "__main__":
    idlist = sys.argv[1]
    seqfile = sys.argv[2]
    pidlist = get_ids(idlist)
    get_seq(pidlist, seqfile)