"""
python script.py *.seq > output
"""

from sys import argv, stdout

d = []
fns = []
for fn in argv[1:]:
    with open(fn) as f:
        line = f.readline().rstrip()
    fns.append(fn[:-4])
    d.append(list(line))


def util(v, j, out=stdout):
    a = ""
    for i in v:
        if i == "N":
            continue
        if a == "":
            a = i
            continue
        if a != i:
            out.write(str(j) + "\t" + "\t".join(v) + "\n")
            return


stdout.write("\t".join(fns) + "\n")

for i in range(len(d[0])):
    util([j[i] for j in d], i+1)
