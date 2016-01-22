#!/usr/bin/env python
# extract data from VCF file
# VCF.py

"""
python  VCF.py  VCF_file_name  >  output
"""



from sys import argv, stdout, stderr


def _util(line, w=stdout):
    items = line.rstrip().split("\t")
    out = [items[1], items[2]]
    for i in items[9:]:
        temp = i.split(":")[0]
        assert len(temp) == 3
        out.append(temp[0])
        out.append(temp[-1])
    w.write("\t".join(out) + "\n")
    w.flush()
    return


def main():
    with open(argv[1]) as f:
        for line in f:
            if not line.startswith("#"):
                try:
                    _util(line)
                except Exception as e:
                    stderr.write(line)
                    exit(1)
                break
        for line in f:
            try:
                _util(line)
            except Exception as e:
                stderr.write(line)
                exit(1)

if __name__ == '__main__':
    main()

