import sys
import subprocess
import shlex

try:
    import pysam
except ImportError:
    sys.stderr.write("Could not import pysam. Not running tests.\n")
    sys.exit(1)
import numpy as np

def get_tags(read):
    ret = {}
    for tag in read.comment.split("\t")[1:]:
        if tag[3] == 'i':
            ret[tag[:2]] = int(tag.split(":")[2])
        elif tag[3] == 'f':
            ret[tag[:2]] = float(tag.split(":")[2])
        elif tag[3] == 'B':
            ret[tag[:2]] = np.array(tag.split(":")[2].split(",")[1:], dtype=np.uint32)
        else:
            raise NotImplementedError("Haven't finshed this.")
    return ret


def main():
    for ex in ["bmftools_db", "bmftools", "bmftools_p"]:
        cstr = "../../%s hashdmp -o hashdmp_test.out hashdmp_test.fq" % ex
        subprocess.check_call(shlex.split(cstr))
        fqh = pysam.FastqFile("hashdmp_test.out")
        try:
            r1 = fqh.next()
        except:
            r1 = next(fqh)  # Python 3
        tags = get_tags(r1)
        assert tags["FM"] == 7
        try:
            assert round(tags["NF"], 2) == 0.14
        except AssertionError:
            sys.stderr.write("Tag for NF: '%f'\n" % tags["NF"])
            raise
        assert tags["RV"] == 2
        assert tags["DR"]
        assert len(r1.name) == 16
        try:
            r1 = fqh.next()
        except:
            r1 = next(fqh)  # Python 3
        tags = get_tags(r1)
        assert tags["FM"] == 1
        assert tags["FP"] == 0
        assert tags["DR"] == 0

if __name__ == "__main__":
    sys.exit(main())
