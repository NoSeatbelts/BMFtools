from MawCluster.BCFastq import FastFisherFlattening, MakeTagComment
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, permuteNucleotides

import numpy as np
import multiprocessing as mp

import argparse
import pysam

def shadeRead(read, BC, head=0):
    if("N" in BC):
        read.comment += "|FP=0|BS=%s" % (read.sequence[1:head + 1] + BC)
    else:
        read.comment += "|FP=1|BS=%s" % (read.sequence[1:head + 1] + BC)
    return read

def worker(Fastq, IndexFq, bcLen, prefix, head=0):
    """
    class based implementation above probably won't work.  switching
    to a function-based one, using a worker function.  I should probably
    have one job that queues workers then a worker function that actually
    executes, going to need to compare pool vs queue.
    """
    fastq = pFastqFile(Fastq)
    indexFq = pFastqFile(IndexFq)
    bcHash = {}
    ifn = indexFq.next
    lenPrefix = len(prefix)
    print("now starting prefix %s" % (prefix))
    for read in fastq:
        BC = ifn().sequence
        if prefix == BC[:lenPrefix]:
            try:
                bcHash[BC].append(shadeRead(read, BC))
            except KeyError:
                bcHash[BC]=[shadeRead(read, BC)]
    fqname = Fastq.split('.')[0]
    output = open(fqname+"."+prefix+".fastq",'w')
    for barcode in bcHash.keys():
        output.write(FastFisherFlattening(bcHash[barcode], barcode))
    output.close()
    if bcHash:
        numBC = len(bcHash.keys())
    else:
        numBC = 0
    del bcHash
    print("prefix %s complete" % (prefix))
    return numBC

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("Fastq1", help="Fastq1", type=str)
    parser.add_argument("Fastq2", help="Fastq2", type=str)
    parser.add_argument("IndexFq", help="Barcode Fastq", type=str)
    parser.add_argument("ncpus", help="number of cpus", type=int)
    parser.add_argument("lenPrefix", help="hacky, temporary length of prefix"
                        " by which we will split the jobs",
                        type=int)
    return parser.parse_args()

if __name__ == '__main__':
    args = getArgs()
    pool = mp.Pool(processes=args.ncpus)
    bcLen = len(pFastqFile(args.IndexFq).next().sequence)
    allPrefixes = permuteNucleotides(4**args.lenPrefix,kmerLen=args.lenPrefix)
    results = [pool.apply_async(worker, args=(args.Fastq1, args.IndexFq, bcLen,
                p)) for p in allPrefixes]
    things = [p.get() for p in results]
    print things
