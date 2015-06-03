# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True
import logging
from operator import attrgetter as oag
import sys
from subprocess import CalledProcessError, check_call

import pysam
import cython
from cytoolz import map as cmap

from . import BCVCF
from .SNVUtils import GetVCFHeader
from .PileupUtils import (pPileupColumn,
                          GetDiscordantReadPairs, PCInfo)
from utilBMF.HTSUtils import (PysamToChrDict, printlog as pl,
                              ParseBed, PopenDispatcher, PopenCall,
                              parseConfig, TrimExt)
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import ThisIsMadness, FunctionCallException
from MawCluster.SNVUtils cimport VCFPos
from utilBMF.HTSUtils cimport pPileupRead
cimport cython
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport MawCluster.PileupUtils as PileupUtils
ctypedef PileupUtils.PCInfo PCInfo_t
ctypedef PileupUtils.pPileupColumn pPileupColumn_t
ctypedef pPileupRead pPileupRead_t
ctypedef VCFPos VCFPos_t


"""
Programs which write VCFs.
Currently: SNVCrawler.

In development: SV

Settled on two major filters for inclusion in a pileup
1. Min FA (number of family members agreeing) [int=2] (I just can't get
   too big a family...)
2. Min Fraction Agreement Within Family [float=0.6667]

"""


@cython.returns(cython.str)
def SNVCrawler(inBAM,
               cython.str bed="default",
               int minMQ=0,
               int minBQ=0,
               cython.str OutVCF="default",
               cython.float MaxPValue=1e-30,
               cython.bint keepConsensus=False,
               cython.str reference="default",
               cython.bint reference_is_path=False,
               cython.str commandStr="default",
               cython.str fileFormat="default",
               cython.str FILTERTags="default",
               cython.str INFOTags="default",
               cython.str FORMATTags="default",
               cython.bint writeHeader=True,
               cython.float minFracAgreed=0.0,
               int minFA=2,
               cython.str experiment="",
               cython.bint parallel=True,
               sampleName="DefaultSampleName",
               conf="default"):
    cdef int NumDiscordantPairs = 0
    cdef cython.str VCFLineString, VCFString
    cdef pysam.calignmentfile.IteratorColumnRegion ICR
    cdef pysam.calignmentfile.IteratorColumnAllRefs ICAR
#   cdef pysam.calignmentfile.AlignmentFile discPairHandle, inHandle
    cdef pysam.cfaidx.FastaFile refHandle
    cdef list line, discReads, VCFLines, bedlines
    cdef pPileupRead_t i
    cdef pysam.calignmentfile.AlignedSegment read
    if(conf != "default"):
        confDict = parseConfig(conf)
    if(bed != "default"):
        pl("Bed file used: {}".format(bed))
        bedSet = True
        bedlines = ParseBed(bed)
    else:
        if("bed" not in confDict.iterkeys()):
            print("no bed file provided...")
            bedSet = False
        else:
            bedlines = ParseBed(confDict["bed"])
    if(isinstance(bed, list)):
        bedlines = bed
    if("confDict" in locals()):
        try:
            minMQ = int(confDict['minMQ'])
        except KeyError:
            pass
        try:
            minBQ = int(confDict['minBQ'])
        except KeyError:
            pass
        try:
            MaxPValue = float(confDict['MaxPValue'])
        except KeyError:
            pass
        try:
            minFA = int(confDict['minFA'])
        except KeyError:
            pass
        try:
            reference = confDict['ref']
        except KeyError:
            pass
        try:
            minFracAgreed = float(confDict['minFracAgreed'])
        except KeyError:
            pass
    refHandle = pysam.FastaFile(reference)
    if(OutVCF == "default"):
        OutVCF = TrimExt(inBAM) + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    if(OutVCF == "stdout"):
        outHandle = sys.stdout
    else:
        outHandle = open(OutVCF, "w")
    pileupCall = inHandle.pileup
    discPairHandle = pysam.AlignmentFile(
        TrimExt(inBAM) + ".discReadPairs.bam", "wb", template=inHandle)
    ohw = outHandle.write
    dpw = discPairHandle.write
    if(writeHeader):
        try:
            ohw(GetVCFHeader(fileFormat=fileFormat, FILTERTags=FILTERTags,
                             commandStr=commandStr, reference=reference,
                             reference_is_path=False, header=inHandle.header,
                             INFOTags=INFOTags, FORMATTags=FORMATTags))
        except ValueError:
            pl("Looks like the RG header wasn't parseable by pysam - that's u"
               "sually an artefact of the clash between pysam and GATK's ways"
               "of working with RG fields.", level=logging.DEBUG)
            ohw(GetVCFHeader(fileFormat=fileFormat,
                             FILTERTags=FILTERTags,
                             commandStr=commandStr,
                             reference=reference,
                             reference_is_path=False,
                             INFOTags=INFOTags,
                             FORMATTags=FORMATTags))
    if(bedSet):
        pl("Bed file provided - iterating through bed columns")
        for line in bedlines:
            pl("Making pileup call for region %s" % line)
            ICR = pileupCall(line[0], line[1],
                             max_depth=200000,
                             multiple_iterators=False)
            PileupIt = ICR.next
            while True:
                try:
                    pPC = pPileupColumn(PileupIt())
                except StopIteration:
                    pl("Finishing iterations for bed line: %s" % repr(
                        line))
                    break
                except ValueError:
                    pl("Pysam's heinous errors in iteratio.")
                    pl("Region: %s" % repr(line))
                    raise ValueError(repr(line))
                posStr = pPileupColToVCFLines(
                    pPC, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
                    minFracAgreed=minFracAgreed, experiment=experiment,
                    MaxPValue=MaxPValue, refHandle=refHandle,
                    keepConsensus=keepConsensus, reference=reference)
                if(len(posStr) != 0):
                    ohw(posStr + "\n")
                if(pPC.reference_pos > line[2]):
                    pl("Bed region complete.")
                    break
    else:
        ICAR = pileupCall(max_depth=200000, multiple_iterators=False)
        PileupIt = ICAR.next
        while True:
            try:
                VCFString = PileupItToVCFLines(
                    PileupIt(), minMQ=minMQ, minBQ=minBQ, minFA=minFA,
                    minFracAgreed=minFracAgreed, MaxPValue=MaxPValue,
                    experiment=experiment, reference=reference,
                    refHandle=refHandle, keepConsensus=keepConsensus)
            except StopIteration:
                pl("Finished iterations.")
                break
            if(len(VCFString) != 0):
                ohw(VCFString + "\n")
    discPairHandle.close()
    return OutVCF


@cython.returns(tuple)
def PileupItToVCFLines(pysam.calignmentfile.PileupColumn PileupCol,
                       int minMQ=-1, int minBQ=-1,
                       cython.str experiment="",
                       int minFA=-1, cython.float minFracAgreed=-1.,
                       cython.float MaxPValue=-1.,
                       cython.bint keepConsensus=False,
                       cython.str reference="default",
                       pysam.cfaidx.FastaFile refHandle=None):
    cdef pPileupColumn_t PileupColumn
    cdef PCInfo_t PC
    cdef list DiscRPs, reads, discReads, DiscRPNames
    cdef pysam.calignmentfile.AlignedSegment read
    cdef pPileupRead_t i
    cdef VCFPos_t pos
    cdef int NumDiscordantPairs
    if(refHandle is None):
        raise ThisIsMadness("refHandle must be provided to write VCF lines!")
    PileupColumn = pPileupColumn(PileupCol)
    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ, experiment=experiment,
                minFracAgreed=minFracAgreed, minFA=minFA)
    pos = VCFPos(PC, MaxPValue=MaxPValue,
                 keepConsensus=keepConsensus,
                 reference=reference,
                 minFracAgreed=minFracAgreed,
                 minFA=minFA, refHandle=refHandle,
                 NDP=len(PC.DiscNames))
    return str(pos)


@cython.returns(tuple)
def pPileupColToVCFLines(pPileupColumn_t PileupColumn,
                         int minMQ=-1, int minBQ=-1,
                         cython.str experiment="",
                         int minFA=-1, cython.float minFracAgreed=-1.,
                         cython.float MaxPValue=-1.,
                         cython.bint keepConsensus=False,
                         cython.str reference="default",
                         pysam.cfaidx.FastaFile refHandle=None):
    cdef PCInfo_t PC
    cdef VCFPos_t pos
    if(refHandle is None):
        raise ThisIsMadness("refHandle must be provided to write VCF lines!")
    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ, experiment=experiment,
                minFracAgreed=minFracAgreed, minFA=minFA,
                experiment=experiment)
    """
    Commenting out this section - I moved the filtering out of
    discordant read pairs and into the PCInfo for some performance.
    PC.Records = [i for i in PC.Records if
                  i.alignment.query_name not in DiscRPNames]
    discReads = []
    for RP in DiscRPs:
        reads = RP.RP.getReads()
        for read in reads:
            read.set_tag("DP", RP.discordanceString, "Z")
            discReads.append(read)
    """
    pos = VCFPos(PC, MaxPValue=MaxPValue,
                 keepConsensus=keepConsensus,
                 reference=reference,
                 minFracAgreed=minFracAgreed,
                 minFA=minFA, refHandle=refHandle,
                 NDP=len(PC.DiscNames))
    return str(pos)


@cython.returns(tuple)
def IteratorColumnRegionToTuple(
        pysam.calignmentfile.IteratorColumnRegion ICR,
        int minMQ=-1, int minBQ=-1, cython.str experiment="",
        int minFA=-1, cython.float minFracAgreed=-1.,
        cython.float MaxPValue=-1., puEnd=-1,
        pysam.cfaidx.FastaFile refHandle=None,
        cython.bint keepConsensus=True,
        cython.str reference="default"):
    cdef pysam.calignmentfile.PileupColumn psPileupColumn
    cdef list VCFLines, discReads, allDiscReads
    if(puEnd < 0):
        raise ThisIsMadness("Need final entry in bed line to do ICR->VCF")
    VCFLines = []
    allDiscReads = []
    for psPileupColumn in ICR:
        posStr, discReads = PileupItToVCFLines(
            psPileupColumn, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
            minFracAgreed=minFracAgreed, experiment=experiment,
            MaxPValue=MaxPValue, refHandle=refHandle,
            keepConsensus=keepConsensus, reference=reference)
        if(len(posStr) != 0):
            VCFLines.append(posStr)
        allDiscReads = list(set(allDiscReads + discReads))
        if(psPileupColumn.pos > puEnd):
            break
    return VCFLines


@cython.returns(cython.str)
def IteratorColumnRegionToStr(
        pysam.calignmentfile.IteratorColumnRegion ICR,
        int minMQ=-1, int minBQ=-1, cython.str experiment="",
        int minFA=-1, cython.float minFracAgreed=-1.,
        cython.float MaxPValue=-1., puEnd=-1,
        pysam.cfaidx.FastaFile refHandle=None,
        cython.bint keepConsensus=True,
        cython.str reference="default"):
    cdef pysam.calignmentfile.PileupColumn psPileupColumn
    cdef list VCFLines, discReads, allDiscReads
    if(puEnd < 0):
        raise ThisIsMadness("Need final entry in bed line to do ICR->VCF")
    VCFLines = []
    allDiscReads = []
    for psPileupColumn in ICR:
        posStr, discReads = PileupItToVCFLines(
            psPileupColumn, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
            minFracAgreed=minFracAgreed, experiment=experiment,
            MaxPValue=MaxPValue, refHandle=refHandle,
            keepConsensus=keepConsensus, reference=reference)
        if(len(posStr) != 0):
            VCFLines.append(posStr)
        allDiscReads = list(set(allDiscReads + discReads))
        if(psPileupColumn.pos > puEnd):
            break
    return "\n".join(VCFLines)


@cython.returns(pysam.calignmentfile.IteratorColumnRegion)
def GetICR(line, pysam.calignmentfile.AlignmentFile inHandle=None):
    """
    Gets an ICR from an alignment file and a bed line.
    """
    try:
        return inHandle.pileup(line[0], line[1], line[2])
    except IndexError:
        try:
            return inHandle.pileup(line[0], line[1])
        except IndexError:
            raise ThisIsMadness("bed line has only one field?")


def GetICRString(bedline, bampath, minMQ=-1, minFA=-1, minBQ=-1,
                 experiment="", minFracAgreed=-1., MaxPValue=-1.,
                 ref="default", keepConsensus=True):
    if(ref == "default"):
        raise ThisIsMadness("Hey, I need a path to an faidx'd "
                            "fasta reference to variant-call.")
    puEnd = bedline[2]
    cStr = ("python -c 'import pysam; af = pysam.AlignmentFile("
            "\"%s\");  b = af.pileup" % bampath +
            "(\"%s\", %s, max_depth=200000," % (bedline[0], bedline[1]) +
            " multiple_iterators=False);"
            "from MawCluster.VCFWriters import IteratorColumnRegionToStr;" +
            "c = IteratorColumnRegionToStr(b, minMQ=%s, minBQ=" % (minMQ) +
            "%s, minFA=%s, experiment=\"%s\", minFr" % (minBQ, minFA,
                                                        experiment) +
            "acAgreed=%s, MaxPValue=%s, puEn" % (minFracAgreed, MaxPValue) +
            "d=%s, refHandle=pysam.FastaFile(\"%s\"), " % (puEnd, ref) +
            "keepConsensus=%s, reference=\"%s\"" % (keepConsensus, ref) +
            "); import sys; sys.stdout.write(c)'")
    return cStr


def GetPopenDispatcherICRs(bed, bampath, minMQ=-1, minFA=-1, minBQ=-1,
                           experiment="", minFracAgreed=-1., MaxPValue=-1.,
                           ref="default", keepConsensus=True, threads=2):
    if isinstance(bed, str):
        bed = ParseBed(bed)
    return PopenDispatcher([GetICRString(bedline, bampath, minMQ=minMQ,
                                         minFA=minFA, minBQ=minBQ,
                                         MaxPValue=MaxPValue, ref=ref,
                                         keepConsensus=keepConsensus,
                                         minFracAgreed=minFracAgreed) for
                            bedline in bed], threads=threads)


def PSNVCall(inBAM, conf="default", threads=-1, outVCF="default"):
        config = parseConfig(conf)
        if("threads" in config.keys() and threads < 0):
            threads = int(config["threads"])
        elif(threads < 0):
            threads = 4  # If threads is set in kwargs, keep it.
        if(outVCF == "default"):
            outVCF = TrimExt(inBAM) + ".psnv.bmf.vcf"
        outHandle = open(outVCF, "w")
        outHandle.write(GetVCFHeader(reference=config["ref"],
                        header=pysam.AlignmentFile(inBAM, "rb").header))
        pl("Splitting BAM file by contig.")
        Dispatcher = HTSUtils.GetBMFsnvPopen(inBAM, config['bed'],
                                             conf=conf,
                                             threads=threads)
        if(Dispatcher.daemon() != 0):
            raise ThisIsMadness("Dispatcher failed somehow.")
        pl("Shell calls completed without errors.")
        for commandString, vcffile in Dispatcher.outstrs.items():
            if(vcffile is None):
                raise FunctionCallException(
                    commandString,
                    "Attempt to cat this vcf to final file failed.", -1)
            check_call("cat %s | grep -v '^#' >> %s" % (vcffile, outVCF),
                       shell=True)
        pl("Filtering VCF by bed file. Pre-filter path: %s" % outVCF)
        bedFilteredVCF = BCVCF.FilterVCFFileByBed(
                    outVCF, bedfile=config["bed"])
        pl("Filtered VCF: %s" % bedFilteredVCF)
        return bedFilteredVCF
