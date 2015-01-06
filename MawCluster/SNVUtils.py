import pysam

from MawCluster.PileupUtils import PCInfo, AlleleAggregateInfo
from utilBMF import HTSUtils
from utilBMF.HTSUtils import printlog as pl

"""
This module contains a variety of tools for calling variants.
Currently, it primarily works with SNPs primarily with experimental
features present for structural variants
TODO: Filter based on variants supported by reads going both ways.
TODO: Make calls for SNPs, not just reporting frequencies.
Reverse Strand Fraction - RSF
Both Strands Support Variant - BS
Fraction of unmerged reads supporting variant - TF
Allele Fraction - AF
DP - Depth (Merged)
"""


class VCFLine:

    def __init__(self,
                 AltAggregateObject,
                 MaxPValue=float("1e-15"),
                 ID=".",
                 DOCMerged="default",
                 DOCTotal="default",
                 TotalFracStr="default",
                 MergedFracStr="default",
                 TotalCountStr="default",
                 MergedCountStr="default",
                 FailedBQReads="default",
                 FailedMQReads="default"):
        if(isinstance(AltAggregateObject, AlleleAggregateInfo) is False):
            raise HTSUtils.ThisIsMadness("VCFLine requires an AlleleAgg"
                                         "regateInfo for initialization")
        if(DOCMerged == "default"):
            raise HTSUtils.ThisIsMadness("DOC (Merged) required!")
        if(DOCTotal == "default"):
            raise HTSUtils.ThisIsMadness("DOC (Total) required!")
        self.CHROM = AltAggregateObject.contig
        self.POS = AltAggregateObject.pos
        self.CONS = AltAggregateObject.consensus
        self.ALT = AltAggregateObject.ALT
        self.QUAL = AltAggregateObject.SumBQScore
        if(AltAggregateObject.BothStrandSupport is True):
            self.QUAL *= 100
            # This is terribly arbitrary...
        self.ID = ID
        try:
            if(float(MaxPValue) > 10 ** float(self.QUAL / -10)):
                self.FILTER = "PASS"
        except TypeError:
            print("TypeError! MaxPValue is: {}".format(MaxPValue))
            raise TypeError("MaxPValue not properly parsed.")
        else:
            self.FILTER = "LowQual"
        if(self.ALT == AltAggregateObject.consensus):
            self.FILTER = "CONSENSUS"
        self.InfoFields = {"AF": AltAggregateObject.MergedReads
                           / float(AltAggregateObject.DOC),
                           "TF": AltAggregateObject.TotalReads
                           / float(AltAggregateObject.DOCTotal),
                           "BS": AltAggregateObject.BothStrandSupport,
                           "RSF": AltAggregateObject.ReverseMergedReads
                           / AltAggregateObject.MergedReads,
                           "MQM": AltAggregateObject.AveMQ,
                           "MQB": AltAggregateObject.AveBQ,
                           "MMQ": AltAggregateObject.minMQ,
                           "MBQ": AltAggregateObject.minBQ,
                           "QA": AltAggregateObject.SumBQScore,
                           "NUMALT": AltAggregateObject.NUMALT,
                           "BQF": FailedBQReads,
                           "MQF": FailedMQReads,
                           "TYPE": "snp",
                           "MVQ": MaxPValue}
        if(TotalCountStr != "default"):
            self.InfoFields["TACS"] = TotalCountStr
        if(TotalFracStr != "default"):
            self.InfoFields["TAFS"] = TotalFracStr
        if(MergedCountStr != "default"):
            self.InfoFields["MACS"] = MergedCountStr
        if(MergedFracStr != "default"):
            self.InfoFields["MAFS"] = MergedFracStr
        self.InfoStr = ";".join(
            ["=".join([key, str(self.InfoFields[key])])
             for key in self.InfoFields.keys()])
        self.FormatFields = {"DP": DOCMerged,
                             "DPA": AltAggregateObject.MergedReads,
                             "DPT": DOCTotal,
                             "QA": AltAggregateObject.SumBQScore}
        self.FormatStr = (":".join(self.FormatFields.keys()) + "\t" +
                          ":".join(str(self.FormatFields[key])
                                   for key in self.FormatFields.keys()))
        self.str = "\t".join([str(i) for i in [self.CHROM,
                                               self.POS,
                                               self.ID,
                                               self.CONS,
                                               self.ALT,
                                               self.QUAL,
                                               self.FILTER,
                                               self.InfoStr,
                                               self.FormatStr]])

    def update(self):
        self.FormatKey = ":".join(self.FormatFields.keys())
        self.FormatValue = ":".join([str(self.FormatFields[key])
                                     for key in self.FormatFields.keys()])
        self.FormatStr = (":".join(self.FormatFields.keys()) + "\t" +
                          ":".join(str(self.FormatFields[key])
                                   for key in self.FormatFields.keys()))
        self.InfoStr = ";".join([key + "=" + str(self.InfoFields[key])
                                 for key in self.InfoFields.keys()])

    def ToString(self):
        self.update()
        self.str = "\t".join([str(i) for i in [self.CHROM,
                                               self.POS,
                                               self.ID,
                                               self.CONS,
                                               self.ALT,
                                               self.QUAL,
                                               self.FILTER,
                                               self.InfoStr,
                                               self.FormatStr]])
        return self.str


class VCFPos:

    def __init__(self, PCInfoObject,
                 MaxPValue=1e-15,
                 keepConsensus=True):
        if(isinstance(PCInfoObject, PCInfo) is False):
            raise HTSUtils.ThisIsMadness("VCFPos requires an "
                                         "PCInfo for initialization")
        self.pos = PCInfoObject.pos
        self.minMQ = PCInfoObject.minMQ
        self.consensus = PCInfoObject.consensus
        self.TotalFracStr = PCInfoObject.TotalFracStr
        self.MergedFracStr = PCInfoObject.MergedFracStr
        self.TotalCountStr = PCInfoObject.TotalCountStr
        self.MergedCountStr = PCInfoObject.MergedCountStr
        # Create VCFLines object using the calculated statistics
        # if(isinstance(MaxPValue, float) is False):
        #     print("repr of MaxPValue: {}".format(repr(MaxPValue)))
        #     raise HTSUtils.ThisIsMadness("MaxPValue must be a float!")
        self.VCFLines = [VCFLine(
            alt, TotalCountStr=self.TotalCountStr,
            MergedCountStr=self.MergedCountStr,
            TotalFracStr=self.TotalFracStr,
            MergedFracStr=self.MergedFracStr,
            DOCMerged=PCInfoObject.MergedReads,
            DOCTotal=PCInfoObject.TotalReads,
            MaxPValue=MaxPValue,
            FailedBQReads=PCInfoObject.FailedBQReads,
            FailedMQReads=PCInfoObject.FailedMQReads)
            for alt in PCInfoObject.AltAlleleData]
        self.keepConsensus = keepConsensus
        if(self.keepConsensus is True):
            self.str = "\n".join([line.ToString()
                                  for line in
                                  self.VCFLines])
        else:
            self.str = "\n".join([line.ToString()
                                  for line in
                                  self.VCFLines
                                  if line.FILTER != "CONSENSUS"])

    def ToString(self):
        [line.update() for line in self.VCFLines]
        if(self.keepConsensus is True):
            self.str = "\n".join([line.ToString()
                                  for line in
                                  self.VCFLines])
        else:
            self.str = "\n".join([line.ToString()
                                  for line in
                                  self.VCFLines
                                  if line.FILTER != "CONSENSUS"])
        return self.str


def SNVCrawler(inBAM,
               bed="default",
               minMQ=0,
               minBQ=0,
               VariantCallingTsv="default",
               MaxPValue=1e-15,
               keepConsensus=False
               ):
    if(isinstance(bed, str)):
        bed = HTSUtils.ParseBed(bed)
        if(VariantCallingTsv == "default"):
            VariantCallingTsv = inBAM[0:-4] + ".vc.tsv"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(VariantCallingTsv, "w")
    outHandle.write("##HEADER___NOT_WRITTEN")
    for line in bed:
        puIterator = inHandle.pileup(line[0], line[1], line[2],
                                     max_depth=30000)
        while True:
            try:
                PC = PCInfo(puIterator.next(), minMQ=minMQ, minBQ=minBQ)
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                pl("Finished iterations.")
                break
            VCFLineSet = VCFPos(PC, MaxPValue=MaxPValue,
                                keepConsensus=keepConsensus)
            # TODO: Check to see if it speeds up to not assign and only write.
            outHandle.write(VCFLineSet.ToString() + "\n")
    outHandle.close()
    inHandle.close()
    return None