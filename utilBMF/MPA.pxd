cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx

cimport utilBMF.HTSUtils

from cpython cimport array as c_array
from numpy cimport ndarray
from cython cimport bint
from utilBMF.Inliners cimport CigarOpToCigarChar, ChrToRefIDInline
from utilBMF.Inliners cimport chrInline, opLenToStr
from utilBMF.PysamUtils cimport PysamToChrInline
from utilBMF.cstring cimport PH2CHR_TRANS
from MawCluster.Math cimport CHI2_FROM_PHRED, INV_CHI2_FROM_PHRED, igamc
from libc.math cimport log10 as c_log10
from libc.stdint cimport int8_t

ctypedef cython.str cystr
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t
ctypedef c_array.array py_array

ctypedef PyLayout PyLayout_t
ctypedef LayoutPos LayoutPos_t
ctypedef utilBMF.HTSUtils.BamTag BamTag_t

cpdef cystr MPA2stdout(cystr inBAM)
cpdef MPA2Bam(cystr inBAM, cystr outBAM=?,
              bint u=?, bint coorsort=?,
              cystr sortMem=?, bint assume_sorted=?,
              bint dry_run=?,
              cystr prepend=?)
cpdef LayoutPos_t MergePositions(LayoutPos_t pos1, LayoutPos_t pos2)

cdef class LayoutPos:
    cdef public int pos, readPos, quality, agreement
    cdef public char operation, base, mergeAgreed
    cdef int8_t oqual
    cdef bint merged
    # cdef public cystr
    cpdef bint ismapped(self)
    cdef bint getMergeAgreed(self)
    cdef bint getMergeSet(self)
    cdef cystr tostring(self)


cdef class PyLayout:
    cdef public list positions
    cdef public dict tagDict
    cdef public int firstMapped, InitPos, flag, pnext, tlen, mapq
    cdef public cystr Name, contig, rnext
    cdef public bint isMerged, is_reverse, mergeAdjusted
    cdef int aend

    cdef py_array getQualStringScores(self)
    cpdef int getAlignmentStart(self)
    cpdef cystr getCigarString(self)
    cpdef cystr getSeq(self)
    cdef py_array getSeqArr(self)
    cdef int cGetRefPosForFirstPos(self)
    cpdef int getRefPosForFirstPos(self)
    cpdef py_array getAgreement(self)
    cdef py_array cGetAgreement(self)
    cdef py_array cGetQual(self)
    cpdef py_array getQual(self)
    cdef cystr cGetQualString(self)
    cpdef cystr getQualString(self)
    cdef int cGetLastRefPos(self)
    cpdef int getLastRefPos(self)
    cdef cystr cGetCigarString(self)
    cdef update_tags(self)
    cdef py_array cGetMergedPositions(self)
    cdef py_array cGetMergeAgreements(self)
    cpdef py_array getMergedPositions(self)
    cpdef py_array getMergeAgreements(self)
    cdef py_array cGetGenomicDiscordantPositions(self)
    cdef py_array cGetReadDiscordantPositions(self)
    cdef bint test_merge_success(self)
    cpdef py_array getReadDiscordantPositions(self)
    cpdef py_array getGenomicDiscordantPositions(self)

cpdef bint LayoutsOverlap(PyLayout_t L1, PyLayout_t L2)
cdef LayoutPos_t cMergePositions(LayoutPos_t pos1, LayoutPos_t pos2)
cdef int getLayoutLen(AlignedSegment_t read)

cdef class ListBool:
    cdef list List
    cdef bint Bool

cdef object oagtag


cdef inline int MergeAgreedQualities(int q1, int q2) nogil:
    return <int>(-10 * c_log10(igamc(2., CHI2_FROM_PHRED(q1) +
                                     CHI2_FROM_PHRED(q2) / 2.0)) + 0.5)


cdef inline int MergeDiscQualities(int q1, int q2) nogil:
    if(q1 > q2):
        return <int>(- 10 * c_log10(igamc(2., INV_CHI2_FROM_PHRED(q2) +
                                          CHI2_FROM_PHRED(q1))) + 0.5)
    else:
        return <int>(- 10 * c_log10(igamc(2., INV_CHI2_FROM_PHRED(q1) +
                                          CHI2_FROM_PHRED(q2))) + 0.5)

cdef inline char c_max(char i, char j) nogil:
    return i if(i > j) else j
