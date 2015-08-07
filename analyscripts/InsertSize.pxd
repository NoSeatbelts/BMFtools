from pysam.calignmentfile cimport AlignmentFile, IteratorColumnRegion
from pysam.calignedsegment cimport PileupColumn, PileupRead
from libc.stdint cimport *
from cpython cimport array as c_array
from cython cimport bint
cimport cython
cimport numpy as np
from numpy cimport ndarray

from utilBMF.HTSUtils cimport pPileupRead
from utilBMF.Inliners cimport Nuc2NumN

# Typedefs
# Cython
ctypedef cython.str cystr
ctypedef c_array.array py_array_t

# Numerical definitions
ctypedef double double_t

# Pysam
ctypedef AlignmentFile AlignmentFile_t
ctypedef PileupColumn PileupColumn_t
ctypedef IteratorColumnRegion IteratorColumnRegion_t
ctypedef PileupRead PileupRead_t

# BMFTools
ctypedef pPileupRead pPileupRead_t


# Class definitions
cdef class SNVAlleleWrangler:
    cdef AlignmentFile_t handle
    cdef public int32_t reference_id, pos, max_insert_size
    cdef public list alleles
    cdef public py_array_t insert_sizes
    cdef public dict insert_size_dict
    cdef public cystr ref
    cdef py_array_t filter_cnf

    cdef inline bint pass_record(self, pPileupRead_t PR)
    cdef inline void fast_forward(self)
    cdef inline void build_insert_size_dict(self)
    cdef inline void get_insert_sizes(self)
    cdef ndarray[int32_t, ndim=2] get_allele_counts(self, dict insert_size_dict)

cdef inline int iabs(int integer) nogil:
    return integer if(integer > 0) else -1 * integer