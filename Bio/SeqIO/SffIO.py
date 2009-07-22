# Copyright 2009 by James Casbon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the Standard Flowgram Format (SFF) produced by 454 Sequencers.

You are expected to use this module via the Bio.SeqIO functions.

This implements the SFF file format as described in the "Genome Sequencer 
Data Analysis Software Manual Software Version 2.0.00, October 2008"" pages 528-531.

The file is split into three sections: 

1. Common Header
2. Reads sections which has pairs of:
    Read Header
    Read Data
3. Index

This parser implements the first two sections. 

Note the 454 analysis software will output fna and qual files which can be read by the QualityIO 
package.  This may be useful to you if:
 * you want to process an sff file on a system without the Roche software
 * you want to access the untrimmed sequences 
 * you want to look at the actual flow data
"""

import struct

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord

class SffException(Exception):
    pass


class SffHeader(object):

    header = [
        '>', # use big endian
        'I', # Magic number, should be 0x2E736666 (779314790)
        '4s',# version number 
        'Q', # index offset
        'I', # index length
        'I', # number of reads
        'H', # header length
        'H', # key length
        'H', # number of flows per read
        'B', # flowgram format code (1)
        '400s', # flows
        '4s', # key sequence
    ]

    header_fmt = ''.join(header)
    
    # add in eight byte padding
    while struct.calcsize(header_fmt) % 8 != 0:
        header_fmt = header_fmt + 'x'
        
    size = struct.calcsize(header_fmt)

    def __init__(self, stream):
        data = struct.unpack(self.header_fmt, stream.read(self.size))
        self.raw = data
        self.magic_number = data[0]
        self.version = data[1]
        self.index_offset = data[2]
        self.index_length = data[3]
        self.number_of_reads = data[4]
        self.header_length = data[5]
        self.key_length = data[6]
        self.number_of_flows = data[7]
        self.flowgram_format_code = data[8]
        self.flows = data[9]
        self.key_sequence = data[10]
        
        try: 
            assert self.magic_number == 779314790, 'wrong magic number'
            assert self.version == '\x00\x00\x00\x01', 'wrong version'
            assert self.flowgram_format_code == 1, 'wrong flowgram format'
        except AssertionError, msg:
            raise SffException(msg)
        

class SffReadHeader(object):
    format = [
        '>'
        'H', # read header length
        'H', # name length
        'I', # number of bases
        'H', # clip qual left
        'H', # clip qual right
        'H', # clip adaptor left
        'H', # clip adaptor right 
        '14s' # accession
    ]
    format = ''.join(format)
    
    # add in eight byte padding
    while struct.calcsize(format) % 8 != 0:
        format = format + 'x'
        
    size = struct.calcsize(format)
    
    def __init__(self, stream):
        data = struct.unpack(self.format, stream.read(self.size))
        self.raw = data
        self.name_length = data[1]
        self.number_of_bases = data[2]
        self.clip_qual_left = data[3] - 1 if data[3] else 0
        self.clip_qual_right = data[4] - 1 if data[4] else self.number_of_bases
        self.clip_adaptor_left = data[5] -1 if data[5] else 0
        self.clip_adaptor_right = data[6] -1 if data[6] else self.number_of_bases
        self.accession = data[7]
        
        try: 
            assert self.name_length == 14, 'read name length is not 14'
        except AssertionError, msg:
            raise SffException(msg) 
            
        

class SffReadData(object):
    
    def __init__(self, number_of_flows, header, stream):
        self.header = header
        number_of_bases = header.number_of_bases
        
        fv_format = '>%sH' % number_of_flows 
        self.flow_values = struct.unpack(fv_format, stream.read(struct.calcsize(fv_format)))
        
        fi_format = '>%sB' % number_of_bases
        self.flow_index = struct.unpack(fi_format, stream.read(struct.calcsize(fi_format)))
        
        b_format = '>%ss' % number_of_bases
        self.bases = struct.unpack(b_format, stream.read(struct.calcsize(b_format)))[0]
        
        qs_format = '>%sB' % number_of_bases
        # add in eight byte padding
        while struct.calcsize(fv_format[1:]+fi_format[1:]+b_format[1:]+qs_format[1:]) % 8 != 0:
            qs_format = qs_format + 'x'
        self.quality_scores = struct.unpack(qs_format, stream.read(struct.calcsize(qs_format)))
        
    
    @property
    def start(self):
        return max(self.header.clip_adaptor_left, self.header.clip_qual_left)
        
    @property
    def end(self):
        return min(self.header.clip_adaptor_right, self.header.clip_qual_right)
        
    @property
    def sequence(self):
        return self.bases[self.start:self.end]
        
    @property
    def qualities(self):
        return self.quality_scores[self.start:self.end]

        
class SffFile(object):
    """ SffFile parser. Create with a stream and iterate to get pairs of (ReadHeader, ReadData)
    """
    
    def __init__(self, stream):
        self.stream = stream 
        self.header = SffHeader(self.stream)
        
    def __iter__(self):
        reads = 0
        while True:
            reads += 1 
            rh = SffReadHeader(self.stream)
            rd = SffReadData(self.header.number_of_flows, rh, self.stream)
            yield (rh, rd)
            
            if reads == self.header.number_of_reads:
                # breaking here to prevent reading into the index
                break


def SffIterator(handle, alphabet = single_letter_alphabet, title2ids = None):
    sfffile = SffFile(handle)
    
    for (rh, rd) in sfffile:
        record = SeqRecord(Seq(rd.sequence, alphabet),
                           id=rh.accession,
                           name=rh.accession,
                           description=rh.accession)
        record.letter_annotations["phred_quality"] = rd.qualities
        yield record
