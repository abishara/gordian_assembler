import os
import sys
import pysam
import cPickle as pickle

import util

class BCBamIndex(object):

  file_suffix = '.bamidx.p'
  #file_suffix = '.bamidx2.p'

  @staticmethod
  def get_index_path(bam_path):
    return bam_path + BCBamIndex.file_suffix

  @property
  def bcode_set(self):
    if self._bcode_set == None:
      self._bcode_set = set(self._bcode_off_map.keys())
    return self._bcode_set

  def __init__(
    self,
    bam_path,
  ):
    
    self.bam_path = bam_path
    self.index_path = self.get_index_path(bam_path)
    self._bcode_set = None
    self._bcode_off_map = None

    if not os.path.isfile(self.index_path):
      self.__build_index__()
    else:
      self.__load_index__()

    self.f_map = None
    self.open()

  def open(self):
    assert self.f_map == None, "fp map already populated"
    self.f_map = {}
    self.f_map[self.bam_path] = pysam.Samfile(self.bam_path, 'rb')
    return self

  def close(self):
    for f in self.f_map.values():
      f.close()
    return

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    self.close()

  def __build_index__(self):  
    self._bcode_off_map = {}
    seen_set = set()
    step = 50000
    fhandle = pysam.Samfile(self.bam_path, 'rb')
    print 'building index for bam', self.bam_path
    fp = None
    pfp = None
    for read in fhandle:
      bcode = get_barcode(read)
      if bcode == None:
        continue
      pfp = fp
      fp = fhandle.tell()
      seen_set.add(bcode)
      if bcode not in self._bcode_off_map:
        self._bcode_off_map[bcode] = pfp
        if len(seen_set) % step == 0:
          print '  - indexed {} barcodes'.format(len(seen_set))
    fhandle.close()

    print 'writing index for bams'
    for bam_path in [self.bam_path]:
      print '  -', bam_path
    util.write_pickle(self.index_path, self._bcode_off_map)

  def __load_index__(self):  
    print 'loading index', self.index_path
    self._bcode_off_map = util.load_pickle(self.index_path)
    print '  - done'

  def get_reads(self, bcode):
    assert self._bcode_off_map != None, 'index {} not loaded'.format(self.index_path)

    if bcode not in self.bcode_set:
      raise StopIteration

    offset = self._bcode_off_map[bcode]
    fhandle = self.f_map[self.bam_path]
    fhandle.seek(offset, 0)

    for read in fhandle:
      r_bcode = get_barcode(read)
      if r_bcode != bcode:
        break
      yield read

    raise StopIteration

def get_barcode(read):
  filt_list = filter(lambda(k, v): k == 'BX', read.tags)
  if filt_list == []: 
    return None
  else:
    k, v = filt_list[0]
    return v
