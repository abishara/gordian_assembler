import os
import sys
import pysam
from itertools import chain, tee, izip
import cPickle as pickle
from collections import defaultdict, Counter, namedtuple

#--------------------------------------------------------------------------
# os
#--------------------------------------------------------------------------
def mkdir_p(dir):
  if not os.path.isdir(dir):
    os.makedirs(dir)

def mktmpdir(prefix='tmp'):
  path = tempfile.mkdtemp(dir='.', prefix=prefix)
  return path

class cd: 
  def __init__(self, newPath):
    self.newPath = newPath

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)

def concat_files(input_list, output_path):
  with open(output_path, 'w') as outf:
    for path in input_list:
      with open(path) as inf:
        for line in inf:
          outf.write(line)
  return

def touch(path, times=None):
  with open(path, 'a'):
    os.utime(path, times)

#--------------------------------------------------------------------------
# pickle
#--------------------------------------------------------------------------
def write_pickle(path, obj):
  f = open(path,'w')
  pickle.dump(
    obj,
    f,  
    pickle.HIGHEST_PROTOCOL
  )
  f.close()

def load_pickle(path):
  f = open(path,'r')
  obj = pickle.load(f)
  f.close()
  return obj

#--------------------------------------------------------------------------
# partitioning
#--------------------------------------------------------------------------
def pairwise(iterable):
  a, b = tee(iterable)
  next(b, None)
  return izip(a, b)
  
#--------------------------------------------------------------------------
# fastq
#--------------------------------------------------------------------------
def grouped(iterator, lines_per_read):
  while True:
    vals = tuple(next(iterator, None) for _ in xrange(lines_per_read))
    if None not in vals:
      yield vals
    else:
      raise StopIteration

def fastq_iter(f):
  for lines in grouped(f, 4):
    qname, bcode = lines[0].strip().split()
    bcode = bcode.split(':')[-1]
    yield bcode, qname, lines
  raise StopIteration

def get_fasta_sizes(fa_path):
  fasta = pysam.FastaFile(fa_path)
  ctg_size_map = {}
  for ctg in fasta.references:
    size = fasta.get_reference_length(ctg)
    ctg_size_map[ctg] = size
  return ctg_size_map

