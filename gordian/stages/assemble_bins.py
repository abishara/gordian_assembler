import os
import sys
import subprocess
from collections import defaultdict
import shutil
import pysam

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

# NOTE must be in path
idbabin_path = 'idba_ud'

wd = os.path.dirname(os.path.abspath(__file__))

class AssembleBinsStep(StepChunk):
  @staticmethod
  def get_steps(options):
    bins_map = util.load_pickle(options.phased_bins_path)
    for idx, bcodes in bins_map.items():
      yield AssembleBinsStep(options, idx, bcodes)
      #break

  def outpaths(self, final=False):
    paths = {}
    paths['contig.fa'] = os.path.join(self.outdir, 'contig.fa')
    paths['binreads.fa'] = os.path.join(self.outdir, 'binreads.fa')
    #paths['shit.file'] = os.path.join(self.outdir, 'shit')
    return paths
 
  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      self.__class__.__name__,
      str(self),
    )
 
  def __str__(self):
    return 'c{}'.format(self.idx)

  def __init__(self, options, idx, bcodes):
    self.options = options
    self.idx = idx
    self.bcodes = bcodes
    util.mkdir_p(self.outdir)

  def run(self):
    scratch_path = self.options.get_bin_scratch(self.idx)
    util.mkdir_p(scratch_path)
    readsfa_path = os.path.join(scratch_path, 'binreads.fa')
    asmdir_path = os.path.join(scratch_path, 'asm')
    contig_path = os.path.join(asmdir_path, 'contig.fa')

    # get all barcoded reads 
    seen_set = set()
    with open(readsfa_path, 'w') as fout:
      with FastqIndex(self.options.fq_path) as idx:
        for bcode in self.bcodes:
          for bcode, qname, lines in idx.get_reads(bcode):
            seq = lines[1].strip()
            fout.write('>{}\t{}\n'.format(qname, bcode))
            fout.write('{}\n'.format(seq))
            seen_set.add(bcode)
      assert self.bcodes.issubset(seen_set), 'not all barcodes loaded'

    # run idba
    cmd = '{} --num_threads 2 -r {} -o {}'.format(
      idbabin_path,
      readsfa_path,
      asmdir_path,
    )
    pp = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    retcode = pp.wait()
    print "::::", retcode
    if not os.path.exists(contig_path):
      if retcode != 0 and  "invalid insert distance" in pp.stderr.read():
        self.logger.log("idba_ud internal error; this is probably a bug in idba_ud, so we have to bail here")
        open(contig_path, "w")
      else:
        error_message = "assembly failed to produce contig.fa"
        self.logger.log(error_message)
        raise Exception()
    elif retcode != 0:
      self.logger.log(
        "something funky happened while running idba_ud (got error code "
        "{}) but there's not much we can do so continuing".format(retcode))
    # check contig.fa is valid
    try:
      fasta = pysam.FastaFile(contig_path)
    except Exception as e:
      self.logger.log('failed to produce correctly formatted contigs in {}'.format(contig_path))
      sys.exit(1)

    shutil.copy(contig_path, os.path.join(self.outdir))
    shutil.copy(readsfa_path, os.path.join(self.outdir))

