import os
import sys
import subprocess
from collections import defaultdict
import shutil
from phaser.phase import phase_barcodes

from .step import StepChunk
from ..mlib import util

class PhaseBarcodesStep(StepChunk):
  @staticmethod
  def get_steps(options):
    yield PhaseBarcodesStep(options)

  def outpaths(self, final=False):
    paths = {}
    paths['phased_bins.p'] = self.options.phased_bins_path
    #paths['shit.file'] = os.path.join(self.outdir, 'shit')
    return paths
 
  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      str(self),
    )
 
  def __str__(self):
    return self.__class__.__name__

  def __init__(self, options):
    self.options = options

  def run(self):
    scratch_path = os.path.join(
      self.options.working_dir,
      str(self),
    )
    util.mkdir_p(scratch_path)

    bcode_whitelist = None
    if 'bcode_whitelist_path' in self.options.phasing_args:
      try:
        whitelist_path = self.options.phasing_args['bcode_whitelist_path']
        bcode_whitelist = util.load_pickle(whitelist_path)
        self.logger.log("using barcode whitelist of {} barcodes".format(
          len(bcode_whitelist)))
      except:
        self.logger.error('malformed barcode whitelist pickle file {}'.format(whitelist_path))
        sys.exit(1)

    K = self.options.phasing_args['K']
    self.logger.log("phasing barcodes for K = {} haps".format(K))
    clusters_map = phase_barcodes(
      self.options.bam_path,
      self.options.vcf_path,
      scratch_path,
      K=K,
      bcodes=bcode_whitelist,
    )
    self.logger.log('  - done')
    util.write_pickle(self.options.phased_bins_path, clusters_map)
      
