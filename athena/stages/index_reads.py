import os

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

class IndexReadsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    yield IndexReadsStep(options)

  def outpaths(self, final=False):
    paths = {}
    paths['pass.file'] = os.path.join(self.outdir, 'pass')
    paths['index.file'] = FastqIndex.get_index_path(self.fq_path)
    #paths['shit.file'] = os.path.join(self.outdir, 'shit')
    return paths
 
  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      self.__class__.__name__,
      str(self),
    )

  def __init__(self, options):
    self.options = options
    self.fq_path = self.options.fq_path
    util.mkdir_p(self.outdir)

  def __fqid(self):
    return os.path.basename(os.path.dirname(os.path.dirname(self.fq_path)))

  def __str__(self):
    return '{}_{}'.format(
      self.__class__.__name__,
      self.__fqid(),
    )

  def run(self):
    #if self.fq_path.endswith('.gz'):
    #  cmd = 'zcat {} > {}'.format(self.fq_path, self.nfq_path)
    #  os.system(cmd)

    self.logger.log('index fastq {}'.format(self.fq_path))
    with FastqIndex(self.fq_path) as idx:
      pass
    passfile_path = os.path.join(self.outdir, 'pass')
    util.touch(passfile_path)
    self.logger.log('done')

