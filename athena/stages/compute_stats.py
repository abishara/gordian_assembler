import os
import pysam
import shutil
import subprocess
from collections import defaultdict, Counter
from bx.intervals.cluster import ClusterTree

from .step import StepChunk
from ..mlib import util

from assemble_groups import AssembleGroupsStep

class ComputeStatsStep(StepChunk):

  @classmethod
  def get_steps(cls, options):
    yield cls(options)

  def outpaths(self, final=False):
    paths = {}
    paths['shit'] = 'shit'
    return paths
 
  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      self.__class__.__name__,
    )

  def __str__(self):
    return self.__class__.__name__

  def __init__(
    self,
    options,
  ):
    self.options = options
    util.mkdir_p(self.outdir)

  def run(self):
    # load all computed overlaps and build graph
    self.logger.log('merge aligned contigs in each bin')
    group_steps = list(AssembleGroupsStep.get_steps(self.options))
    inputbam_paths = map(
      lambda(s): s.outpaths()['align.contig.bam'],
      group_steps,
    )
    inputfa_paths = map(
      lambda(s): s.outpaths()['contig.fa'],
      group_steps,
    )
    inputpaths_str = ' '.join(inputbam_paths)
    mergedbam_path = os.path.join(self.outdir, 'merged-contigs.bam')
    cmd = 'samtools merge {} {}'.format(
      mergedbam_path,
      inputpaths_str,
    )
    #print 'cmd', cmd
    #subprocess.check_call(cmd, shell=True)
    cmd = 'samtools index {}'.format(mergedbam_path)
    #print 'cmd', cmd
    #subprocess.check_call(cmd, shell=True)

    ctg_size_map = {}
    for fa_path in inputfa_paths:
      # FIXME remove
      break
      ctg_size_map.update(util.get_fasta_sizes(fa_path))

    with util.cd(self.outdir):
      #self.logger.log('filter for covering contigs')
      #filter_contigs(
      #  ctg_size_map,
      #  mergedbam_path,
      #  'athena-filt.bam',
      #)
      #cmd = 'samtools index athena-filt.bam'
      #subprocess.check_call(cmd, shell=True)

      #self.logger.log('get athena coverage stats for roi')
      #get_coverage_stats(
      #  'athena-filt.bam',
      #  self.options.regions_bed_path,
      #  'athena-hits.bed',
      #)

      #self.logger.log('filter for covering contigs in allpaths')
      #allpathsfa_path = '/scratch/PI/serafim/abishara/data/na12878/allpaths-lg/contigs.fa'
      #allpathsbam_path = '/scratch/PI/serafim/abishara/data/na12878/allpaths-lg/align-contig.bam'
      #ctg_size_map = util.get_fasta_sizes(allpathsfa_path)
      #filter_contigs2(
      #  self.options.regions_bed_path,
      #  ctg_size_map,
      #  allpathsbam_path,
      #  'allpaths-filt.bam',
      #)
      #cmd = 'samtools index allpaths-filt.bam'
      #subprocess.check_call(cmd, shell=True)

      #self.logger.log('get allpaths coverage stats for roi')
      #get_coverage_stats(
      #  'allpaths-filt.bam',
      #  self.options.regions_bed_path,
      #  'allpaths-hits.bed',
      #)

      self.logger.log('filter for covering contigs in supernova')
      supernovafa_path  = '/scratch/PI/serafim/abishara/data/na12878/supernova/NA12878_megabubbles.fasta'
      supernovabam_path = '/scratch/PI/serafim/abishara/data/na12878/supernova/align-contig.bam'
      ctg_size_map = util.get_fasta_sizes(supernovafa_path)
      filter_contigs2(
        self.options.regions_bed_path,
        ctg_size_map,
        supernovabam_path,
        'supernova-filt.bam',
      )
      cmd = 'samtools index supernova-filt.bam'
      subprocess.check_call(cmd, shell=True)

      self.logger.log('get supernova coverage stats for roi')
      get_coverage_stats(
        'supernova-filt.bam',
        self.options.regions_bed_path,
        'supernova-hits.bed',
      )

    self.logger.log('done')

def filter_contigs2(
  roi_bed_path,
  size_map,
  inbam_path,
  outbam_path,
):

  fin = pysam.Samfile(inbam_path, 'rb')
  fout = pysam.Samfile(outbam_path, 'wb', template=fin)

  # filter contigs where at least 95% of the query sequence maps
  seen_set = set()
  seen2_set = set()
  for (ctg, b, e) in iter_bed(roi_bed_path):
    for read in fin.fetch(ctg, b, e):
      if read.is_unmapped or read.is_secondary:
        continue
      #if read.query_alignment_length < 0.95 * size_map[read.qname]:
      #  continue
      #if read.qname in seen_set:
      #  continue
      if (read.qname, read.pos) in seen2_set:
        continue
      fout.write(read)
      seen_set.add(read.qname)
      seen2_set.add((read.qname, read.pos))

  fin.close()
  fout.close()

def filter_contigs(size_map, inbam_path, outbam_path):

  fin = pysam.Samfile(inbam_path, 'rb')
  fout = pysam.Samfile(outbam_path, 'wb', template=fin)

  # filter contigs where at least 95% of the query sequence maps
  for read in fin:
    if read.is_unmapped:
      continue
    if read.query_alignment_length < 0.95 * size_map[read.qname]:
      continue
    fout.write(read)

  fin.close()
  fout.close()

def get_coverage_stats(
  inbam_path,
  roi_bed_path,
  roi_hits_bed_path,
):

  fin = pysam.Samfile(inbam_path, 'rb')

  with open(roi_hits_bed_path, 'w') as fout:
    for ctg, b, e in iter_bed(roi_bed_path):
      regions = ClusterTree(1,1)
      for i, read in enumerate(fin.fetch(ctg, b, e)):
        tpos = max(b, read.pos)
        taend = min(e, read.aend)
        regions.insert(tpos, taend, i)
      for ob, oe, _ in regions.getregions():
        fout.write('{}\t{}\t{}\n'.format(ctg, ob, oe))

  fin.close()

def iter_bed(bed_path):
  with open(bed_path) as fin:
    for line in fin:
      if line.startswith('#'):
        continue
      words = line.split('\t')
      ctg = words[0]
      b = int(words[1])
      e = int(words[2])
      yield (ctg, b, e)

  raise StopIteration

