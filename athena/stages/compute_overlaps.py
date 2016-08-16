import os
import pysam
import shutil
from collections import defaultdict

from .step import StepChunk
from ..mlib import util
import haplotype_reads

MIN_BARCODE_OVERLAP = 20

class ComputeOverlapsStep(haplotype_reads.HaplotypeReadsStep):

  def run(self):
    self.logger.log('determine all vs all bin overlaps')

    binid_bcodes_map = dict(map(
      lambda(k, v1, v2): (k, v1),
      self.get_bins(),
    ) )
    binid_olaps_map = defaultdict(list)
    for cbinid, cbcodes, _ in self.get_all_bins(self.options):
      for binid, bcodes in binid_bcodes_map.items():
        # only compute each comparison once
        if cbinid <= binid:
          continue

        olap = len(bcodes & cbcodes)
        if olap >= MIN_BARCODE_OVERLAP:
          binid_olaps_map[binid].append((cbinid, olap))

    stats_path = self.outpaths()['stats']
    util.write_pickle(stats_path, binid_olaps_map)
    self.logger.log('done')

  def get_bins(self):
    hap_step = haplotype_reads.HaplotypeReadsStep(
      self.options,
      self.ctg,
      self.begin,
      self.end,
    )
    hapstats_path = hap_step.outpaths()['stats']
    assert os.path.isfile(hapstats_path), "haplotyper not yet run on window {}".format(str(self))
    (cluster_info_map, _) = util.load_pickle(hapstats_path)
    unassn_bcodes = cluster_info_map[None][1]
    for cidx, info in cluster_info_map.items():
      numreads, bcode_set, _, assm = info
      # add in unassigned
      all_bcode_set = unassn_bcodes | bcode_set
      binid = (self.ctg, self.begin, self.end, cidx)
      if assm and cidx != None:
        yield (binid, bcode_set, all_bcode_set)

  @classmethod
  def get_all_bins(cls, options):
    steps = list(cls.get_steps(options))
    mstep = max(1, len(steps) / 100)
    print 'iterating bins from {} windows'.format(len(steps))
    for i, step in enumerate(steps):
      if i % mstep == 0:
        print '  - {}/ {}'.format(i, len(steps))
      for binid, bcodes, all_bcodes in step.get_bins():
        yield binid, bcodes, all_bcodes
        
  def get_overlaps(self):
    stats_path = self.outpaths()['stats']
    assert os.path.isfile(stats_path), "overlaps not yet run on window {}".format(str(self))
    binid_olaps_map = util.load_pickle(stats_path)
    return binid_olaps_map

  @classmethod
  def get_all_overlaps(cls, options):
    steps = list(cls.get_steps(options))
    mstep = max(1, len(steps) / 100)
    print 'iterating overlaps from {} windows'.format(len(steps))
    for i, step in enumerate(steps):
      if i % mstep == 0:
        print '  - {}/ {}'.format(i, len(steps))
      yield step.get_overlaps()



