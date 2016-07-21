import os
import pysam
import shutil
from collections import defaultdict, Counter
import networkx as nx

from .step import StepChunk
from ..mlib import util

import haplotype_reads
from compute_overlaps import ComputeOverlapsStep

MIN_BARCODE_OVERLAP = 20

class GroupBinsStep(StepChunk):

  @classmethod
  def get_steps(cls, options):
    yield cls(options)

  def outpaths(self, final=False):
    paths = {}
    #paths['shit'] = 'shit'
    paths['groups.p'] = self.options.groups_pickle_path
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
    self.logger.log('read in computed overlaps and bin overlap graph')
    G = nx.Graph()
    binid_bcodes_map = {}
    for binid, bcode_set in ComputeOverlapsStep.get_all_bins(self.options):
      G.add_node(binid)
      binid_bcodes_map[binid] = bcode_set

    for binid_olaps_map in ComputeOverlapsStep.get_all_overlaps(self.options):
      for binid, olaps in binid_olaps_map.items():
        for cbinid, olap in olaps:
          G.add_edge(binid, cbinid, olap=olap)

    self.logger.log('  - done')
    self.logger.log('  - {} bin nodes and {} overlap edges'.format(
      G.number_of_nodes(),
      G.number_of_edges(),
    ))

    def get_bcode(read):
      filt_list = filter(lambda(k, v): k == 'BX', read.tags)
      if filt_list == []: 
        return None
      else:
        k, v = filt_list[0]
        return v

    def is_overlap_edge(e):
      bid1, bid2 = e
      ctg1, b1, e1 = bid1[:-1]
      ctg2, b2, e2 = bid2[:-1]
      if ctg1 != ctg2:
        return False
      return max(b1, b2) <= min(e1, e2)

    # iterate through nodes and remove edges to the same window (but
    # different bins) if there is enough of a barcode overlap difference

    def get_window(bid):
      return bid[:-1]

    numdup = 0
    remove_edges = set()
    self.logger.log('remove weak barcode overlap edges')
    for bid1 in G.nodes():
      window_counts = Counter(map(get_window, G.neighbors(bid1)))
      dup_windows = set(filter(
        lambda(w): window_counts[w] > 1,
        window_counts.keys(),
      ))
      if len(dup_windows) == 0:
        continue

      window_n_map = defaultdict(list)
      for bid2 in G.neighbors(bid1):
        window = get_window(bid2)
        if window in dup_windows:
          window_n_map[window].append(bid2)

      for window in dup_windows:
        olaps = sorted(map(
          lambda(bid2): (G[bid1][bid2]['olap'], bid2),
          window_n_map[window],
        ), reverse=True)
        if olaps[0][0] > 2 * olaps[1][0]:
          # remove other barcode overlap edges
          for _, bid2 in olaps[1:]:
            remove_edges.add((bid1, bid2))

    G.remove_edges_from(remove_edges)
    self.logger.log('  - {} bin nodes and {} overlap edges'.format(
      G.number_of_nodes(),
      G.number_of_edges(),
    ))

    def get_bin_str(bid):
      ctg, b, e, cidx = bid
      return '{}.{}-{}.c{}'.format(ctg, b, e, cidx)

    # create subgraph with only nonduplicate edges
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    break_nodes = set()
    dup_edges = set()

    seen_set = set()
    self.logger.log('create subgraph with nonduplicate window edge connections')
    for bid1 in G.nodes():
      n_window_map = dict(map(
        lambda(bid2): (bid2, get_window(bid2)),
        G.neighbors(bid1),
      ))
      window_counts = Counter(map(get_window, G.neighbors(bid1)))
      dup_windows = set(filter(
        lambda(w): window_counts[w] > 1,
        window_counts.keys(),
      ))
      if len(window_counts) > 2 or len(dup_windows) > 0:
        break_nodes.add(bid1)
      for bid2, w in n_window_map.items():
        olap = G[bid1][bid2]['olap']
        if w not in dup_windows:
          H.add_edge(bid1, bid2, olap=olap)
        else:
          dup_edges.add((bid1, bid2))
          self.logger.log('  dup edge {} to {}, olap {}'.format(
            get_bin_str(bid1),
            get_bin_str(bid2),
            olap,
          ))
    H.remove_edges_from(dup_edges)
    for n in H.nodes():
      if len(H.neighbors(n)) == 0:
        print 'singleton node', n

    ccs = nx.connected_component_subgraphs(H)
    seen_set = set()
    groups = []
    for gid, cc in enumerate(sorted(ccs,key=len,reverse=True)):
      group = set(cc.nodes())
      bcodes = set()
      for bid in group:
        bcodes |= binid_bcodes_map[bid]

      assert len(seen_set & group) == 0, "node duplicated in subgraph"
      seen_set |= group
      self.logger.log("group {}, num bins {}".format(gid, len(group)))
      for bid in sorted(group):
        self.logger.log("  - node {}".format(str(bid)))

      groups.append((gid, group, bcodes))

    util.write_pickle(self.options.groups_pickle_path, groups)

    self.logger.log('done')

