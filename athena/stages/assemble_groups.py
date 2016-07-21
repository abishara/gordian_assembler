import os
import pysam
import subprocess
from collections import defaultdict
import glob
import shutil

from ..assembler_tools.architect import architect
from ..assembler_tools.architect import local_assembly

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

# NOTE must be in path
idbabin_path = 'idba_ud'

wd = os.path.dirname(os.path.abspath(__file__))
architect_scripts_path = os.path.join(
  wd,
  '..',
  'assembler_tools/architect/scripts/',
)
edgesbin_path = os.path.join(architect_scripts_path, 'pe-connections.py')
containmentbin_path = os.path.join(architect_scripts_path, 'bam_to_containment.py')

#--------------------------------------------------------------------------
# base class
#--------------------------------------------------------------------------
class AssembleBaseStep(StepChunk):

  # NOTE to be defined in subclass
  #@staticmethod
  #def get_steps(options):
  #  pass
  #def run(self):
  #def __init__(self):

  @property
  def outdir(self):
    return self.options.get_group_dir(self.gid, final=True)

  def outpaths(self, final=False):
    paths = {}
    paths['contig.fa'] = os.path.join(self.outdir, 'contig.fa')
    #paths['local-asm-merged.fa'] = os.path.join(self.outdir, 'local-asm-merged.fa')
    paths['architect-scaff.fasta'] = os.path.join(self.outdir, 'scaff.fasta')
    return paths

  def clean_working(self):
    bin_path = self.options.get_group_dir(self.gid)
    self.logger.log('removing bin directory {}'.format(bin_path))
    shutil.rmtree(bin_path)
    return 

  def assemble(self):
    self.logger.log('assembling barcoded reads for this bin')

    fqdir_path = self.options.get_group_fq_dir(self.gid)
    allfa_path = os.path.join(fqdir_path, 'tenxreads-bcoded.fa')
    asmdir_path = self.options.get_group_asm_dir(self.gid)
    # initial idba assembly
    cmd = '{} -r {} -o {}'.format(
      idbabin_path,
      allfa_path,
      asmdir_path,
    )
    subprocess.check_call(cmd, shell=True)
    
    with util.cd(asmdir_path):
      assert os.path.isfile('contig.fa'), 'idba failed to generate contig.fa'

      # index initial contigs
      cmd = 'bwa index contig.fa'
      subprocess.check_call(cmd, shell=True)

      # align reads to contigs
      self.logger.log('aligning reads to contigs')
      cmd = 'bwa mem -p contig.fa {} > align.on-contig.sam'.format(
        '../fqs/tenxreads-bcoded.fa'
      )
      subprocess.check_call(cmd, shell=True)
      cmd = 'cat align.on-contig.sam | samtools view -bS - | samtools sort - align.on-contig.sorted'
      subprocess.check_call(cmd, shell=True)

      # filter bad reads, create idba *fa reads with only filtered reads
      allbam_path = 'align.on-contig.sorted.bam'
      filtbam_path = 'align.on-contig.sorted.filt.bam'
      filter_alignments(
        allbam_path,
        filtbam_path,
      )
      cmd = 'samtools index {}'.format(filtbam_path)
      subprocess.check_call(cmd, shell=True)

      # create edges.tsv and containment files from architect
      self.logger.log('generating local reassembly inputs')
      cmd = 'python {} -b {} -f {} -e {}'.format(
        edgesbin_path,
        filtbam_path,
        'contig.fa',
        'edges.tsv',
      )
      subprocess.check_call(cmd, shell=True)
      cmd = 'python {} -b {} -c {}'.format(
        containmentbin_path,
        filtbam_path,
        'containment',
      )
      subprocess.check_call(cmd, shell=True)

      # extract barcoded filtered alignments
      filtfa_path = 'tenxreads-bcoded-filt.fa'
      get_filtered_fa(
        filtbam_path,
        '../fqs/tenxreads-bcoded.fa',
        filtfa_path,
      )

      # run architect for local assemblies
      #local_assembly.setup_cwd()
      #self.logger.log('performing local reassembly')
      #architect.do_local_reassembly(
      #  'contig.fa',
      #  'edges.tsv',
      #  'containment',
      #  filtfa_path,
      #)

      architect.do_scaffolding(
        'contig.fa',
        'edges.tsv',
        'containment',
        0,
        'scaff',
      )

    def copyfinal(src, dest):
      shutil.copyfile(
        os.path.join(self.options.get_group_asm_dir(self.gid), src),
        os.path.join(self.options.get_group_dir(self.gid, final=True), dest),
      )

    self.logger.log('copying deliverables to final')
    copyfinal('contig.fa', 'contig.fa')
    #copyfinal('local-assemblies/local-asm-merged.fa', 'local-asm-merged.fa')
    copyfinal('scaff.fasta', 'scaff.fasta')

    self.logger.log('done')

#--------------------------------------------------------------------------
# assemble specified barcoded reads step
#--------------------------------------------------------------------------
class AssembleSpecReadsStep(AssembleBaseStep):

  @staticmethod
  def get_steps(options):
    yield AssembleSpecReadsStep(options)

  def __init__(
    self,
    options,
  ):
    self.gid = None
    self.options = options
    fqdir_path = self.options.get_group_fq_dir(self.gid)
    util.mkdir_p(self.outdir)
    util.mkdir_p(fqdir_path)

  def __str__(self):
    return self.__class__.__name__

  def run(self):
    self.logger.log('copy input fq to scratch directory')
    fqdir_path = self.options.get_group_fq_dir(self.gid)
    shutil.copy(
      self.options.tenxfq_path, 
      os.path.join(fqdir_path, 'tenxreads.fq'),
    )
    self.assemble()

#--------------------------------------------------------------------------
# assemble barcoded reads binned by buckets and haplotyping of all reads
#--------------------------------------------------------------------------
class AssembleGroupsStep(AssembleBaseStep):

  @classmethod
  def get_steps(cls, options):
    groups = util.load_pickle(options.groups_pickle_path)

    for i, (gid, bins, bcodes) in enumerate(groups):
      yield cls(options, gid, bcodes)

  def __init__(
    self,
    options,
    gid,
    bcodes,
  ):
    self.options = options
    self.gid = gid
    self.bcodes = bcodes
    fqdir_path = self.options.get_group_fq_dir(self.gid)
    util.mkdir_p(self.outdir)
    util.mkdir_p(fqdir_path)

  def __str__(self):
    return '{}.g{}'.format(
      self.__class__.__name__,
      self.gid,
    )

  def run(self):
    self.logger.log('extract fastq fragments for all barcodes in this group')
    fqdir_path = self.options.get_group_fq_dir(self.gid)
    readsfa_path = os.path.join(fqdir_path, 'tenxreads-bcoded.fa')
    rootfq_path = self.options.longranger_fqs_path
    tenxfq_paths = list(glob.glob(rootfq_path + '/chnk*/files/*fastq'))
    get_bcode_reads(
      tenxfq_paths,
      readsfa_path,
      self.bcodes,
    )

    self.assemble()

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def filter_alignments(in_path, out_path):
  def get_clip_info(read):
    (c, l) = read.cigar[0]
    lclip = (c in [4,5])
    (c, l) = read.cigar[-1]
    rclip = (c in [4,5])
  
    return (lclip, rclip)
  
  fin = pysam.Samfile(in_path, 'rb')
  fullrid_set = set()
  for read in fin:
    rid = (read.qname, read.is_read1)
    if read.is_unmapped:
      continue
    (lclip, rclip) = get_clip_info(read)
    if not (lclip or rclip):
      fullrid_set.add(rid)
  fin.close()
  
  fin = pysam.Samfile(in_path, 'rb')
  fout = pysam.Samfile(out_path, 'wb', template=fin)
  for read in fin:
    rid = (read.qname, read.is_read1)
    prid = (read.qname, not read.is_read1)
    if rid in fullrid_set or prid in fullrid_set:
      fout.write(read)
    #if rid in fullrid_set and prid in fullrid_set:
    #  fout.write(read)
  fin.close()
  fout.close()
  return

def get_filtered_fa(
  filtbam_path,
  infa_path,
  outfa_path,
):
  # load all filtered query names to accept
  qname_set = set()
  fin = pysam.Samfile(filtbam_path, 'rb')
  for read in fin:
    qname_set.add(read.qname)
  fin.close()

  with open(outfa_path, 'w') as fout:
    for e in util.fa_iter(infa_path):
      if e.qname_full in qname_set:
        fout.write(e.txt)

def get_lrhints_fa(
  ctgs,
  infa_path,
  outfa_path,
):
  with open(outfa_path, 'w') as fout:
    ctg_fasta = pysam.FastaFile(infa_path)
    for ctg in ctgs:
      seq = str(ctg_fasta.fetch(ctg).upper())
      for _ in xrange(10):
        fout.write('>{}\n{}\n'.format(ctg, seq))

def get_bcode_reads(infq_paths, outfa_path, bcode_set):
  seen_set = set()
  with open(outfa_path, 'w') as fout:
    for fq_path in infq_paths:
      with FastqIndex(fq_path) as idx:
        for bcode in bcode_set & idx.bcode_set:
          for e in idx.get_reads(bcode):
            # tag qname with barcode
            nqname = '{}${}'.format(bcode, e.qname)
            fout.write('>{}\n'.format(nqname))
            fout.write('{}\n'.format(e.seq1))
            fout.write('>{}\n'.format(nqname))
            fout.write('{}\n'.format(e.seq2))
            seen_set.add(e.bcode)
  assert seen_set.issubset(bcode_set), 'not all barcodes loaded'
  return
