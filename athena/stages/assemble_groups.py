import os
import pysam
import subprocess
from collections import defaultdict, deque
import glob
import shutil
import random
from bx.intervals.intersection import IntervalTree

from ..assembler_tools.architect import architect
from ..assembler_tools.architect import local_assembly

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex
from ..mlib.bam_idx import BCBamIndex

# NOTE must be in path
idbabin_path = 'idba_ud'
supernovabin_path = '/scratch/users/abishara/sources/supernova-1.0.0/supernova-cs/1.0.0/bin/supernova'

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
    paths['align.contig.bam'] = os.path.join(self.outdir, 'align.contig.bam')
    paths['align.contig.bam.bai'] = os.path.join(self.outdir, 'align.contig.bam.bai')
    #paths['local-asm-merged.fa'] = os.path.join(self.outdir, 'local-asm-merged.fa')
    #paths['architect-scaff.fasta'] = os.path.join(self.outdir, 'scaff.fasta')
    return paths

  def clean_working(self):
    bin_path = self.options.get_group_dir(self.gid)
    self.logger.log('removing bin directory {}'.format(bin_path))
    shutil.rmtree(bin_path)
    return 

  def assemble(self):
    self.logger.log('assembling barcoded reads for this bin')

    fqdir_path = self.options.get_group_fq_dir(self.gid)
    readsfa_path = os.path.join(fqdir_path, 'tenxreads-bcoded.fa')
    asmdir_path = self.options.get_group_asm_dir(self.gid)
    # initial idba assembly
    cmd = '{} -r {} -o {} --maxk 120'.format(
      idbabin_path,
      readsfa_path,
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
      cmd = 'cat align.on-contig.sam | samtools view -bS - | samtools sort -o align.on-contig.sorted.bam -'
      subprocess.check_call(cmd, shell=True)

      # filter bad reads, create idba *fa reads with only filtered reads
      allbam_path = 'align.on-contig.sorted.bam'
      filtbam_path = 'align.on-contig.sorted.filt.bam'
      #shutil.copyfile(allbam_path, filtbam_path)
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

      architect.do_scaffolding(
        'contig.fa',
        'edges.tsv',
        'containment',
        0,
        'scaff',
      )

    # add group uid to contig and scaffold names
    contig_path = os.path.join(asmdir_path, 'contig.fa')
    ncontig_path = os.path.join(asmdir_path, 'contig.rename.fa')
    with open(ncontig_path, 'w') as fout:
      fasta = pysam.FastaFile(contig_path)
      for ctg in sorted(fasta.references):
        nctg = '{}.g{}'.format(ctg.strip().split()[0], self.gid)
        seq = str(fasta.fetch(ctg).upper())
        fout.write('>{}\n'.format(nctg))
        fout.write(str(seq) + '\n')

    # align renamed contigs to reference
    self.logger.log('aligning renamed contigs to reference')
    with util.cd(asmdir_path):
      cmd = 'bwa mem -p {} contig.rename.fa > align.contig.sam'.format(
        self.options.ref_fasta,
      )
      subprocess.check_call(cmd, shell=True)
      cmd = 'cat align.contig.sam | samtools view -bS - | samtools sort -o align.contig.bam -'
      subprocess.check_call(cmd, shell=True)
      cmd = 'samtools index align.contig.bam'
      subprocess.check_call(cmd, shell=True)

    def copyfinal(src, dest):
      shutil.copyfile(
        os.path.join(self.options.get_group_asm_dir(self.gid), src),
        os.path.join(self.options.get_group_dir(self.gid, final=True), dest),
      )

    self.logger.log('copying deliverables to final')
    copyfinal('contig.rename.fa', 'contig.fa')
    copyfinal('align.contig.bam', 'align.contig.bam')
    copyfinal('align.contig.bam.bai', 'align.contig.bam.bai')
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
      #if gid != 999:
      #  continue

      yield cls(options, gid, bins, bcodes)

  def outpaths(self, final=False):
    paths = {}
    paths['supernova.contig.fa'] = os.path.join(self.outdir,'supernova.contig.fa')
    paths['supernova.scaff.fa']  = os.path.join(self.outdir,'supernova.scaff.fa')
    paths['supernova.bam']       = os.path.join(self.outdir,'supernova.bam')
    #paths['shit'] = 'shit'
    return paths

  def __init__(
    self,
    options,
    gid,
    bins,
    bcodes,
  ):
    self.options = options
    self.gid = gid
    self.bins = bins
    self.bcodes = bcodes
    fqdir_path = self.options.get_group_fq_dir(self.gid)

    print 'bins', len(bins)
    print 'bcodes', len(bcodes)
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
    supernova_dir_path = os.path.join(fqdir_path, 'supernova-input-fqs')
    util.mkdir_p(supernova_dir_path)
    supernova_ra_fq_path = os.path.join(supernova_dir_path,
      'read-RA_si-CCTGGAGA_lane-001-chunk-000.fastq')
    supernova_si_fq_path = os.path.join(supernova_dir_path,
      'read-I1_si-CCTGGAGA_lane-001-chunk-000.fastq')
    idba_fa_path = os.path.join(fqdir_path, 'idba-bcoded-reads.fa')

    bcodes = self.bcodes
    self.logger.log('obtaining read alignments for barcodes of interest')

    readsbampre_path = os.path.join(fqdir_path, 'reads-pre.bam')
    readsbam_path = os.path.join(fqdir_path, 'reads.bam')
    bcbam_path = self.options.longranger_bcbam_path
    fhandle = pysam.Samfile(bcbam_path, 'rb')
    fhandle_out = pysam.Samfile(readsbampre_path, 'wb', template=fhandle)
    with BCBamIndex(bcbam_path) as idx:
      for bcode in bcodes:
        if bcode not in idx.bcode_set:
          self.logger.log('  - alignments for barcode {} not found'.format(
            bcode,
          ))
        for read in idx.get_reads(bcode):
          fhandle_out.write(read)

    fhandle.close()
    fhandle_out.close()

    self.logger.log('  - sorting')
    cmd = 'samtools sort -o {} {}'.format(
      readsbam_path,
      readsbampre_path,
    )
    subprocess.check_call(cmd, shell=True)
    cmd = 'samtools index {}'.format(readsbam_path)
    subprocess.check_call(cmd, shell=True)

    self.logger.log('  - filtering reads based on alignment positions')
    filtbam_path = os.path.join(fqdir_path, 'filt-reads.bam')
    get_filtered_bam(
      self.bins,
      readsbam_path,
      filtbam_path,
    )
    cmd = 'samtools index {}'.format(filtbam_path)
    subprocess.check_call(cmd, shell=True)

    passqname_set = set()
    fhandle = pysam.Samfile(filtbam_path, 'rb')
    for read in fhandle.fetch(until_eof=True):
      passqname_set.add(read.qname)
    fhandle.close()

    self.logger.log('  - get reads from fastq')
    get_bcode_reads(
      readsbam_path,
      passqname_set,
      supernova_ra_fq_path,
      supernova_si_fq_path,
      idba_fa_path,
    )

    supernova_contig_path, supernova_scaff_path = self.supernova_assemble()
    #idba_contig_path = self.idba_assemble()

    # align contigs to reference
    asmdir_path = self.options.get_group_asm_dir(self.gid)
    supernova_bam_path = os.path.join(asmdir_path, 'supernova.bam')
    cmd = 'bwa mem -p {} {} | samtools view -bS - | samtools sort -o {} -'.format(
      self.options.ref_fasta,
      supernova_contig_path,
      supernova_bam_path,
    )
    subprocess.check_call(cmd, shell=True)
    cmd = 'samtools index {}'.format(supernova_bam_path)
    subprocess.check_call(cmd, shell=True)

    #idba_bam_path = os.path.join(asmdir_path, 'idba.bam')
    #cmd = 'bwa mem -p {} {} | samtools view -bS - | samtools sort -o {} -'.format(
    #  self.options.ref_fasta,
    #  idba_contig_path,
    #  idba_bam_path,
    #)
    #subprocess.check_call(cmd, shell=True)
    #cmd = 'samtools index {}'.format(idba_bam_path)
    #subprocess.check_call(cmd, shell=True)

    def copyfinal(src, dest):
      shutil.copyfile(
        os.path.join(self.options.get_group_asm_dir(self.gid), src),
        os.path.join(self.options.get_group_dir(self.gid, final=True), dest),
      )

    self.logger.log('copying deliverables to final')
    copyfinal(supernova_contig_path, 'supernova.contig.fa')
    copyfinal(supernova_scaff_path, 'supernova.scaff.fa')
    copyfinal(supernova_bam_path, 'supernova.bam')
    copyfinal(supernova_bam_path+'.bai', 'supernova.bam.bai')
    #copyfinal(idba_contig_path, 'idba.fa')
    #copyfinal(idba_bam_path, 'idba.bam')
    #copyfinal(idba_bam_path+'.bai', 'idba.bam.bai')
    return

  def supernova_assemble(self):
    fqdir_path = self.options.get_group_fq_dir(self.gid)
    supernova_fqdir_path = os.path.join(fqdir_path, 'supernova-input-fqs')
    supernova_ra_fq_path = os.path.join(supernova_fqdir_path,
      'read-RA_si-CCTGGAGA_lane-001-chunk-000.fastq')
    supernova_si_fq_path = os.path.join(supernova_fqdir_path,
      'read-I1_si-CCTGGAGA_lane-001-chunk-000.fastq')

    # gzip reads
    supernova_ra_fqgz_path = supernova_ra_fq_path + '.gz'
    supernova_si_fqgz_path = supernova_si_fq_path + '.gz'
    cmd = 'gzip {}'.format(supernova_ra_fq_path)
    subprocess.check_call(cmd, shell=True)
    cmd = 'gzip {}'.format(supernova_si_fq_path)
    subprocess.check_call(cmd, shell=True)

    asmdir_path = self.options.get_group_asm_dir(self.gid)
    supernova_asmdir_path = os.path.join(
      self.options.get_group_asm_dir(self.gid),
      'supernova-asm',
    )
    util.mkdir_p(supernova_asmdir_path)
    with util.cd(supernova_asmdir_path):
      cmd = '{} run --id=run --fastqs={}'.format(
        supernovabin_path,
        supernova_fqdir_path,
      )
      subprocess.check_call(cmd, shell=True)
      cmd = '''{} mkfasta \
--asmdir={} \
--outprefix={} \
--style={}'''.format(
        supernovabin_path,
        'run/outs/assembly',
        'supernova',
        #'megabubbles',
        'pseudohap',
      )
      subprocess.check_call(cmd, shell=True)
    scaff_path = os.path.join(
      supernova_asmdir_path,
      'supernova.fasta',
    )
    assert os.path.isfile(scaff_path), "supernova contigs not generated"
    contig_path = os.path.join(
      supernova_asmdir_path,
      'supernova.contig.fasta',
    )
    split_scaffs(scaff_path, contig_path)
    return contig_path, scaff_path

  def idba_assemble(self):
    fqdir_path = self.options.get_group_fq_dir(self.gid)
    idba_fa_path = os.path.join(fqdir_path, 'idba-bcoded-reads.fa')
    idba_asmdir_path = os.path.join(
      self.options.get_group_asm_dir(self.gid),
      'idba-asm',
    )
    cmd = '{} -r {} -o {} --maxk 120'.format(
      idbabin_path,
      idba_fa_path,
      idba_asmdir_path,
    )
    subprocess.check_call(cmd, shell=True)
    contig_path = os.path.join(idba_asmdir_path, 'contig.fa')
    return contig_path

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

def get_bcode_reads(
  inbam_path,
  passqname_set,
  supernova_ra_fq_path,
  supernova_si_fq_path,
  idba_fa_path,
):
  seen_set = set()
  fhandle = pysam.Samfile(inbam_path, 'rb')

  def get_random_bcode():
    chars = list('ACGT')
    bcode = ''.join([random.choice(chars) for _ in xrange(16)]) + '-1'
    return bcode

  def get_bcode(read):
    filt_list = filter(lambda(k, v): k == 'BX', read.tags)
    if filt_list == []: 
      return None
    else:
      k, v = filt_list[0]
      return v

  def reverse_complement(bases):
    lut = {
      'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n',
      'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N',
    }
    nbases = reversed([lut[b] for b in bases])
    nbases = ''.join(nbases)
    return nbases

  r1_map = {}
  r2_map = {}
  for read in fhandle.fetch(until_eof=True):
    if read.qname not in passqname_set:
      continue
    if read.is_secondary:
      continue
    qname = read.qname
    seq = read.query_sequence
    qual = read.query_qualities
    qual = read.qual
    bcodepre = get_bcode(read)
    if bcodepre == None:
      #continue
      bcodepre = get_random_bcode()
    bcode = bcodepre.strip('-1')
    if read.is_reverse:
      qual = qual[::-1]
      seq = reverse_complement(seq)
    if read.is_read1:
      r1_map[qname] = (seq, qual, bcode)
    else:
      r2_map[qname] = (seq, qual)
    seen_set.add(qname)
  fhandle.close()

  # FIXME uncomment
  #assert passqname_set.issubset(seen_set), "not all reads extracted from *bam"
  # FIXME make ==
  assert len(r1_map) <= len(r2_map), "not all read pairs properly extracted from *bam"

  with open(supernova_ra_fq_path, 'w') as fout_ra, \
       open(supernova_si_fq_path, 'w') as fout_si, \
       open(idba_fa_path, 'w') as fout_idba:
    for i, qname in enumerate(r1_map):
      #if i == 1000:
      #  break
      r1seq, r1qual, bcode = r1_map[qname]
      r2seq, r2qual = r2_map[qname]
      fout_ra.write('@{} 1:N:0:0\n{}\n+\n{}\n'.format(
        qname,
        bcode + 'NNNNNNN' + r1seq,
        'K' * (16 + 7) + r1qual,
      ))
      fout_ra.write('@{} 3:N:0:0\n{}\n+\n{}\n'.format(
        qname,
        r2seq,
        r2qual,
      ))
      fout_si.write('@{} 2:N:0:0\n{}\n+\n{}\n'.format(
        qname,
        'CCTGGAGA',
        'AAFFFKKK',
      ))

      # tag qname with barcode
      nqname = '{}${}'.format(bcode, qname)
      fout_idba.write('>{}\n'.format(nqname))
      fout_idba.write('{}\n'.format(r1seq))
      fout_idba.write('>{}\n'.format(nqname))
      fout_idba.write('{}\n'.format(r2seq))

  return

def get_filtered_bam(
  bins,
  inbam_path,
  outbam_path,
):

  SLACK_LEN = 10000000
  SLACK2_LEN = 1000
  intervals_map = defaultdict(lambda: IntervalTree(1,1))
  for i, (ctg, b, e, _) in enumerate(bins):
    intervals_map[ctg].insert(b, e, i)

  fhandle = pysam.Samfile(inbam_path, 'rb')
  fhandle_out = pysam.Samfile(outbam_path, 'wb', template=fhandle)

  def update_neighbors(neighbors, read):
    assert not read.is_unmapped
    currctg = read.tid if neighbors == [] else neighbors[0].tid
    if currctg == read.tid:
      neighbors.append(read)
      clipidx = 0
      for _read in neighbors:
        if read.pos - _read.pos < SLACK2_LEN:
          break
        clipidx += 1
      neighbors = neighbors[clipidx:]
      assert len(neighbors) > 0
    else:
      neighbors = []
    return neighbors

  neighbors = []
  for i, read in enumerate(fhandle.fetch(until_eof=True)):
    if i % 1000000 == 0:
      print '  ', i
    # save all unmapped reads
    if read.is_unmapped:
      fhandle_out.write(read)
      continue
    ctg = fhandle.getrname(read.tid)
    # save any reads within 10mb of target
    hits = intervals_map[ctg].find(read.pos, read.aend)
    if len(hits) > 0:
      fhandle_out.write(read)
      continue

    ## if any remaining read is mq60 and its neighborhood coverage is < 4x
    ## then filter
    #neighbors = update_neighbors(neighbors, read)
    ##print map(lambda(r): r.pos, neighbors)
    #if (
    #  read.mapq >= 60 and 
    #  len(neighbors) * 140. / SLACK2_LEN < 3
    #):
    #  continue

    #fhandle_out.write(read)
  
  fhandle.close()
  fhandle_out.close()

def split_scaffs(scaff_path, contig_path):
  fasta = pysam.FastaFile(scaff_path)
  with open(contig_path, 'w') as fout:
    for i, scaff in enumerate(fasta.references):
      seq = str(fasta.fetch(scaff).upper())
      cidx = 0
      for subseq in seq.split('N'):
        if len(subseq) > 0:
          fout.write('>{}.{}\n'.format(scaff, cidx))
          fout.write('{}\n'.format(subseq))
          cidx +=1
  return

