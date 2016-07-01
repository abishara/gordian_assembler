--------------------------------------------------------------------------
install
--------------------------------------------------------------------------

install with develop to modify have the changes in the source reflected in
the run

% python setup.py develop

The following must be visible in $PATH:
* IDBA (https://github.com/loneknightpy/idba)
* samtools (assumes format of version:0.1.19)
* BWA

--------------------------------------------------------------------------
run
--------------------------------------------------------------------------

% python main.py <path/to/scratch-dir/config.json>

The config.json should be placed in a scratch directory where all logs,
intermediate files, and results will be generated.

At the end of the ref-asm pipeline the assembly results will be in:

<path/to/scratch-dir/results/bins/*/>
  contig.fa
  local-asm-merged.fa  
  scaff.fasta

* contig.fa: idba assembled contigs
* local-asm-merged.fa: contigs produced from local reassembly using the
  barcodes.  (this is a work in progress, most of the contigs are
  identitcal to those in contig.fa, but some may be extended and assembled
  furhter).
* scaff.fasta: scaffolded contigs from contig.fa.  (this is also a work in
  progress, and others have found it sometimes inverts the contigs).

The pipeline can be run with three different modes:

* local: run on the local host
* multiprocessing: run tasks using the python multiprocessing module (on a multicore
  interactive node)
* IPCluster: run tasks using ipython cluster helper ontop of SGE, SLURM
  (tested), torque, and lsf (untested).

