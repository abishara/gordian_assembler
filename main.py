import os
import sys
import collections
import logging

from athena.mlib import util

from athena import pipeline
from athena.options import Options

from athena.stages import index_reads
from athena.stages import phase_barcodes
from athena.stages import assemble_bins

logging.basicConfig(format='%(message)s', level=logging.DEBUG)

def run(options):
    """
    1. create output directories
    2. collect args for each stage
    3. check which stages need to run
    4. iterate through stages and submit jobs
    5. validate that we're done running
    """

    util.mkdir_p(options.output_dir)
    util.mkdir_p(options.results_dir)
    util.mkdir_p(options.working_dir)
    util.mkdir_p(options.log_dir)

    stages = get_stages(options)
    runner = pipeline.Runner(options)

    for stage_name, stage in stages.items():
        runner.run_stage(stage, stage_name)

def clean(options):
  stages = get_stages(options)

  for stage_name, stage in stages.items():
    stage.clean_all_steps(options)

def get_stages(options):
    stages = collections.OrderedDict()

    stages["index_reads"]  = index_reads.IndexReadsStep
    stages["phase_barcodes"] = phase_barcodes.PhaseBarcodesStep
    stages["assemble_bins"] = assemble_bins.AssembleBinsStep

    return stages

def clean_up():
  junk = filter(
    lambda(f): (
      f.startswith('SLURM_controller') or 
      f.startswith('SLURM_engine') or 
      f.startswith('sge_controller') or 
      f.startswith('sge_engine') or 
      f.startswith('bcbio-')
    ),
    os.listdir('.'),
  )
  map(lambda(f): os.remove(f), junk)

def main():
  """
  1. process command-line arguments
  2. run
  """

  help_str = '''
  usage: athena.py <path/to/config.json>

  NOTE: dirname(config.json) specifies root output directory
  '''
  argv = sys.argv
  if len(argv) != 2 and len(argv) != 3:
    print help_str
    sys.exit(1)

  config_path = argv[1]
  options = Options.deserialize(config_path)

  #clean(options)
  run(options)

  clean_up()

if __name__ == '__main__':
  main()

