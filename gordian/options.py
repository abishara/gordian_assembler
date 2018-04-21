import abc
import os
import json

from gordian.mlib import util

class classproperty(object):
    def __init__(self, f):
        self.f = f
    def __get__(self, obj, owner):
        return self.f(owner)

class ClusterSettings(object):
    def __init__(self):

        # local
        self.cluster_type = "local"
        self.processes = 1
        self.cluster_options = {}

    @staticmethod
    def deserialize(options_dict):
        settings = ClusterSettings()

        if "processes" in options_dict:
            settings.processes = options_dict["processes"]
        if "cluster_type" in options_dict:
            settings.cluster_type = options_dict["cluster_type"]
        if "cluster_options" in options_dict:
            settings.cluster_options = options_dict["cluster_options"]

        return settings

#--------------------------------------------------------------------------
# options base
#--------------------------------------------------------------------------
class Options(object):
    __metaclass__ = abc.ABCMeta

    @classproperty
    def pipe_type(self): return None

    @classproperty
    def required(self):
      return [
        'bam_path',
        'vcf_path',
        'fq_path',
      ]

    @classproperty
    def optional(self):
      return [
        (
          'phasing_args', {
            'K' : 10,
          },
        )
      ]

    def __init__(self, options_path, **kwdargs):
        self.options_path = options_path
        if os.path.dirname(options_path) == '':
          self._output_dir = os.getcwd()
        else:
          self._output_dir = os.path.dirname(self.options_path)

        # set required attributes to None
        for opt in self.required:
          setattr(self, opt, None)

        # set optional to default
        for opt, val in self.optional:
          setattr(self, opt, val)

        self.cluster_settings = ClusterSettings()

    @classmethod
    def deserialize(cls, options_path):
      # load json config
      with open(options_path) as f:
        options_dict = json.load(f)

      options = cls(options_path)
      # required
      for opt in cls.required:
        assert opt in options_dict, 'required option "{}" missing'.format(opt)
        setattr(options, opt, options_dict[opt])

      # optional
      for opt, val in cls.optional:
        setattr(options, opt, options_dict.get(opt, val))

      # cluster settings
      options.cluster_settings = ClusterSettings.deserialize(
        options_dict.get("cluster_settings", {}))

      return options

    @property
    def output_dir(self):
        return self._output_dir
    
    @property
    def results_dir(self):
        return os.path.join(self.output_dir, "results")

    @property
    def working_dir(self):
        return os.path.join(self.output_dir, "working")

    @property
    def log_dir(self):
        return os.path.join(self.output_dir, "logs")

    @property
    def phased_bins_path(self):
        return os.path.join(self.results_dir, "clusters.p")

    def get_bin_scratch(self, idx):
      return os.path.join(self.working_dir, 'c{}'.format(idx))

    def __getstate__(self):
        """
        allows pickling of Options instances, necessary for ipyparallel
        """
        state = self.__dict__.copy()
        return state

