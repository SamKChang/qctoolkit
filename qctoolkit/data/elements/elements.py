import yaml, os, re
import qctoolkit as qtk

class Elements(object):
  path = re.sub('pyc', 'yml', os.path.realpath(__file__))
  data = yaml.safe_load(open(path))
  def __init__(self):
    pass

  @classmethod
  def ve_list(cls):
    return {e:q['valence_electrons']\
            for e, q in cls.data.iteritems()}

  @classmethod
  def z_list(cls):
    return {e:q['atomic_number']\
            for e, q in cls.data.iteritems()}

  @classmethod
  def type_list(cls):
    return {q['atomic_number']:e\
            for e, q in cls.data.iteritems()}

  @classmethod
  def mass_list(cls):
    return {e:q['atomic_weight']\
            for e, q in cls.data.iteritems()}
