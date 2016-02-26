import cpmd 

pp_setting = {
               'program': 'cpmd',
               'xc': 'pbe',
               'method': 'goedecker',
             }
class PP(object):
  def __init__(self, path=None, **kwargs):
    self.param = {}
    self.setting = {}
    self.info = ''

    for key, value in pp_setting.iteritems():
      if key in kwargs:
        self.setting[key] = kwargs[key]
      else:
        self.setting[key] = value
    if path:
      self.read(path)

  def read(self, path):
    if self.setting['program'] == 'cpmd':
      cpmd.read(self, path)

    else:
      qtk.exit('program %s is not implemented for PP'\
        % self.setting['program'])

  def write(self, name=None):
    if self.setting['program'] == 'cpmd':
      cpmd.write(self, name)

    else:
      qtk.exit('program %s is not implemented for PP'\
        % self.setting['program'])
