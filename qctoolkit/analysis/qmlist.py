import qctoolkit as qtk

class QMList(object):
  def __init__(self, file_list, **kwargs):
    self.data = []
    self.name_list = []
    self.data_dict = {}
    for out in file_list:
      qmout = qtk.QMOut(out, program = kwargs['program'])
      self.data_dict[qmout.name] = qmout
      self.name_list.append(qmout.name)
      self.data.append(qmout)

  def __repr__(self):
    return str(self.data)

  def __getitem__(self, key):
    if type(key) is int:
      return self.data[key]
    elif type(key) is str:
      match_list = filter(lambda x: key in x, self.name_list)
      if len(match_list) > 0:
        return [self.data_dict[name] for name in match_list]
      
    
