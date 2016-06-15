import qctoolkit as qtk
import copy

class QMList(object):
  def __init__(self, file_list = None, **kwargs):
    self.data = []
    self.name_list = []
    self.data_dict = {}
    if file_list:
      if isinstance(file_list, qtk.QMList):
        self.data = file_list.data
        self.name_list = file_list.name_list
        self.data_dict = file_list.data_dict
      else:
        for out in file_list:
          qmout = qtk.QMOut(out, program = kwargs['program'])
          self.data_dict[qmout.name] = qmout
          self.name_list.append(qmout.name)
          self.data.append(qmout)

  def __repr__(self):
    return str(self.data)

  def __getitem__(self, key):
    out = qtk.QMList()
    if type(key) is str:
      match_list = filter(lambda x: key in x, self.name_list)
      if len(match_list) > 0:
        out.name_list = match_list
        for name in out.name_list:
          out.data.append(self.data_dict[name])
          out.data_dict[name] = self.data_dict[name]
    else:
      out = qtk.QMList()
      out.data = self.data[key]
      out.name_list = self.name_list[key]
      for i in range(len(out.data)):
        out.data_dict[out.name_list[i]] = out.data[i]
    return out
