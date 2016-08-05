import qctoolkit as qtk
import pkgutil
import grid_points as gp
from libxc_dict import xc_dict
import os


selfPath = os.path.realpath(__file__)
selfPath = os.path.split(selfPath)[0]
xcpath = os.path.join(selfPath, 'libxc_exc.so')
xc_found = os.path.exists(xcpath)
if xc_found:
  from libxc_exc import libxc_exc
  from libxc_vxc import libxc_vxc
else:
  pass

def libxc_report(self, xc_id, flag):
  for k, v in xc_dict.iteritems():
    if v == xc_id:
      key = k
      qtk.report("libxc_%s" % flag, "xc: %s, id: %d\n" % (key, xc_id))
      break

def exc(self, xcFlag=1, rhoSigma = None, report=True):
  if type(xcFlag) is int:
    if xcFlag not in xc_dict.values():
      qtk.exit("libxc functional id number %d is not valid" % xcFlag)
    else:
      xc_id = xcFlag
  elif type(xcFlag) is str:
    if xcFlag not in xc_dict:
      qtk.exit("libxc functional id %s is not valid" % xcFlag)
    else:
      xc_id = xc_dict[xcFlag]
  if report:
    self.libxc_report(xc_id, 'exc')
  coords = self.grid.points
  if rhoSigma is None:
    rho = gp.getRho(self, coords)
    sigma = gp.getSigma(self, coords)
  else:
    rho, sigma = rhoSigma
  return libxc_exc(rho, sigma, len(coords), xc_id)

def vxc(self, xcFlag=1, rhoSigma = None, report=True):
  if type(xcFlag) is int:
    if xcFlag not in xc_dict.values():
      qtk.exit("libxc functional id number %d is not valid" % xcFlag)
    else:
      xc_id = xcFlag
  elif type(xcFlag) is str:
    if xcFlag not in xc_dict:
      qtk.exit("libxc functional id %s is not valid" % xcFlag)
    else:
      xc_id = xc_dict[xcFlag]
  if report:
    self.libxc_report(xc_id, 'vxc')
  coords = self.grid.points
  if rhoSigma is None:
    rho = gp.getRho(self, coords)
    sigma = gp.getSigma(self, coords)
  else:
    rho, sigma = rhoSigma
  return libxc_vxc(rho, sigma, len(coords), xc_id)
