#!/usr/bin/python

import qctoolkit as qtk
import re

print 'VOID:'
print ' name: VOID'
print ' period: 0'
print ' group: 0'
print ' atomic_number: 0'
print ' valence_electrons: 0'
print ' atomic_weight: 0\n'

pattern = re.compile('\d[spdf]')
for e in qtk.element:
  ve = 0
  ecfg_lst = e.eleconfig.split()
  ecfg_str = filter(pattern.match, ecfg_lst)
  for s in ecfg_str[-2:]:
    n = re.sub(pattern, '', s)
    try:
      ve = ve + int(n)
    except ValueError:
      ve = ve + 1

  if e.period > 2:
    if e.group == 18:
      ve = 8
    elif e.group < 3:
      ve = ve + 8
    else:
      if e.period == 4:
        if e.group < 11:
          ve = ve + 8
      elif e.period == 5:
        if e.group < 8:
          ve = ve + 8
      elif e.period == 6:
        if e.group < 5:
          ve = ve + 8
    
  elif e.period > 1:
    if e.group < 13:
      ve = ve + 2

  print '%s:' % e.symbol
  print ' name: %s' % e.name
  print ' period: %d' % e.period
  print ' group: %d' % e.group
  print ' atomic_number: %d' % e.number
  print ' valence_electrons: %d' % ve
  print ' atomic_weight: %f\n' % e.mass

#  msg = "%-2s %3d: %-22s %2d (%d, %2d)" %\
#    (e.symbol, e.number, e.eleconfig, ve, e.period, e.group)
#  print msg
