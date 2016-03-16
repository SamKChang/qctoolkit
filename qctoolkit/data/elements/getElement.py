#!/usr/bin/python

import urllib2
import re

url = 'https://en.wikipedia.org/wiki/List_of_elements'
page = urllib2.urlopen(url).readlines()

ind = []
for i in range(1, 94):
  s = '<td>%d</td>' % i
  s_str = filter(lambda x: s in x, page)[0]
  ind.append(page.index(s_str))

for i in range(len(ind) - 1):
  j = ind[i]
  k = ind[i+1]
  print 'element %d:' % (i + 1)
  for l in range(j, k):
    print page[l],
