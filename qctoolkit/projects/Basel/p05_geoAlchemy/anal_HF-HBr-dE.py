#!/usr/bin/python

import qctoolkit as qtk
import matplotlib.pyplot as pl

Ed_ver = qtk.QMData('03*', '*-d??-*', 'cpmd').kcal()
Er_ver = qtk.QMData('03*', '*-ref*', 'cpmd').kcal()
Ed_non = qtk.QMData('04*', '*-d??-*', 'cpmd').kcal()
Er_non = qtk.QMData('04*', '*-ref*', 'cpmd').kcal()

Ed_ver.extract('.*-l(...)\.out')
Er_ver.extract('.*-l(...)\.out')
Ed_non.extract('.*-l(...)\.out')
Er_non.extract('.*-l(...)\.out')

Ed_ver.index_div(100)
Er_ver.index_div(100)
Ed_non.index_div(100)
Er_non.index_div(100)

dE_ver = Ed_ver - Er_ver
dE_non = Ed_non - Er_non

fig = pl.figure('DE')
ver = dE_ver.E.plot(marker='o', color='black')
non = dE_non.E.plot(marker='s', color='red')
pl.legend(['vertical', 'non-vertical'], numpoints=1)
pl.xlabel('$\lambda$',fontsize=18)
pl.ylabel('$\Delta E$ [kcal/mol]', fontsize=18)
pl.show()
