import pkgutil

eggs_loader = pkgutil.find_loader('horton')
found = eggs_loader is not None
if found:
  try:
    from horton import GOBasisFamily
    import qctoolkit as qtk
    import os
    package_path = os.path.split(qtk.__file__)[0]
    module_path = 'data/basis_set'
    path = os.path.join(package_path, module_path)

    Ne_def2_tzvp_uncontracted = GOBasisFamily(
      'H_He_basis', 
      filename='%s/Def2-TZVP/Ne_uncontracted.nwchem' % path
    )

    Ne_def2_svp_uncontracted = GOBasisFamily(
      'H_He_basis', 
      filename='%s/Def2-SVP/Ne_uncontracted.nwchem' % path
    )

    Ne_sto2g_uncontracted = GOBasisFamily(
      'H_He_basis', 
      filename='%s/sto-2g/Ne_uncontracted.nwchem' % path
    )
  except:
    pass
else:
  pass



