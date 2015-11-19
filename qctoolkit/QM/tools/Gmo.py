#!/usr/bin/python

import numpy as np
import scipy.special as sp
import sys
import pdb
import os
import os.path
from copy import deepcopy

##########################################
#  Solid harmonic transformation matrix  #
##########################################
def Pl(l,I):
  ##########################
  #  NOT COMPLETED!!!!!!   #
  ##########################
  if l==2:
    blk = np.array([
     [-0.315392,-0.315392, 0.630783, 0.0     , 0.0     , 0.0     ],\
     [ 0.0     , 0.0     , 0.0     , 0.0     , 1.092548, 0.0     ],\
     [ 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 1.092548],\
     [ 0.546274,-0.546274, 0.0     , 0.0     , 0.0     , 0.0     ],\
     [ 0.0     , 0.0     , 0.0     , 1.092548, 0.0     , 0.0     ]
    ])
  for i in range(0,len(AOMAP)):
    if (MAP[AOMAP[i]-1]==I) and (sum(lmMAP[AOMAP[i]-1])==l):
      print "MATCH!"
      
  return np.dot(np.transpose(blk),blk)


##########################
#  matlab style blkdiag  #
##########################
def blkdiag(*arrs):
  arrs = [np.asarray(a) for a in arrs]
  shapes = np.array([a.shape for a in arrs])
  out = np.zeros(np.sum(shapes, axis=0))

  r, c = 0, 0
  for i, (rr, cc) in enumerate(shapes):
      out[r:r + rr, c:c + cc] = arrs[i]
      r += rr
      c += cc
  return out

#########################
#  fractorial function  #
#########################
def factorial(n):
  return np.prod(range(1,n+1))

#######################
#  factorial2 source  #
#######################
def factorial2(n, exact=False):
    if exact:
        if n < -1:
            return 0
        if n <= 0:
            return 1
        val = 1
        for k in xrange(n,0,-2):
            val *= k
        return val
    else:
        n = asarray(n)
        vals = zeros(n.shape,'d')
        cond1 = (n % 2) & (n >= -1)
        cond2 = (1-(n % 2)) & (n >= -1)
        oddn = extract(cond1,n)
        evenn = extract(cond2,n)
        nd2o = oddn / 2.0
        nd2e = evenn / 2.0
        place(vals,cond1,gamma(nd2o+1)/sqrt(pi)*pow(2.0,nd2o+0.5))
        place(vals,cond2,gamma(nd2e+1) * pow(2.0,nd2e))
        return vals

#######################################
#  Gaussian overlap for ECP operator  #
#######################################
def PP_ij(ai,aj,alpha,PP_I):

  def dPEXP(t,q,P,C):
    px = -q*(P-C)**2
    if t==0:
      return np.exp(px)
    elif t==1:
      factor = q*(P-C)
      return 2*np.exp(px)*factor
    elif t==2:
      factor = 2*(q*(P-C))**2 - 1
      return 2*np.exp(px)*factor
    else:
      sys.exit("pseudo potentials for d orbitals is not implemented")

  Sab = 0
  AOi  = AOMAP[ai] - 1
  AOj  = AOMAP[aj] - 1
  I = MAP[AOi] - 1
  J = MAP[AOj] - 1
  A = COORD[I]
  B = COORD[J]
  C = COORD[PP_I]
  AB= A-B
  m = lmMAP[ai]
  n = lmMAP[aj]
  L0 = m[0]+n[0]
  L1 = m[1]+n[1]
  L2 = m[2]+n[2]
  # loop through gaussian components for each AO
  for gi in range(AGMAP[AOi], AGMAP[AOi+1]):
    for gj in range(AGMAP[AOj], AGMAP[AOj+1]):
      a  = PEXP[gi]
      b  = PEXP[gj]
      ci = CTRC[gi]*N(a,m[0],m[1],m[2])
      cj = CTRC[gj]*N(b,n[0],n[1],n[2])
      cij= ci*cj
      p  = a+b
      q  = p*alpha/(p+alpha)
      p2 = 2*p
      mu = a*b/float(p)
      P  = (a*A+b*B)/float(p)
      PA = P - A
      PB = P - B
      for t in range(0, L0+1):
        Ex = E(m[0],n[0],t,p2,mu,AB[0],PA[0],PB[0])
        dPx = dPEXP(t,q,P[0],C[0])
        for u in range(0, L1+1):
          Ey = E(m[1],n[1],u,p2,mu,AB[1],PA[1],PB[1])
          dPy = dPEXP(u,q,P[1],C[1])
          for v in range(0, L2+1):
            # print E, R function format:
            #  def E(m,n,t,p2,mu,ABi,PAi,PBi)
            #  def R(t,u,v,n,PC,p2,x)
            Ez = E(m[2],n[2],v,p2,mu,AB[2],PA[2],PB[2])
            dPz = dPEXP(v,q,P[2],C[2])
            Sab += cij*(Ex*Ey*Ez)*(dPx*dPy*dPz)*(np.pi/float(p+alpha))
  return Sab

####################################
#  Hermite expansion coefficients  #
####################################
def E(m,n,t,p2,mu,Xab,Xpa,Xpb):
  if((t>m+n)|(t<0)|(m<0)|(n<0)):
    Emn = 0
  elif((m>=n)&(m>0)):
    E0 = E(m-1,n,t-1,p2,mu,Xab,Xpa,Xpb)/float(p2)
    E1 = E(m-1,n,t,p2,mu,Xab,Xpa,Xpb)*Xpa
    E2 = E(m-1,n,t+1,p2,mu,Xab,Xpa,Xpb)*(t+1)
    Emn = E0 + E1 + E2
  elif((n>m)&(n>0)):
    E0 = E(m,n-1,t-1,p2,mu,Xab,Xpa,Xpb)/float(p2)
    E1 = E(m,n-1,t,p2,mu,Xab,Xpa,Xpb)*Xpb
    E2 = E(m,n-1,t+1,p2,mu,Xab,Xpa,Xpb)*(t+1)
    Emn = E0 + E1 + E2
  else:
    Emn = np.float(np.exp(-mu*(Xab**2)))
  return Emn

############################
#  Normalization constant  #
############################
def N(a,li,lj,lk):
  Ni = factorial2(2*li-1,exact=True)
  Nj = factorial2(2*lj-1,exact=True)
  Nk = factorial2(2*lk-1,exact=True)
  Nijk = Ni*Nj*Nk
  nijk = li+lj+lk
  Num1 = (2*a/np.pi)**(3/float(2))
  Num2 = (4*a)**nijk
  return (Num1*Num2/float(Nijk))**(0.5)

################################
#  Coulomb-Gaussian integrals  #
################################
def coulomb_Gij(ai,aj,C,Zc):

  # Boys function  
  def F(n,x):
    if(x>0):
      G1 = sp.gamma(0.5 + n)
      G2 = sp.gammainc(0.5 + n, x)
      fac= 0.5 * x**(-n-0.5)
      return fac*G1*G2
    else:
      return 1/float(1+2*n)

  # Hermite-Coulomb integral coefficient
  def R(t,u,v,n,Rpc,p2,x):
    #pdb.set_trace()
    Rtuv = 0
    if((t==0)&(u==0)&(v==0)):
      Rtuv += F(n,x)*((-p2)**n)
    elif(((t>=u)|(t>=v))&(t>0)):
      Rtuv += R(t-1,u,v,n+1,Rpc,p2,x)*Rpc[0]
      Rtuv += R(t-2,u,v,n+1,Rpc,p2,x)*(t-1)
    elif((u>=v)&(u>0)):
      Rtuv += R(t,u-1,v,n+1,Rpc,p2,x)*Rpc[1]
      Rtuv += R(t,u-2,v,n+1,Rpc,p2,x)*(u-1)
    elif(v>0):
      Rtuv += R(t,u,v-1,n+1,Rpc,p2,x)*Rpc[2]
      Rtuv += R(t,u,v-2,n+1,Rpc,p2,x)*(v-1)
    return Rtuv


  # main function
  # loading parameters, coordinates, etc
  # additional variables passed:
  #  AOMAP, MAP, COORD, lmMAP, PEXP, CTRC
  Vp = 0
  AOi  = AOMAP[ai] - 1
  AOj  = AOMAP[aj] - 1
  I = MAP[AOi] - 1
  J = MAP[AOj] - 1
  A = COORD[I]
  B = COORD[J]
  AB= A-B
  m = lmMAP[ai]
  n = lmMAP[aj]
  L0 = m[0]+n[0]
  L1 = m[1]+n[1]
  L2 = m[2]+n[2]
  # loop through gaussian components for each AO
  for gi in range(AGMAP[AOi], AGMAP[AOi+1]):
    for gj in range(AGMAP[AOj], AGMAP[AOj+1]):
      a  = PEXP[gi]
      b  = PEXP[gj]
      ci = CTRC[gi]*N(a,m[0],m[1],m[2])
      cj = CTRC[gj]*N(b,n[0],n[1],n[2])
      cij= ci*cj
      p  = a+b
      p2 = 2*p
      mu = a*b/float(p)
      P  = (a*A+b*B)/float(p)
      PA = P - A
      PB = P - B
      PC = P - C
      PC2= np.sum(np.square(PC))
      x  = np.float(p*PC2)
      #Exs = E(m[0],n[0],0,p2,mu,AB[0],PA[0],PB[0])
      #Eys = E(m[1],n[1],0,p2,mu,AB[1],PA[1],PB[1])
      #Ezs = E(m[2],n[2],0,p2,mu,AB[2],PA[2],PB[2])
      #Rab = np.sum(np.square(AB))
      for t in range(0, L0+1):
        Ex = E(m[0],n[0],t,p2,mu,AB[0],PA[0],PB[0])
        for u in range(0, L1+1):
          Ey = E(m[1],n[1],u,p2,mu,AB[1],PA[1],PB[1])
          for v in range(0, L2+1):
            # print E, R function format:
            #  def E(m,n,t,p2,mu,ABi,PAi,PBi)
            #  def R(t,u,v,n,PC,p2,x)
            Ez = E(m[2],n[2],v,p2,mu,AB[2],PA[2],PB[2])
            Rtuv=R(t,u,v,0,PC,p2,x)
            Vp += cij*(Ex*Ey*Ez*Rtuv)*(2*np.pi/float(p))

  # returning final integral
  return Zc*Vp
########### END of coulomb_Gij ###########

##############################
#  reading from perl output  #
##############################
crt_i = 0
for line in open("COORD.pltmp", 'r'):
  if(crt_i == 0):
    COORD = np.atleast_2d(np.array(map(np.float, line.split( ))))
    crt_i += 1
  else:
    crd = np.array(map(np.float, line.split( )))
    COORD = np.atleast_2d(np.vstack([COORD, crd]))

lm_i = 0
for line in open("lmMAP.pltmp", 'r'):
  if(lm_i == 0):
    lmMAP = np.array(map(int, line.split( )))
    lm_i += 1
  else:
    crd = np.array(map(int, line.split( )))
    lmMAP = np.vstack([lmMAP, crd])



# CMNB format:
# charge, multiplicity, Ne, N-alpha, N-beta, N-AObasis
CMNB   = np.loadtxt('CMNBasis.pltmp', delimiter='\n').astype(int)
ACHARGE= np.atleast_1d(np.loadtxt( 'ACHARGE.pltmp', delimiter='\n').astype(int))
AGMAP  = np.loadtxt(   'AGMAP.pltmp', delimiter='\n').astype(int)
AOMAP  = np.loadtxt(   'AOMAP.pltmp', delimiter='\n').astype(int)
TYPE   = np.loadtxt(    'TYPE.pltmp', delimiter='\n').astype(int)
NSHELL = np.loadtxt(  'NSHELL.pltmp', delimiter='\n').astype(int)
MAP    = np.loadtxt(     'MAP.pltmp', delimiter='\n').astype(int)
PEXP   = np.loadtxt(    'PEXP.pltmp', delimiter='\n')
CTRC   = np.loadtxt(    'CTRC.pltmp', delimiter='\n')
EW_a   = np.loadtxt(    'EW_a.pltmp', delimiter='\n')
EV_a   = np.loadtxt(    'EV_a.pltmp', delimiter='\n')

ECPFlag=0
ECPPath='./ECPList.pltmp'
if os.path.isfile(ECPPath) and os.access(ECPPath, os.R_OK):
  ECPFlag=1

  talk= sys.argv[1].replace(".crd",".alk")
  texp= sys.argv[1].replace(".crd",".exp")
  tkf = sys.argv[1].replace(".crd",".kf")
  tki = sys.argv[1].replace(".crd",".ki")
  tlmx= sys.argv[1].replace(".crd",".lmx")
  tppl= sys.argv[1].replace(".crd",".ppl")
  tza = sys.argv[1].replace(".crd",".za")

  TEXP = np.loadtxt(texp, delimiter='\n')
  TAlk = np.loadtxt(talk, delimiter='\n')
  TKi  = np.loadtxt(tki, delimiter='\n')
  TKf  = np.loadtxt(tkf, delimiter='\n')
  TPPL = np.loadtxt(tlmx, delimiter='\n')
  TLMX = np.loadtxt(tppl, delimiter='\n')

  alk    = np.loadtxt(    'ECPa.pltmp', delimiter='\n')
  Alk    = np.loadtxt(    'ECPA.pltmp', delimiter='\n')
  Ki     = np.loadtxt(   'ECPKi.pltmp', delimiter='\n')
  Kf     = np.loadtxt(   'ECPKf.pltmp', delimiter='\n')
  PPSkip = np.loadtxt( 'ECPList.pltmp', delimiter='\n')
  LMax   = np.loadtxt( 'ECPLMax.pltmp', delimiter='\n')
  RZa    = np.loadtxt(   'ECPZa.pltmp', delimiter='\n')


############### END of reading from perl output ############


##################
#  system check  #
##################
NATOMS = np.array(ACHARGE).size
if (len(EV_a)/len(EW_a) == len(EW_a)):
  EVa = np.array(EV_a).reshape(len(EW_a),len(EW_a))
else:
  sys.exit("degenerated AO's is not yet implemented!")

spherical=0
s2c=[[1]]
for i in range(1,len(TYPE)):
  if(TYPE[i]>=0):
    nf = factorial(2+TYPE[i])/(factorial(2)*factorial(TYPE[i]))
    blk = np.eye(nf)
  else:
    sys.exit("spherical harmonics is not yet implemented!")
    spherical=1
    if(TYPE[i]==-1):
      sys.exit("Gaussian SP type orbitals are not yet implemented!")
    elif(TYPE[i]==-2):
      # WARNING: specific format for Gaussian 09!
      # Q: spherical to cartisian transformation?
      blk = np.transpose(np.array([
      [-0.315392,-0.315392, 0.630783, 0.0     , 0.0     , 0.0     ],\
      [ 0.0     , 0.0     , 0.0     , 0.0     , 1.092548, 0.0     ],\
      [ 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 1.092548],\
      [ 0.546274,-0.546274, 0.0     , 0.0     , 0.0     , 0.0     ],\
      [ 0.0     , 0.0     , 0.0     , 1.092548, 0.0     , 0.0     ]
     ]))

  s2c = blkdiag(s2c, blk)

if(spherical==1):
  EVa = np.dot(s2c, EVa)
  (NAO, tmp) = s2c.shape
else:
  NAO = len(EW_a)



HOMO = CMNB[3]-1
LUMO = CMNB[3]

ma  = int(len(EW_a) - LUMO)
tcrd= sys.argv[1]


##################################
#  construct coordinate systems  #
##################################

# readin target coordinates and charge #
if ECPFlag==0:
  crt_i = 0
  for line in open(sys.argv[1], 'r'):
    if(crt_i == 0):
      TCRD = np.atleast_2d(np.array(map(np.float, line.split( )))[1:4])
      TATOM= np.atleast_1d(np.array(map(np.float, line.split( )))[0])
      crt_i += 1
    else:
      crd = np.array(map(np.float, line.split( )))
      TCRD = np.atleast_2d(np.vstack([TCRD, crd[1:4]]))
      TATOM = np.vstack([TATOM, crd[0]])
else:
  crt_i = 0
  for line in open(sys.argv[1], 'r'):
    if(crt_i == 0):
      TCRD = np.atleast_2d(np.array(map(np.float, line.split( )))[1:4])
      TATOM= np.atleast_1d(np.array(map(np.float, line.split( )))[0])
      crt_i += 1
    else:
      crd = np.array(map(np.float, line.split( )))
      TCRD = np.atleast_2d(np.vstack([TCRD, crd[1:4]]))
      TATOM = np.vstack([TATOM, crd[0]])

  crt_i = 0
  for line in open(tza, 'r'):
    if(crt_i == 0):
      TZa= np.atleast_1d(np.array(map(np.float, line.split( )))[0])
      crt_i += 1
    else:
      crd = np.array(map(np.float, line.split( )))
      TZa = np.vstack([TZa, crd[0]])

TCHARGE=TATOM.tolist()

# check all the atoms with same coordinates and     #
# construct list of atoms with different coordinate #
match_i=0
match=[]
pcrd_i=0
for ri in range(0,len(ACHARGE)):
  for tj in range(0,len(TCHARGE)):
    if (COORD[ri]==TCRD[tj]).all():
      if(match_i==0):
        match=[ri,tj]
        match_i += 1
      else:
        match.append(ri)
        match.append(tj)
      if ECPFlag==0:
        if (ACHARGE[ri] != TATOM[tj]):
          if(pcrd_i==0):
            PCRD = np.atleast_2d(
                   np.hstack([ACHARGE[ri]-TATOM[tj], COORD[ri]]))
            pcrd_i += 1
          else:
            pcrd = np.atleast_2d(
                   np.hstack([ACHARGE[ri]-TATOM[tj], COORD[ri]]))
            PCRD = np.vstack([pcrd,PCRD])
      # compare nuclear number instead of charge #
      else:
        if (RZa[ri] != TZa[tj]):
          if(pcrd_i==0):
            PCRD = np.atleast_2d(
                   np.hstack([ACHARGE[ri]-TATOM[tj], COORD[ri]]))
            PPP = np.array([ri, tj])
            pcrd_i += 1
          else:
            pcrd = np.atleast_2d(
                   np.hstack([ACHARGE[ri]-TATOM[tj], COORD[ri]]))
            PCRD = np.vstack([pcrd,PCRD])
            PPP = np.vstack([[ri,tj],[PPP]])


# cat all reference atoms with different coordinate #
mi=0
for ri in range(0,len(ACHARGE)):
  if(mi < len(match)/2):
    if (ri != match[2*mi]):
      if(pcrd_i == 0):
        PCRD = np.atleast_2d(np.hstack([ACHARGE[ri], COORD[ri]]))
        if ECPFlag==1:
          PPP = np.array([ri, -1])
        pcrd_i += 1
      else:
        pcrd = np.hstack([ACHARGE[ri], COORD[ri]])
        PCRD = np.atleast_2d(np.vstack([pcrd,PCRD]))
        if ECPFlag==1:
          PPP = np.vstack([PPP, [ri, -1]])
    else:
      mi += 1
  else:
    if(pcrd_i == 0):
      PCRD = np.atleast_2d(np.hstack([ACHARGE[ri], COORD[ri]]))
      if ECPFlag==1:
        PPP = np.array([ri, -1])
      pcrd_i += 1
    else:
      pcrd = np.hstack([ACHARGE[ri], COORD[ri]])
      PCRD = np.atleast_2d(np.vstack([pcrd,PCRD]))
      if ECPFlag==1:
        PPP = np.vstack([PPP, [ri, -1]])

# cat all target atoms with different coordinate #
mi=0
for tj in range(0,len(TCHARGE)):
  if(mi < len(match)/2):
    if (tj != match[2*mi+1]):
      if(pcrd_i == 0):
        PCRD = np.atleast_2d(np.hstack([-TATOM[tj], TCRD[tj]]))
        if ECPFlag==1:
          PPP = np.array([-1, tj])
        pcrd_i += 1
      else:
        pcrd = np.hstack([-TATOM[tj], TCRD[tj]])
        PCRD = np.atleast_2d(np.vstack([pcrd,PCRD]))
        if ECPFlag==1:
          PPP = np.vstack([PPP, [-1,tj]])
    else:
      mi += 1
  else:
    if(pcrd_i == 0):
      PCRD = np.atleast_2d(np.hstack([-TATOM[tj], TCRD[tj]]))
      if ECPFlag==1:
        PPP = np.array([-1,tj])
      pcrd_i += 1
    else:
      pcrd = np.hstack([-TATOM[tj], TCRD[tj]])
      PCRD = np.atleast_2d(np.vstack([pcrd,PCRD]))
      if ECPFlag==1:
        PPP = np.vstack([PPP, [-1,tj]])

if(pcrd_i == 0):
  Npcrd = 0
else:
  Npcrd=len(PCRD[:,0])

print PPP

###############################
#  Calculate integral matrix  #
###############################
if ECPFlag==0:
  itrG = 0
  for a in range(0,Npcrd):
    new0 = []
    for i in range(0,NAO):
      new1 = []
      for j in range(0,NAO):
        # coulomb_Gij format
        #  def coulomb_Gij(ai,aj,C,Zc)
        gij_a = coulomb_Gij(i,j,PCRD[a][1:4],int(PCRD[a][0]))
        new1.append(gij_a);
      new0.append(new1);
    if(itrG == 0):
      G = np.array(new0)
      itrG += 1
    else:
      G += np.array(new0)

else:
  itrG = 0
  itrP = 0
  for a in range(0,Npcrd):
    ppi = PPP[a][0]
    ppf = PPP[a][1]
    new0 = []
    pp0  = []
    for i in range(0,NAO):
      new1 = []
      pp1  = []
      for j in range(0,NAO):

        ############################################
        #  Coulomb integral with effective charge  #
        ############################################
        # coulomb_Gij format
        #  def coulomb_Gij(ai,aj,C,Zc)
        gij_a = coulomb_Gij(i,j,PCRD[a][1:4],int(PCRD[a][0]))
        new1.append(gij_a);

        ########################################
        #  Pseudo integral of reference atoms  #
        ########################################
        ppl = 0
        if ppi>0:
          for l in range(0,int(LMax[ppi])):
            ki = int(Ki[3*l + ppi])
            kf = int(Kf[3*l + ppi])
            for k in range(ki,kf+1):
              #print k,l,alk[k]
              #<AOi|Alk*Exp[-alk*rI]|AOj>
              ppl -= Alk[k]*PP_ij(i,j,alk[k],ppi)
              #spherical projector is not yet implemented
              #<AOi|Alk*Exp[-alk*rI] (sum|lm><lm|) |AOj>
          #pp1.append(ppl)
          #print pp1

        #####################################
        #  Pseudo integral of target atoms  #
        #####################################
        if ppf>0:
          for l in range(0,int(TLMX[ppf])):
            ki = int(TKi[3*l + ppf])
            kf = int(TKf[3*l + ppf])
            for k in range(ki,kf+1):
              #print k,l,alk[k]
              #<AOi|Alk*Exp[-alk*rI]|AOj>
              ppl += TAlk[k]*PP_ij(i,j,TEXP[k],ppf)
              #spherical projector is not yet implemented
              #<AOi|Alk*Exp[-alk*rI] (sum|lm><lm|) |AOj>

        pp1.append(ppl)
      new0.append(new1)
      if ppi>0:
        pp0.append(pp1)
    if(itrG == 0):
      G = np.array(new0)
      itrG += 1
    else:
      G += np.array(new0)

    if ppi>0 or ppf>0:
      if itrP == 0:
        Gpp = np.array(pp0)
        print Gpp.shape
        itrP += 1
      else:
        Gpp += np.array(pp0)
        print Gpp.shape
  G += Gpp

#print G



#  Gaij.append(new0);
#G = np.array(Gaij)

#new0 = []
#for i in range(0,NAO):
#  new1 = []
#  for j in range(0,NAO):
#    new1.append(PP_ij(i,j,0,1))
#  new0.append(new1)
#PP = np.array(new0)
#print PP

if(Npcrd > 0):
  Gmo = np.dot(np.dot(EVa,G),np.transpose(EVa))
else:
  Gmo = np.zeros((NAO,NAO))
GMOout = Gmo.reshape(Gmo.size)
gmoOut = sys.argv[1].replace(".crd",".gmo")
gmoout = open(gmoOut,'w')
for i in range(0,Gmo.size):
#  gstring = str(Gout[i]) + "\n"
#  gout.write(str(gstring))
  print >>gmoout, '{0:14.7E}'.format(GMOout[i])
gmoout.close
