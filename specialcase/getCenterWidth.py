#!/usr/bin/python
import numpy as np

f=open("DOS0","r")
lines=f.readlines()
base=np.float64([lines[i].split() for i in range(1,len(lines))])
fermi=abs(base[:,0]).argmin()
if abs(base[fermi][1])>1E-2:
    spin=1
else:
    spin=2

for element_idx in [2,3]:
  d_orbital=spin
  element=elements[element_idx],
  g=open("DOS"+str(element_idx),"r")
  lines=g.readlines()
  basic=np.float64([lines[i].split() for i in range(0,len(lines))])
  base=np.float64([[item[0],item[11]+item[13]+item[15]+item[17],abs(item[12]+item[14]+item[16]+item[18])] for item in basic])
  fermi=abs(base[:,0]).argmin()
  startingpoint_rev=np.argmax(abs(upperpart[:,d_orbital])>0.15)
  lowercount_rev=np.argmax(abs(upperpart[startingpoint+3:,d_orbital])<0.01)
  startingpoint=fermi-startingpoint_rev
  lowercount=fermi-lowercount_rev
  uppercount=np.argmax(abs(base[startingpoint:fermi,d_orbital])<0.01)
  width=base[uppercount,0]-base[lowercount,0]
  print width
  
