#!/usr/bin/python
import os
import numpy as np
from subprocess import call
import re
import sys
import scipy
import scipy.optimize as optimize
import imp
#np.set_printoptions(threshold=np.nan,suppress=True)
SCRIPTDIR=os.path.dirname(os.path.realpath(__file__))
#generic read POSCAR

master_data_structure={}
f=open("POSCAR","r")
master_data_structure['lines']=f.readlines()
master_data_structure['base']=[master_data_structure['lines'][i].split()[0:3] for i in range(8,18)]
master_data_structure['cell']=[master_data_structure['lines'][i].split() for i in range(2,5)]
master_data_structure['element']=[]
for tmp_nelement in range(0,len(master_data_structure['lines'][5].split())):
    for tmp_natom_perelement in range(0,int(master_data_structure['lines'][6].split()[tmp_nelement])):
        master_data_structure['element'].append(master_data_structure['lines'][5].split()[tmp_nelement])
if len(master_data_structure['base'][-1])==0:
 print 'Warning: last line of POSCAR should not be empty, watch it! Removing the last line...'
 master_data_structure['base'].pop(-1)
if len(master_data_structure['base'][-1])==0:
 print 'Error: last line of POSCAR still empty! '
 exit(-1)
master_data_structure['base']=np.float64(master_data_structure['base'])
master_data_structure['cell']=np.float64(master_data_structure['cell'])
#2.image to supercell
master_data_structure['pos_imaged']=[]
master_data_structure['rpos_imaged']=[]
for i in [0,1,-1]:
 for j in [0,1,-1]:
  for k in [0,1,-1]:
   tmp_image_shift=np.float64([i,j,k])
   for id_base in range(0,len(master_data_structure['base'])):
    ele_pos_imaged=tmp_image_shift+master_data_structure['base'][id_base]
    master_data_structure['rpos_imaged'].append(ele_pos_imaged)
    ele_pos_imaged=np.dot(ele_pos_imaged,master_data_structure['cell'])
    master_data_structure['pos_imaged'].append(ele_pos_imaged)
master_data_structure['pos_imaged']=np.float64(master_data_structure['pos_imaged'])
master_data_structure['rpos_imaged']=np.float64(master_data_structure['rpos_imaged'])
master_data_structure['pos_original']=np.dot(master_data_structure['base'],master_data_structure['cell'])
master_data_structure['pos_original']=np.float64(master_data_structure['pos_original'])
#3.calculate distances
#calculate dist_qui
master_data_structure['dist_qui']=[]
changefromto=[]
for id1_base in range(0,len(master_data_structure['base'])):
 for id2_base in range(0,len(master_data_structure['base'])):
  for id2_pos_imaged in [id2_base+tmp_image_shift*len(master_data_structure['base']) for tmp_image_shift in range(0,27)]:
   tmp_dist=np.linalg.norm(master_data_structure['pos_original'][id1_base]-master_data_structure['pos_imaged'][id2_pos_imaged])
   if abs(tmp_dist)<0.1:
    continue
   #THE 'IF' JUDGEMENT THAT ADD 0.4 TO disTANCES WHEN NECESSARY. i IS THE INDEX OF FIRST ATOM, STARTING FROM 0. j IS THAT OF THE SECOND. k IS THE INDEX OF IMAGED SECOND ATOM.
   #if i>1 and j>1:
    #dist=dis`t+0.4
   id1_pos_imaged=id1_base
   master_data_structure['dist_qui'].append([id1_pos_imaged,id2_pos_imaged,id1_base,id2_base,tmp_dist])
master_data_structure['dist_qui']=np.float64(master_data_structure['dist_qui'])
master_data_structure['dist_qui']=master_data_structure['dist_qui'][np.argsort(master_data_structure['dist_qui'][:,4])]

#RB(C,O)
f=open("POSCAR","r")
elements=f.readlines()[5].split()
EleB=elements[1]
EleC=elements[2]
EleO=elements[3]
RB=float(os.popen(SCRIPTDIR+'/radius '+EleB).read())
RC=float(os.popen(SCRIPTDIR+'/radius '+EleC).read())
RO=float(os.popen(SCRIPTDIR+'/radius '+EleO).read())
#MAGMOM
magmoms=[abs(float(x)) for x in os.popen(SCRIPTDIR+'/mag').read().split()]
MAGMOM=max(magmoms)
#UB,UC
UB=float(os.popen(SCRIPTDIR+'/U '+EleB).read().split()[0])-float(os.popen(SCRIPTDIR+'/U '+EleB).read().split()[1])
UC=float(os.popen(SCRIPTDIR+'/U '+EleC).read().split()[0])-float(os.popen(SCRIPTDIR+'/U '+EleC).read().split()[1])
#StericBondParameterB(C)
rBO=float(master_data_structure['dist_qui'][0][4])
rCO=float(master_data_structure['dist_qui'][12][4])
StericBondParameterB=abs(rBO-RB-RO)/(RB+RO)
StericBondParameterC=abs(rCO-RC-RO)/(RC+RO)
#IsHalfMetal,HalfMetalGapWidth
tmp_spintronics=os.popen(SCRIPTDIR+'/spintronics').read().split()
IsHalfMetal=0
HalfMetalGapWidth=0
if tmp_spintronics[0] == 'HM':
    IsHalfMetal=1
    HalfMetalGapWidth=float(tmp_spintronics[1])
#CompoundName, NumberElectronBD(CD), NumberElectronB, NumberElectronC, AvalB, AvalC
tmpgen=os.popen(SCRIPTDIR+'/spec_edi/gen generic notex spin splitdos').read().split()
CompoundName=elements[0]+'2'+elements[1]+elements[2]
NumberElectronBD=float(tmpgen[8])
NumberElectronCD=float(tmpgen[10])
NumberElectronB=float(tmpgen[2])
NumberElectronC=float(tmpgen[3])
NumberElectronO=float(tmpgen[4])
AvalB=float(tmpgen[0])
AvalC=float(tmpgen[1])
#IsConventionalIntegerItyB(C),NumberElectronRealBD(CD)
NumberElectronRealCD=NumberElectronCD-(8-NumberElectronO)/15*5
NumberElectronRealBD=NumberElectronBD-(8-NumberElectronO)/15*5
IsConventionalIntegerItyB=abs(NumberElectronBD-np.floor(NumberElectronBD)-0.5)
IsConventionalIntegerItyC=abs(NumberElectronCD-np.floor(NumberElectronCD)-0.5)

print CompoundName,RB,RC,UB,UC,MAGMOM,StericBondParameterB,StericBondParameterC, NumberElectronBD, NumberElectronCD, NumberElectronB, NumberElectronC, AvalB, AvalC,IsConventionalIntegerItyB,IsConventionalIntegerItyC,NumberElectronRealBD,NumberElectronRealCD,IsHalfMetal,HalfMetalGapWidth,
#EnergyBD1,WidthBD1,EnergyBD2,WidthBD2,EnergyCD1,WidthCD1,EnergyCD2,WidthCD2,EnergyOD1,WidthOD1,EnergyOD2,WidthOD2,

#getCenterWidth
#s p_y p_z p_x d_xy d_yz d_z2 d_xz d_x^2-y^2
os.popen('~/src/qchem_scripts/vtstscripts/split_dos 2>/dev/null')
element_idx=3
for spin in [1,2]:
    g=open("DOS"+str(element_idx),"r")
    g.seek(0,0)
    lines=g.readlines()
    basic=np.float64([lines[i].split() for i in range(1,len(lines))])
    base=np.float64([[item[0],abs(item[3]+item[5]+item[7]),abs(item[4]+item[6]+item[8])] for item in basic])
    fermi=abs(base[:,0]).argmin()
    if abs(base[fermi][spin])>1E-1:
        tophalf=base[fermi:0:-1]
        bothalf=base[fermi:len(base)]
        top1=np.where(abs(tophalf[:,spin])<0.05)[0][0]
        bot1=np.where(abs(bothalf[:,spin])<0.05)[0][0]
        WidthBD=-tophalf[top1,0]+bothalf[bot1,0]
        EnergyBD=(tophalf[top1,0]+bothalf[bot1,0])/2
        print WidthBD,EnergyBD,
    else:
        d_orbital=spin
        upperpart=base[fermi:0:-1]
        startingpoint_rev=np.where(abs(upperpart[:,d_orbital])>0.3)[0][0]
        startingpoint=fermi-startingpoint_rev
        lowercount_rev=np.where(upperpart[startingpoint_rev:,spin]<0.1)[0][0]
        lowercount_rev=lowercount_rev+startingpoint_rev
        lowercount=fermi-lowercount_rev
        uppercount=np.where(base[startingpoint:fermi,d_orbital]<0.1)[0][0]
        uppercount=uppercount+startingpoint
        width=base[uppercount,0]-base[lowercount,0]
        WidthBD=base[uppercount,0]-base[lowercount,0]
        EnergyBD=(base[uppercount,0]+base[lowercount,0])/2
        print WidthBD,EnergyBD,
element_idx=4
for spin in [1,2]:
    g=open("DOS"+str(element_idx),"r")
    g.seek(0,0)
    lines=g.readlines()
    basic=np.float64([lines[i].split() for i in range(1,len(lines))])
    base=np.float64([[item[0],abs(item[3]+item[5]+item[7]),abs(item[4]+item[6]+item[8])] for item in basic])
    fermi=abs(base[:,0]).argmin()
    if abs(base[fermi][spin])>1E-1:
        tophalf=base[fermi:0:-1]
        bothalf=base[fermi:len(base)]
        top1=np.where(abs(tophalf[:,spin])<0.05)[0][0]
        bot1=np.where(abs(bothalf[:,spin])<0.05)[0][0]
        WidthCD=-tophalf[top1,0]+bothalf[bot1,0]
        EnergyCD=(tophalf[top1,0]+bothalf[bot1,0])/2
        print WidthCD,EnergyCD,
    else:
        d_orbital=spin
        upperpart=base[fermi:0:-1]
        startingpoint_rev=np.where(abs(upperpart[:,d_orbital])>0.3)[0][0]
        startingpoint=fermi-startingpoint_rev
        lowercount_rev=np.where(upperpart[startingpoint_rev:,spin]<0.1)[0][0]
        lowercount_rev=lowercount_rev+startingpoint_rev
        lowercount=fermi-lowercount_rev
        uppercount=np.where(base[startingpoint:fermi,d_orbital]<0.1)[0][0]
        uppercount=uppercount+startingpoint
        width=base[uppercount,0]-base[lowercount,0]
        WidthCD=base[uppercount,0]-base[lowercount,0]
        EnergyCD=(base[uppercount,0]+base[lowercount,0])/2
        print WidthCD,EnergyCD,
element_idx=6
for spin in [1,2]:
    g=open("DOS"+str(element_idx),"r")
    g.seek(0,0)
    lines=g.readlines()
    basic=np.float64([lines[i].split() for i in range(1,len(lines))])
    base=np.float64([[item[0],abs(item[3]+item[5]+item[7]),abs(item[4]+item[6]+item[8])] for item in basic])
    fermi=abs(base[:,0]).argmin()
    if abs(base[fermi][spin])>1E-1:
        tophalf=base[fermi:0:-1]
        bothalf=base[fermi:len(base)]
        top1=np.where(abs(tophalf[:,spin])<0.05)[0][0]
        bot1=np.where(abs(bothalf[:,spin])<0.05)[0][0]
        WidthOD=-tophalf[top1,0]+bothalf[bot1,0]
        EnergyOD=(tophalf[top1,0]+bothalf[bot1,0])/2
        print WidthOD,EnergyOD,
    else:
        d_orbital=spin
        upperpart=base[fermi:0:-1]
        startingpoint_rev=np.where(abs(upperpart[:,d_orbital])>0.3)[0][0]
        startingpoint=fermi-startingpoint_rev
        lowercount_rev=np.where(upperpart[startingpoint_rev:,spin]<0.1)[0][0]
        lowercount_rev=lowercount_rev+startingpoint_rev
        lowercount=fermi-lowercount_rev
        uppercount=np.where(base[startingpoint:fermi,d_orbital]<0.1)[0][0]
        uppercount=uppercount+startingpoint
        width=base[uppercount,0]-base[lowercount,0]
        WidthOD=base[uppercount,0]-base[lowercount,0]
        EnergyOD=(base[uppercount,0]+base[lowercount,0])/2
        print WidthOD,EnergyOD,
