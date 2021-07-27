import ENDF6
import matplotlib.pyplot as plt
import re,os
import numpy as np
import sys

treshold = 100000 #ev

def convert_endf(nType,nMF,nMT):
 folder = 'data/Type'+str(nType)+'/'
 for filename in os.listdir(folder):
  name=folder+filename
  f = open(name,"r")
  lines = f.readlines()
  sec = ENDF6.find_section(lines, MF=nMF, MT=nMT)  # total cross-section
  x, y = ENDF6.read_table(sec)
  if y[0]>sys.float_info.epsilon:
   y = np.insert(y,0,0.)
   x = np.insert(x,0,x[0]-treshold)
   
  ZA=[int(s) for s in re.split('_|-|',name) if s.isdigit()]
  outname='IAEA_combiner/data/'+'inel'+str(ZA[0])+'_'+str(ZA[1])
  f = open(outname,"w")
  f.write(str(x[1]/pow(10,6))+' '+ str(x[-1]/pow(10,6))+' '+ str(len(x))+'\n')
  f.write(str(len(x))+'\n')
  i=0
  while i < len(x):
   f.write(str(x[i]/pow(10,6))+' '+str(round(y[i]*1.e-22,34))+'\n')
   i+=1

def get_x_y_size_as(y,lines,nMF,nMT):
 sect = ENDF6.find_section(lines, MF=nMF, MT=nMT) 
 xt, yt = ENDF6.read_table(sect) 
 padded_yt=np.pad(yt,(len(y)-len(yt),0),'constant')
 return np.array(xt),np.array(padded_yt)
 
def convert_endf_typeIV():
 folder = 'data/Type4/'
 for filename in os.listdir(folder):
  name=folder+filename
  f = open(name,"r")
  lines = f.readlines()
  sec = ENDF6.find_section(lines, MF=3, MT=5)  # absorption cross-section
  x, y = ENDF6.read_table(sec)
  x=np.array(x)
  y=np.array(y)
  xt,yt = get_x_y_size_as(y,lines,3,201)
  y+=np.array(yt)
  xt,yt = get_x_y_size_as(y,lines,3,203)
  y+=np.array(yt)
  xt,yt = get_x_y_size_as(y,lines,3,204)
  y+=np.array(yt)
  xt,yt = get_x_y_size_as(y,lines,3,205)
  y+=np.array(yt)
  xt,yt = get_x_y_size_as(y,lines,3,206)
  y+=np.array(yt)

  if y[0]>sys.float_info.epsilon:
   y = np.insert(y,0,0.)
   x = np.insert(x,0,x[0]-treshold) 
  ZA=[int(s) for s in re.split('_|-|',name) if s.isdigit()]

  outname='IAEA_combiner/data/'+'inel'+str(ZA[0])+'_'+str(ZA[1])
 
  f = open(outname,"w")
  f.write(str(x[0]/pow(10,6))+' '+ str(x[-1]/pow(10,6))+' '+ str(len(x))+'\n')
  f.write(str(len(x))+'\n')
  i=0
  while i < len(x):
   f.write(str(x[i]/pow(10,6))+' '+str(round(y[i]*1.e-22,34))+'\n')
   i+=1

convert_endf(0,3,3)
convert_endf(1,3,5)
convert_endf(2,3,5)
convert_endf(3,3,3)
convert_endf_typeIV()
convert_endf(5,3,1)
