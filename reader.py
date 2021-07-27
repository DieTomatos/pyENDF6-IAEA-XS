import os,re
import numpy as np

def Read_Two_Column_File(file_name):
    with open(file_name, 'r') as data:
        x = []
        y = []
        start = 0 
        for i,line in enumerate(data):
         if i==0:
          t = line.split()
          start = float(t[0])
         elif i>1:
          p = line.split()
          x.append(float(p[0]))
          y.append(float(p[1]))
    return start,x, y

isab=([0.0759,0.9241],[1.],[1.],[1.],[1.],[1.],[1.],[0.7899],[1.],[0.922297],[0.9493],[0.7578,0.2422],[1.],[0.932581,0.067302],[0.96941],[1.],[0.0825,0.0744,0.7372,0.0541,0.0518],[1.],[0.83789],[1.],[0.05845, 0.91754, 0.02119, 0.00282],[1.],[0.680769,0.262231,0.011399,0.036345,0.009256],[0.6917,0.3083],[0.4863,0.2790,0.0410,0.1875,0.0062],[0.2084,0.2754,0.0773,0.3628,0.0761],[1.],[0.4961],[0.0056,0.0986,0.07,0.8258],[1.],[0.5145],[1.],[0.1484,0.0925,0.1592,0.1668,0.0955,0.2413,0.0963],[0.0187],[1.],[0.0102,0.1114,0.2233,0.2733,0.2646,0.1172],[0.51839,0.48161],[0.0125,0.0089,0.1249,0.128,0.2413,0.1222,0.2873,0.0749],[0.9571],[0.2422])
names=(['inel3_6','inel3_7'],['inel4_9'],['inel6_12'],['inel7_14'],['inel8_16'],['inel9_19'],['inel11_23'],['inel12_24'],['inel13_27'],['inel14_28'],['inel16_32'],['inel17_35','inel17_37'],['inel18_40'],['inel19_39','inel19_41'],['inel20_40'],['inel21_45'],['../inel22_46_test','inel22_47','inel22_48','inel22_49','inel22_50'],['inel23_51'],['inel24_52'],['inel25_55'],['inel26_54','inel26_56','inel26_57','inel26_58'],['inel27_59'],['inel28_58','inel28_60','inel28_61','inel28_62','inel28_64'], ['inel29_63','inel29_65'], ['inel30_64','../inel30_66_test','inel30_67','../inel30_68_test','inel30_70'],['inel32_70','inel32_72','inel32_73','inel32_74','inel32_76'],['inel33_75'],['inel34_80'],['inel38_84','inel38_86','inel38_87','inel38_88'],['inel39_89'],['inel40_90'],['inel41_93'],['inel42_92','inel42_94','inel42_95','inel42_96','inel42_97','inel42_98','inel42_100'],['inel44_98'],['inel45_103'],['inel46_102','inel46_104','inel46_105','inel46_106','inel46_108','inel46_110'],['inel47_107','inel47_109'],['inel48_106','inel48_108','inel48_110','inel48_111','inel48_112','inel48_113','inel48_114','inel48_116'],['inel49_115'],['inel50_118'])

for k in range(len(names)):
 start = []
 x = []
 y = []
 for i in range(len(isab[k])):
  tstart,tx, ty = Read_Two_Column_File('data/output/'+names[k][i]) 
  start.append(tstart)
  x.append(tx)
  y.append(ty)
 max_start = max(start)
 max_index = start.index(max_start)
 num = len(x[max_index])
 xout = np.zeros(num)
 yout = np.zeros(num)
 for j in range(len(isab[k])):
  st = int(round((max_start-start[j])/0.5))
  if st != 0:
   del x[j][0:st]
   del y[j][0:st]
 
 for i in range(num):
  xout[i] = x[max_index][i]
  for j in range(len(isab[k])):
   if abs(sum(isab[k])-1.0)>0.01:
    yout[i] += y[j][i]/isab[k][j]
   else:
    if xout[i] == x[j][i]:
     yout[i] += isab[k][j]*y[j][i]
    else:
     print("WARNING!! Different energy points!")
     print(xout[i], x[j][i])
 ZA=[int(s) for s in re.split('_|-|l',names[k][0]) if s.isdigit()]
 outname='data/output/'+'inel'+str(ZA[0])
 f = open(outname,"w")
 f.write(str(xout[1])+' '+ str(xout[-1])+' '+ str(len(xout))+'\n')
 f.write(str(len(xout))+'\n')
 i=0
 while i < len(xout):
  f.write(str(xout[i])+' '+str(round(yout[i],32))+'\n')
  i+=1

#if start2>start1:
 #st = (start2-start1)/0.5
 #for xi in x2:

