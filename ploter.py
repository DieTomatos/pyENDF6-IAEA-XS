import ENDF6
import matplotlib.pyplot as plt

#f = open('g_4-Be-9_0425.endf')
f = open('g_6-C-12_0625.endf')

#f = open('g_39-Y-89_3925.endf')
lines = f.readlines()
sec = ENDF6.find_section(lines, MF=3, MT=3)  # total cross-section
x, y = ENDF6.read_table(sec)

fig=plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, y)
print(y)
ax.set_xlabel('Photon energy [eV]')
ax.set_ylabel('Cross-section [barn]')
plt.show()
