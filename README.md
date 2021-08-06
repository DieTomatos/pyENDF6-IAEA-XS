# pyENDF6

This python module provides a minimal set of function for reading tabulated data from ENDF nuclear data files (see
https://www-nds.iaea.org/public/endf/).

#### Example
```python
import ENDF6
f = open('C006.endf')
lines = f.readlines()
sec = ENDF6.find_section(lines, MF=3, MT=3)  # total cross-section
x, y = ENDF6.read_table(sec)

figure()
plot(x, y)
xlabel('Photon energy [eV]')
ylabel('Cross-section [barn]')
show()
```

# Additions
 pyENDF6 was implemented to convert Photonuclear XS IAEA data files to txt files that can be uploaded with retrive method of G4Physics vector

# converterAllType.py
 Converts all types of file from IAEA database to txt format. Additional point with zero value is added to the left boundary of the cross section. Data taken from:
data/Type* 

# IAEA_combiner/IAEA_combiner.cc
 Create elementary cross-section from isotope XS data files from IAEA database. If there is no isotope data for the element the result is taken from the CHIPS model. 

The input is taken from IAEA_combiner/data
The output is placed at IAEA_combiner/combined_DATA

Adjustable parameters:
1) Linear connection implemented between IAEA photonuclear data and CHIPS model the bound of this linear can be changed: 
  const G4double lTransitionBound = 130.*MeV; 
  const G4double rTransitionBound = 150.*MeV;
2) In the case of Glauberg-Gribov model implementation additional bound is added to the CHIPS model. By default, it is not implemented and the value is equal to right liniar region bound:
  const G4double ChipsGGtransitionBound = 150.*MeV;
3) Some isotopes in the database don't have sufficient abundance to produce cross-section on the element. Threshold of total presented isotope abundances is adjustable:
  const G4double tAbundTreshold = 11e-3; 

# RetrievePlotter/RetrievePlotter.cc
Draw the resulted cross-section data from the txt


