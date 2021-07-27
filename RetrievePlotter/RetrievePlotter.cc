
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     HadrHN.cc
//
//      Author:        V.Ivanchenko
//
//      Creation date: 14 August 2018
//
//      Modifications:
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"

#include "Histo.hh"
#include "G4Timer.hh"

#include "G4HadronNucleonXsc.hh"
#include "G4ComponentSAIDTotalXS.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"
#include "G4ChipsNeutronInelasticXS.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsProtonInelasticXS.hh"
#include "G4ChipsPionPlusElasticXS.hh"
#include "G4ChipsPionPlusInelasticXS.hh"
#include "G4ChipsPionMinusElasticXS.hh"
#include "G4ChipsPionMinusInelasticXS.hh"
#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonMinusElasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4EmParameters.hh"
#include "G4ComponentGGNuclNuclXsc.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

int main(int argc, char** argv)//int argc, char** argv)
{
  // Control on input
  if(argc < 4) {
    G4cout << "Input parameters are not specified!\nUsage : ./RetrievePlotter Z A atomName\nExit" << G4endl;
    return 1;
  }
  G4cout << "=======================================================" << G4endl;
  G4cout << "======   Hadron Nucleon IAEA elementary cross-section ========" << G4endl;
  G4cout << "=======================================================" << G4endl;

  // ------- Initialisation
  G4int verboseLevel = 0;
  std::ostringstream ost;
  G4int Z=atoi(argv[1]),A=atoi(argv[2]);
  G4String atomName=argv[3];
  if(A != 0) ost << "../data/output/inel"<< Z << "_" << A;			    
  else ost << "../IAEA_combiner/combined_Data/inel"<< Z;
  G4PhysicsVector* v = nullptr;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
  } else {
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by RetrievePlotter" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsVector();
    if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4RetrievePlotter::RetrieveVector(..)","had015",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  }
  Histo    histo;
  G4String hname;
  std::ostringstream ostname;
  if(A!=0) ostname << "test/gamma_"<<atomName<<"_"<<"IAEA"<<"_"<<Z<<"_"<<A;
  else ostname << "test/gamma_"<<atomName<<"_"<<"IAEA"<<"_"<<Z;
  hname = ostname.str();
  const G4int nbins = 481;
  G4double lxmin = 0*MeV;
  G4double lxmax = 200*MeV;
  G4double dlx = (lxmax - lxmin)/G4double(nbins-1);

  histo.Add1D("0","Inelastic",nbins,lxmin-dlx*0.5,lxmax+dlx*0.5);
  histo.SetFileName(hname);
  histo.Book();
  G4cout << "Histograms are booked output file <" << hname << "> "
	 << G4endl;
  G4double lxsi[nbins];
  G4double lx0 = lxmin;
  for (G4int i=0; i<nbins; ++i) {
    G4double lp = lx0;
    lx0 += dlx;
    G4double le = lp;
    lxsi[i] = v->Value(le)/CLHEP::millibarn;
    histo.Fill(0,lx0,lxsi[i]);
  }
  
  G4cout << "------------------------------------------------------------------------------------"
	 << G4endl;
  
  histo.Save();
}
