
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
//      File name:     IAEA_combiner.cc
//
//      Author:        Kutsenko Bogdan
//
//      Creation date: 14 July 2021
//
//      Modifications:
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <iostream>
#include <experimental/filesystem>

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
#include "G4Exp.hh"
#include "G4Log.hh"

#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"

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
#include "G4PhotoNuclearCrossSection.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

std::string path = "data/";

G4double stepEval(std::multimap <int, int> & ZA, G4int Z){
  G4double step=0.5;
  G4PhysicsVector* v = nullptr;
  std::ostringstream inname,outname;
  for (auto itr = ZA.begin(); itr != ZA.end(); itr++){
    if (itr -> first == Z){
      inname.clear();
       inname.str("");
       inname << path<<"inel"<< itr -> first << "_" << itr -> second;
       std::ifstream filein;
       filein.open(inname.str().c_str());
      v = new G4PhysicsVector();
      if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << inname.str().c_str()
	 << "> is not retrieved!";
    }
      filein.close();
    G4int size = v->GetVectorLength();
      for(G4int i = 1; i < size;i++) {
	G4double Ei=v->Energy(i);
	G4double Ep=v->Energy(i-1);
	if(step>(Ei-Ep) && (Ei-Ep)>0 && v->Value(Ep)!=0) step = Ei-Ep;
      }
    }
  }
  
  return step;
}

void vectorUntilBound(G4PhysicsVector* v ,std::vector <G4double>& cs, std::vector <G4double>& e, const G4double lTransitionBound){
  //DONE Copy G4physics vector values to double c++ vectors until the lTransitionBound
  G4int Npoints= v->GetVectorLength();
  G4double ei;
  cs.clear();
  e.clear();
  for(G4int i = 0; i < Npoints;i++){
    ei = v->Energy(i);
    if(ei>lTransitionBound && e.back()<lTransitionBound) {
    e.push_back(lTransitionBound);
    cs.push_back(v->Value(lTransitionBound));    
    break;
    }
    else if(ei>lTransitionBound) break;
    e.push_back(ei);
    cs.push_back(v->Value(ei));
  }
}


void transitionRegionCHIPS(std::vector <G4double>& cs, std::vector <G4double>& e, const G4double lTransitionBound,const G4double rTransitionBound, G4int Z, const G4double ChipsGGtransitionBound){
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  G4DynamicParticle dParticle(G4Gamma::Gamma(),aDirection,rTransitionBound);
  G4VCrossSectionDataSet* ggXsection = new G4PhotoNuclearCrossSection();
  //  G4Element* el = G4NistManager::Instance()->FindOrBuildElement(Z);
  //const G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_"+el->GetName());
  G4double rxs = ggXsection->GetElementCrossSection(&dParticle, Z, 0);
  G4double lxs = cs.back();
  G4double xs=0;
  
  for(G4double ei = lTransitionBound+1; ei <= rTransitionBound ; ei++){
      xs = lxs + (ei - lTransitionBound)*(rxs - lxs)/(rTransitionBound-lTransitionBound);
      e.push_back(ei);
      cs.push_back(xs);
  }
  if(rTransitionBound!=lTransitionBound) cs.push_back(rxs);
  for(G4double ei = rTransitionBound+1; ei <= ChipsGGtransitionBound; ei++){
    e.push_back(ei);
    dParticle.SetKineticEnergy(ei);
    cs.push_back(ggXsection->GetElementCrossSection(&dParticle, Z, 0));
  }
  
}

void transitionRegionISOCHIPS(std::vector <G4double>& cs, std::vector <G4double>& e, const G4double lTransitionBound,const G4double rTransitionBound, G4int Z, G4int A, const G4double ChipsGGtransitionBound){
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  G4DynamicParticle dParticle(G4Gamma::Gamma(),aDirection,rTransitionBound);
  G4VCrossSectionDataSet* ggXsection = new G4PhotoNuclearCrossSection();
  //  G4Element* el = G4NistManager::Instance()->FindOrBuildElement(Z);
  //const G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_"+el->GetName());
  G4double rxs = ggXsection->GetElementCrossSection(&dParticle, Z, 0);
  G4double lxs = cs.back();
  G4double xs=0;
  for(G4double ei = lTransitionBound+1; ei <= rTransitionBound ;ei++){
      xs = lxs + (ei - lTransitionBound)*(rxs - lxs)/(rTransitionBound-lTransitionBound);
      e.push_back(ei);
      cs.push_back(xs);
  }
  if(rTransitionBound!=lTransitionBound) cs.push_back(rxs);
  for(G4double ei = rTransitionBound+1; ei <= ChipsGGtransitionBound; ei++){
    e.push_back(ei);
    dParticle.SetKineticEnergy(ei);
    cs.push_back(ggXsection->GetIsoCrossSection(&dParticle, Z, A));
  }
  
}

void clearZero(std::vector <G4double>& cs, std::vector <G4double>& e, G4int isIsotope){

  G4double eraseNum = std::distance(cs.begin(), std::find_if(cs.begin(),cs.end(), [](G4double x) { return x != 0; }));
  cs.erase(cs.begin(),cs.begin()+eraseNum-isIsotope);
  e.erase(e.begin(),e.begin()+eraseNum-isIsotope);
  }

void writeFile(std::vector <G4double>& cs, std::vector <G4double>& e, std::ostringstream& outname){
  std::ofstream fileout;
  G4int esize = e.size();
  fileout.open(outname.str().c_str());       
  fileout<< e[0]<<" "<< e.back()<<" "<<esize<<"\n";
  fileout << esize<<"\n";
  for(G4int i = 0; i<esize; i++) fileout<<e[i]<<" "<<cs[i]<<"\n";
   
  fileout.close();
	   
  }

int main()//int argc, char** argv)
{
  /*
  // Control on input
  if(argc < 3) {
    G4cout << "Input parameters are not specified! Exit" << G4endl;
    return 1;
  }
  */
  G4cout << "=======================================================" << G4endl;
  G4cout << "======   Hadron Nucleon IAEA elementary cross-section combiner ========" << G4endl;
  G4cout << "=======================================================" << G4endl;

  // ------- Initialisation
  const G4double lTransitionBound = 130.*MeV;
  const G4double rTransitionBound = 130.*MeV;
  const G4double ChipsGGtransitionBound = 130.*MeV;
  const G4double tAbundTreshold = 11e-3; 
  //const G4double tAbundTreshold = 11e-3; // C, O, N, Ar, V, La, Ta, U - An elementary cross-section is taken as a copy of the most common isotope cross-section
  const G4int freeVectorException[] = {4, 6, 7, 8, 27, 39, 45, 65, 67, 69, 73};
  G4int flagLighIsotopeGenerate = 0;
  G4int flagCHIPSdata = 1;
  
  G4int verboseLevel = 0;
  std::ostringstream inname,outname;
  G4PhysicsVector* v = nullptr;
  G4VCrossSectionDataSet* inel = nullptr;
  std::string Zoutput, Aoutput, token, delimeter = "_";
  size_t pos = 0;
  std::multimap <int, int> ZA;
  G4int Z,Zmax=0;
  G4double totalAbundancy;
  
  
  G4NistElementBuilder* builder = new G4NistElementBuilder(0);
  
  for (const auto & entry : std::experimental::filesystem::directory_iterator(path)){
    Zoutput = std::regex_replace(entry.path().string() , std::regex("[^0-9]*([0-9]+).*"),std::string("$1"));       
  
    Aoutput = entry.path().string(); 
    while((pos = Aoutput.find(delimeter)) != std::string::npos){
      token = Aoutput.substr(0,pos);
      Aoutput.erase(0,pos+delimeter.length());
    }
    if(Zmax < std::stoi(Zoutput)) Zmax = std::stoi(Zoutput);
    ZA.insert(std::make_pair(std::stoi(Zoutput),std::stoi(Aoutput)));
  }
  
  for(Z = 1;Z<Zmax+1; Z++){
    G4int CHIPSwrite = flagCHIPSdata;
    totalAbundancy = 0.;
    G4int Aonly = 0;
    G4double isoAbund;
    G4double Hstep = 0.1;
    G4int  npoints;
    G4double step = stepEval(ZA,Z);
    if(Z==1 && ChipsGGtransitionBound>144.5*MeV ) npoints = round((ChipsGGtransitionBound-144.5*MeV)/Hstep)+1; // 144.5 MeV is a delta resonance treshold
    else{
    if(step<(0.001*MeV)) step = 0.001*MeV;
    npoints = round(lTransitionBound/step)+1;
    }
    //G4cout<<"step "<<stepEval(ZA,Z)<<" npoints "<<npoints<<G4endl;
    std::vector <G4double> cs;
    std::fill_n(std::back_inserter(cs), npoints, 0.);
    std::vector <G4double> e;
    std::fill_n(std::back_inserter(e), npoints, 0.);
    for (auto itr = ZA.begin(); itr != ZA.end(); itr++){
    if (itr -> first == Z && Z != 1){
      if(verboseLevel>0) G4cout << "Z = " << itr -> first << "  A = "<< itr -> second <<"  Abundancy = "<<builder->GetIsotopeAbundance(itr->first, itr->second)<< G4endl;
       inname.clear();
       inname.str("");
       inname << path<<"inel"<< itr -> first << "_" << itr -> second;
       std::ifstream filein;
       filein.open(inname.str().c_str());
       
     if (!(filein)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << inname.str().c_str()
	 << "> is not opened!";
     } else {
       if(verboseLevel > 0) {
	 G4cout << "File " << inname.str() 
	     << " is opened by IAEA_combiner" << G4endl;
       }
     }
      v = new G4PhysicsVector();
      if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << inname.str().c_str()
	 << "> is not retrieved!";
    }
      filein.close();
      isoAbund = builder->GetIsotopeAbundance(itr->first,itr->second);
      if(std::find(std::begin(freeVectorException), std::end(freeVectorException), Z) != std::end(freeVectorException) && isoAbund > 1 - tAbundTreshold) {
	Aonly = itr -> second;
      }
      totalAbundancy+=isoAbund;
      for(G4int i=0;i<npoints; i++ ){
	cs[i] += (v->Value(i*step*MeV)*builder->GetIsotopeAbundance(itr->first,itr->second));
	e[i] = (i*step);
      }

    }
    }

 
    
 auto itr= ZA.find(Z);
 itr++; 
 /*if(Z!=itr->first && totalAbundancy >= tAbundTreshold ) {// The case when for the element we have only one isotop file at the database
   std::ifstream filein;
   filein.open(inname.str().c_str());
   v = new G4PhysicsVector();
   v->Retrieve(filein, true);
   filein.close();
   vectorUntilBound(v, cs, e, lTransitionBound);
   G4cout<<Z<<" step = "<<step<<G4endl;
   }*/
 if(Aonly!=0){ // Aonly is the case when one isotope with hig habundancy Abund > 1-tAbundTreshold taken as element cross-section with non linear parametrisation due to narrow resonances
   //G4NistManager::Instance()->PrintElement(Z);
   inname.clear();
   inname.str("");
   inname << path<<"inel"<< Z << "_" << Aonly;
   std::ifstream filein;
   filein.open(inname.str().c_str());
   v = new G4PhysicsVector();
   v->Retrieve(filein, true);
   filein.close();
   vectorUntilBound(v, cs, e, lTransitionBound);

 }
  else if(Z==1){ // The proton case // Hstep = 0.1 MeV
    
    inel = new G4PhotoNuclearCrossSection();
    G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
    G4DynamicParticle dParticle(G4Gamma::Gamma(),aDirection,0);
    if(ChipsGGtransitionBound > 144.5){
    for(G4int i=0.;i<npoints; i++ ){  
      dParticle.SetKineticEnergy(144.5*MeV+i*Hstep);
      e[i]=144.5*MeV+i*Hstep;
      cs[i]=inel->GetElementCrossSection(&dParticle,Z,0);
      }

    if(ChipsGGtransitionBound != rTransitionBound){
    G4double logMinE=G4Log(ChipsGGtransitionBound+Hstep);
    G4double logMaxE=G4Log(1e+4);
    G4double dlE = (logMaxE-logMinE)/(1.e+02);
    
    for(G4double lei = logMinE; lei < logMaxE;lei += dlE){
      dParticle.SetKineticEnergy(G4Exp(lei));
      e.push_back(G4Exp(lei));
      cs.push_back(inel->GetElementCrossSection(&dParticle,Z,0));
    }
    }
    }
    else{
      e.push_back(lTransitionBound);
      dParticle.SetKineticEnergy(lTransitionBound);
      cs.push_back(inel->GetElementCrossSection(&dParticle,Z,0));
      }
    CHIPSwrite=CHIPSwrite - 1;
    }
 
 else if(totalAbundancy < tAbundTreshold){ // The case if we dont have sufficient isotope data  
 //DONE Implement CHIPS model if there is no IAEA data
    //bdkq// For CHIPS linear parametrisation used with 0.5 step
   //G4NistManager::Instance()->PrintElement(Z);
   inel = new G4PhotoNuclearCrossSection();
    G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
    G4DynamicParticle dParticle(G4Gamma::Gamma(),aDirection,0);
     for(G4int i=0;i<npoints; i++ ){  
      dParticle.SetKineticEnergy(i*step);
      e[i]=i*step;
      cs[i]=inel->GetElementCrossSection(&dParticle,Z,0);
      }

     CHIPSwrite=CHIPSwrite - 1;
 }
 else {// Case when we have isotope data with total abundancy more than tAbunanceTreshold presented at IAEA data library
   
     for(G4int i=0;i<npoints; i++ ){
	cs[i] = cs[i]/totalAbundancy;
      }
 }

 transitionRegionCHIPS(cs, e, lTransitionBound,lTransitionBound,Z,ChipsGGtransitionBound);//No need for transition region here, only CHIPS model is used
 
 outname.clear();
 outname.str("");
 outname << "combined_Data/inel"<< Z;
 clearZero(cs, e, 0);
 if(Z == 1 && ChipsGGtransitionBound<144.5*MeV){
   cs.clear();
   e.clear();
   e.push_back(0.);
   cs.push_back(0.);
   e.push_back(lTransitionBound);
   cs.push_back(0.);   
 }
 else{
  e.insert(e.begin(),e.front()-step);
  cs.insert(cs.begin(),0);
 }
 if(CHIPSwrite >= 0) writeFile(cs ,e , outname);
  }
 
 
//Next block is relevant to add transition region and CHIPS data to all isotopes
for (auto itr = ZA.begin(); itr != ZA.end(); itr++){
  std::vector <G4double> cs;
  std::vector <G4double> e;
  Z = itr -> first;
  G4int A = itr -> second;
  inname.clear();
  inname.str("");
  inname << path<<"inel"<< Z << "_" << A;
  std::ifstream filein;
  filein.open(inname.str().c_str());
       
  if (!(filein)) {
   G4ExceptionDescription ed;
   ed << "Data file <" << inname.str().c_str()
      << "> is not opened!";
     } else {
       if(verboseLevel > 0) {
	 G4cout << "File " << inname.str() 
	     << " is opened by IAEA_combiner" << G4endl;
       }
     }
  v = new G4PhysicsVector();
  if(!v->Retrieve(filein, true)) {
  G4ExceptionDescription ed;
  ed << "Data file <" << inname.str().c_str()
     << "> is not retrieved!";
    }
  filein.close();
  vectorUntilBound(v, cs, e, lTransitionBound);
    
  outname.clear();
  outname.str("");
  outname << "combined_Data/inel"<< Z <<"_"<<A;
  transitionRegionISOCHIPS(cs, e, lTransitionBound, rTransitionBound, Z, A , ChipsGGtransitionBound);
  auto upitr = itr;
  auto ditr = itr;
  upitr++;
  ditr--;
  clearZero(cs,e, 1);
  if(ditr->first == itr->first || upitr->first == itr->first) writeFile(cs , e , outname);  // No need to save isotope data if it is already saved as cross-section on element
  if(itr->first == 1 && itr->second == 2) writeFile(cs , e , outname);// H-2
  if(itr->first == 2 && itr->second == 3) writeFile(cs , e , outname);// He-3
  if(flagLighIsotopeGenerate == 1){
  if(itr->first == 1 && itr->second == 2) {
    inel = new G4PhotoNuclearCrossSection();
 G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
 G4DynamicParticle dParticle(G4Gamma::Gamma(),aDirection,0);

 G4double logMinE=G4Log(lTransitionBound);
 G4double logMaxE=G4Log(10.*GeV);
 G4double dlE = (logMaxE-logMinE)/(1.e+02);
    
 for(G4double lei = logMinE; lei <= logMaxE;lei += dlE){
     dParticle.SetKineticEnergy(G4Exp(lei));
     e.push_back(G4Exp(lei));
     cs.push_back(inel->GetIsoCrossSection(&dParticle,1,2));
    }
 clearZero(cs, e, 1);
 writeFile(cs , e , outname);// deutron
  }
  
 }
 }
 
 if(flagLighIsotopeGenerate == 1){
 std::vector <G4double> cs;
 std::vector <G4double> e;
 inel = new G4PhotoNuclearCrossSection();
 G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
 G4DynamicParticle dParticle(G4Gamma::Gamma(),aDirection,0);

 G4double logMinE=G4Log(5.*MeV);
 G4double logMaxE=G4Log(10.*GeV);
 G4double dlE = (logMaxE-logMinE)/(1.e+02);
    
 for(G4double lei = logMinE; lei <= logMaxE+dlE;lei += dlE){
     dParticle.SetKineticEnergy(G4Exp(lei));
     e.push_back(G4Exp(lei));
     cs.push_back(inel->GetIsoCrossSection(&dParticle,1,3));
    }

 outname.clear();
 outname.str("");
 outname << "combined_Data/inel"<< 1 <<"_" << 3;
 clearZero(cs, e, 1);
 writeFile(cs , e , outname);
 }

}
