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
// $Id: TntScintSD.hh 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Tnt/include/TntScintSD.hh
/// \brief Definition of the TntScintSD class
//
//
#ifndef TntScintSD_h
#define TntScintSD_h 1

#include "TntScintHit.hh"
#include "TntDataRecordTree.hh"

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

class TntScintSD : public G4VSensitiveDetector
{
  public:

    TntScintSD(G4String name, G4String Light);
    virtual ~TntScintSD();

    virtual void Initialize(G4HCofThisEvent* );
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* );
    virtual void EndOfEvent(G4HCofThisEvent* );
    virtual void clear();
    virtual void DrawAll();
    virtual void PrintAll();

//by Shuya 160407
  G4double ConvertToLight(G4String theName, G4double theCharge, G4double edep, G4String Light);
 
  private:

    TntScintHitsCollection* fScintCollection;

//by Shuya 160407
    //TntScintHitsCollection* Tnt:HitsCollection;
    TntDataRecordTree* TntDataOutEV;

  G4String Light_Conv;
 
};

#endif
