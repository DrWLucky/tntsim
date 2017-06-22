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
// $Id: TntPrimaryGeneratorAction.hh 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Tnt/include/TntPrimaryGeneratorAction.hh
/// \brief Definition of the TntPrimaryGeneratorAction class
//
//
#ifndef TntPrimaryGeneratorAction_h
#define TntPrimaryGeneratorAction_h 1
#include <memory>
#include "G4VUserPrimaryGeneratorAction.hh"

//by Shuya 160407
#include "TntDataRecordTree.hh"

#include "TntReactionGenerator.hh"
#include "TntBeamEmittance.hh"
#include "TntNeutronDecay.hh"

class G4ParticleGun;
class G4Event;

class TntPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	TntPrimaryGeneratorAction();
	virtual ~TntPrimaryGeneratorAction();
 
public:
	virtual void GeneratePrimaries(G4Event* anEvent);

protected:
//by Shuya 160407
	TntDataRecordTree* TntDataOutPG;

	G4ParticleGun* fParticleGun;
//by Shuya 160510
	G4String BeamType;
};



class TntPGAReaction : public TntPrimaryGeneratorAction {
public:
	TntPGAReaction();
	virtual ~TntPGAReaction();
	virtual void GeneratePrimaries(G4Event* anEvent);

protected:
	G4String fReacFile;
	G4int fA[4], fZ[4];
	std::auto_ptr<TntReactionGenerator> fReac;
	std::auto_ptr<TntNeutronDecay> fDecay;
	std::auto_ptr<TntRng> fRngEbeam, fRngEx3, fRngEx4, fRngTheta, fRngPhi;
	std::auto_ptr<TntBeamEmittance> fEmX, fEmY;
};

#endif
