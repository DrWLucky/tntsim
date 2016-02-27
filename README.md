//------------------------------------------------------------------
//------------------------------------------------------------------
// 
// Margot neutron detector simulation for GEANT4 (version 4.8.2)
//
// First Created 14/02/07 - Brian Roeder, LPC Caen
//                        - email - roeder@lpccaen.in2p3.fr
//
// Variation on the original Serra simulation for a Demon module
// Switches from using "Primitive Scorers" for Detector sensitivity
// to "Sensitive Detectors" that measure and output "Hit" positions
// 
// Envisioned as a first step toward multi-detector - tracking sims.
//
//
//------------------------------------------------------------------
// History
//-----------------------------------------------------------------
//
// - 14/02/07 - copied files from Serra and switched names over.
//            - started going to "Sensitive Detector" Readout