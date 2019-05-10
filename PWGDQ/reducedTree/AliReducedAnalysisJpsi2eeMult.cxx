//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisJpsi2eeMult.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TIterator.h>
#include <TList.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"
#include <TRandom3.h>

ClassImp(AliReducedAnalysisJpsi2eeMult);


//___________________________________________________________________________
AliReducedAnalysisJpsi2eeMult::AliReducedAnalysisJpsi2eeMult() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
  fOptionRunMixing(kTRUE),
  fOptionRunRotation(kFALSE),
  fNRotations(1),
  fOptionRunPairing(kTRUE),
  fOptionRunOverMC(kFALSE),
  fOptionRunLikeSignPairing(kTRUE),
  fOptionLoopOverTracks(kTRUE),
  fOptionRunPrefilter(kTRUE),
  fOptionStoreJpsiCandidates(kFALSE),
  fEventCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts(),
  fPosTracks(),
  fNegTracks(),
  fPrefilterPosTracks(),
  fPrefilterNegTracks(),
  fJpsiCandidates(),
  fLegCandidatesMCcuts(),
  fJpsiMotherMCcuts(),
  fJpsiElectronMCcuts()
{
  //
  // default constructor
  //
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2eeMult::AliReducedAnalysisJpsi2eeMult(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
  fOptionRunMixing(kTRUE),
  fOptionRunRotation(kFALSE),
  fNRotations(1),
  fOptionRunPairing(kTRUE),
  fOptionRunOverMC(kFALSE),
  fOptionRunLikeSignPairing(kTRUE),
  fOptionLoopOverTracks(kTRUE),
  fOptionRunPrefilter(kTRUE),
  fOptionStoreJpsiCandidates(kFALSE),
  fEventCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts(),
  fPosTracks(),
  fNegTracks(),
  fPrefilterPosTracks(),
  fPrefilterNegTracks(),
  fJpsiCandidates(),
  fLegCandidatesMCcuts(),
  fJpsiMotherMCcuts(),
  fJpsiElectronMCcuts()
{
  //
  // named constructor
  //
   fEventCuts.SetOwner(kTRUE);
   fTrackCuts.SetOwner(kTRUE);
   fPreFilterTrackCuts.SetOwner(kTRUE);
   fPairCuts.SetOwner(kTRUE);
   fPreFilterPairCuts.SetOwner(kTRUE);
   fPosTracks.SetOwner(kFALSE);
   fNegTracks.SetOwner(kFALSE);
   fPrefilterPosTracks.SetOwner(kFALSE);
   fPrefilterNegTracks.SetOwner(kFALSE);
   fJpsiCandidates.SetOwner(kTRUE);
   fLegCandidatesMCcuts.SetOwner(kTRUE);
   fJpsiMotherMCcuts.SetOwner(kTRUE);
   fJpsiElectronMCcuts.SetOwner(kTRUE);
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2eeMult::~AliReducedAnalysisJpsi2eeMult() 
{
  //
  // destructor
  //
   fEventCuts.Clear("C"); fTrackCuts.Clear("C"); fPreFilterTrackCuts.Clear("C"); fPreFilterPairCuts.Clear("C"); fPairCuts.Clear("C");
   fPosTracks.Clear("C"); fNegTracks.Clear("C"); fPrefilterPosTracks.Clear("C"); fPrefilterNegTracks.Clear("C");
   fJpsiCandidates.Clear("C");
   if(fHistosManager) delete fHistosManager;
   if(fMixingHandler) delete fMixingHandler;
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::AddTrackCut(AliReducedInfoCut* cut) {
   //
   // Add a new cut
   //
   fTrackCuts.Add(cut); 
   fMixingHandler->SetNParallelCuts(fMixingHandler->GetNParallelCuts()+1);
   TString histClassNames = fMixingHandler->GetHistClassNames();
   histClassNames += Form("PairMEPP_%s;", cut->GetName());
   histClassNames += Form("PairMEPM_%s;", cut->GetName());
   histClassNames += Form("PairMEMM_%s;", cut->GetName());
   fMixingHandler->SetHistClassNames(histClassNames.Data());
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if(fEventCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
    if(values) { if(!cut->IsSelected(event, values)) return kFALSE; }
    else { if(!cut->IsSelected(event)) return kFALSE; }
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsTrackSelected(AliReducedTrackInfo* track, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if(fTrackCuts.GetEntries()==0) return kTRUE;
  track->ResetFlags();
  
  for(Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
    if(values) { if(cut->IsSelected(track, values)) track->SetFlag(i); }
    else { if(cut->IsSelected(track)) track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsTrackPrefilterSelected(AliReducedTrackInfo* track, Float_t* values/*=0x0*/) {
   //
   // apply track prefilter cuts
   //
   if(fPreFilterTrackCuts.GetEntries()==0) return kTRUE;
   
   for(Int_t i=0; i<fPreFilterTrackCuts.GetEntries(); ++i) {
      // if there are more cuts specified, we apply an AND on all of them
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fPreFilterTrackCuts.At(i);
      if(values) { if(!cut->IsSelected(track, values)) return kFALSE; }
      else { if(!cut->IsSelected(track)) return kFALSE; }
   }
   return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsPairSelected(Float_t* values) {
  //
  // apply pair cuts
  //
  if(fPairCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
    if(!cut->IsSelected(values)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsPairPreFilterSelected(Float_t* values) {
   //
   // apply event cuts
   //
   if(fPreFilterPairCuts.GetEntries()==0) return kTRUE;
   // loop over all the cuts and make a logical OR between all cuts in the list
   for(Int_t i=0; i<fPreFilterPairCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fPreFilterPairCuts.At(i);
      if(cut->IsSelected(values)) return kTRUE;
   }
   return kFALSE;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::Init() {
  //
  // initialize stuff
  //
   AliReducedVarManager::SetDefaultVarNames();
   // make sure variables needed to create jpsi candidate objects are filled
   if(fOptionStoreJpsiCandidates) {
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPt);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPhi);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kEta);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kMass);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPairType);
   }
   fHistosManager->SetUseDefaultVariableNames(kTRUE);
   fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
   
   fMixingHandler->SetHistogramManager(fHistosManager);
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::Process() {
  //
  // process the current event
  //  
  if(!fEvent) return;
//   AliReducedEventInfo* eventInfo = NULL;
//   if(fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;
//   else {
//      cout << "ERROR: AliReducedAnalysisJpsi2eeMult::Process() needs AliReducedEventInfo events" << endl;
//      return;
//   }
  if(fOptionRunOverMC) {
     if(fEventCounter%10000==0) 
        cout << "Event no. " << fEventCounter << endl;
  }
  else {
    if(fEventCounter%100000==0) 
       cout << "Event no. " << fEventCounter << endl;
  }
  fEventCounter++;
  
  AliReducedVarManager::SetEvent(fEvent);
  
  // reset the values array, keep only the run wise data (LHC and ALICE GRP information)
  // NOTE: the run wise data will be updated automatically in the VarManager in case a run change is detected
  for(Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i) fValues[i]=-9999.;
  
  // fill event information before event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  fHistosManager->FillHistClass("Event_BeforeCuts", fValues);
  for(UShort_t ibit=0; ibit<64; ++ibit) {
     AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
     fHistosManager->FillHistClass("EventTag_BeforeCuts", fValues);
  }
  for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_BeforeCuts", fValues);
  }
  
  LoopOverMCTracks(1, false);
  
  // apply event selection
  if(!IsEventSelected(fEvent , fValues )) return;
  if(fOptionRunOverMC) FillMCTruthHistograms();
  
  // select tracks
  if(fOptionLoopOverTracks)
    RunTrackSelection();
    
  // Run the prefilter  
  // NOTE: Pair each track from the selected tracks list with all selected tracks in the prefilter track list
  //         If the created pair fails the pair prefilter criteria, then the selected trak is removed from the track list
  //          and further pairing
  //FillTrackHistograms("Track_BeforePrefilter");
  //RunSameEventPairing("PairPrefilterSE");
  if(fOptionLoopOverTracks && fOptionRunPrefilter)
    RunPrefilter();
  
  if(fOptionLoopOverTracks) {
    fValues[AliReducedVarManager::kNtracksPosAnalyzed] = fPosTracks.GetEntries();
    fValues[AliReducedVarManager::kNtracksNegAnalyzed] = fNegTracks.GetEntries();
    fValues[AliReducedVarManager::kNtracksAnalyzed] = fValues[AliReducedVarManager::kNtracksNegAnalyzed]+fValues[AliReducedVarManager::kNtracksPosAnalyzed];
    fValues[AliReducedVarManager::kEvAverageTPCchi2] /= (fPosTracks.GetEntries()+fNegTracks.GetEntries()>0 ? fValues[AliReducedVarManager::kNtracksAnalyzed] : 1.0);
    fValues[AliReducedVarManager::kEvAveragegoldenchi2] /= (fPosTracks.GetEntries()+fNegTracks.GetEntries()>0 ? fValues[AliReducedVarManager::kNtracksAnalyzed] : 1.0);
    fValues[AliReducedVarManager::kEvAverageITSchi2] /= (fPosTracks.GetEntries()+fNegTracks.GetEntries()>0 ? fValues[AliReducedVarManager::kNtracksAnalyzed] : 1.0); 
  }
  
  // Fill track histograms
  if(fOptionLoopOverTracks)
    FillTrackHistograms();
  
  // Feed the selected tracks to the event mixing handler 
  if(!fOptionRunOverMC && fOptionRunMixing && fOptionRunPairing)
    fMixingHandler->FillEvent(&fPosTracks, &fNegTracks, fValues, AliReducedPairInfo::kJpsiToEE);
  
  // Do the same event pairing
  if(fOptionRunPairing)
    RunSameEventPairing();
 
  // fill event info histograms after cuts
  fHistosManager->FillHistClass("Event_AfterCuts", fValues);
  for(UShort_t ibit=0; ibit<64; ++ibit) {
     AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
     fHistosManager->FillHistClass("EventTag_AfterCuts", fValues);
  }
  for(UShort_t ibit=0; ibit<64; ++ibit) {
     AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
     fHistosManager->FillHistClass("EventTriggers_AfterCuts", fValues);
  }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillTrackHistograms(TString trackClass /*= "Track"*/) {
   //
   // Fill all track histograms
   //
   for(Int_t i=0;i<36; ++i) fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+i] = 0.;
   AliReducedTrackInfo* track=0;
   TIter nextPosTrack(&fPosTracks);
   for(Int_t i=0;i<fPosTracks.GetEntries();++i) {
      track = (AliReducedTrackInfo*)nextPosTrack();
      //Int_t tpcSector = TMath::FloorNint(18.*track->Phi()/TMath::TwoPi());
      fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
      // reset track variables
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kEMCALmatchedEOverP; ++i) fValues[i]=-9999.;
      
      AliReducedVarManager::FillTrackInfo(track, fValues);
      FillTrackHistograms(track, Form("%s+", trackClass.Data()) );
      FillTrackHistograms(track, trackClass);
   }
   TIter nextNegTrack(&fNegTracks);
   for(Int_t i=0;i<fNegTracks.GetEntries();++i) {
      track = (AliReducedTrackInfo*)nextNegTrack();
      //Int_t tpcSector = TMath::FloorNint(18.*track->Phi()/TMath::TwoPi());
      fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
      // reset track variables
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kEMCALmatchedEOverP; ++i) fValues[i]=-9999.;
      
      AliReducedVarManager::FillTrackInfo(track, fValues);
      FillTrackHistograms(track, Form("%s-", trackClass.Data()) );
      FillTrackHistograms(track, trackClass);
   }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillTrackHistograms(AliReducedTrackInfo* track, TString trackClass /*="Track"*/) {
   //
   // fill track level histograms
   //
   UInt_t mcDecisionMap = 0;
   if(fOptionRunOverMC) mcDecisionMap = CheckReconstructedLegMCTruth(track);      
   
   for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
      if(track->TestFlag(icut)) {
         fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
         if(mcDecisionMap) {       // Fill histograms for tracks identified as MC truth
            for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
               if(mcDecisionMap & (UInt_t(1)<<iMC))
                  fHistosManager->FillHistClass(Form("%s_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(), 
                                                                  fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
            }
         }
         
         if(track->IsA() != AliReducedTrackInfo::Class()) continue;
         
         AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(track);
         if(!trackInfo) continue;
         
         for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
            AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
            fHistosManager->FillHistClass(Form("%sStatusFlags_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sStatusFlags_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                                     fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
         for(Int_t iLayer=0; iLayer<6; ++iLayer) {
            AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
            fHistosManager->FillHistClass(Form("%sITSclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sITSclusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                                     fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
         for(Int_t iLayer=0; iLayer<8; ++iLayer) {
            AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
            fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                                     fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
      } // end if(track->TestFlag(icut))
   }  // end loop over cuts
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillPairHistograms(ULong_t mask, Int_t pairType, TString pairClass /*="PairSE"*/, UInt_t mcDecisions /* = 0*/) {
   //
   // fill pair level histograms
   // NOTE: pairType can be 0,1 or 2 corresponding to ++, +- or -- pairs
   TString typeStr[3] = {"PP", "PM", "MM"};
   for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
      if(mask & (ULong_t(1)<<icut)) {
         fHistosManager->FillHistClass(Form("%s%s_%s", pairClass.Data(), typeStr[pairType].Data(), fTrackCuts.At(icut)->GetName()), fValues, icut);
         if(mcDecisions && pairType==1) {
            for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
               if(mcDecisions & (UInt_t(1)<<iMC))
                  fHistosManager->FillHistClass(Form("%s%s_%s_%s", pairClass.Data(), typeStr[pairType].Data(), fTrackCuts.At(icut)->GetName(), fLegCandidatesMCcuts.At(iMC)->GetName()), fValues, icut);
            }
         }
      }
   }  // end loop over cuts
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::RunTrackSelection() {
   //
   // select electron candidates and prefilter tracks
   //
   // clear the track arrays
   fPosTracks.Clear("C"); fNegTracks.Clear("C"); fPrefilterPosTracks.Clear("C"); fPrefilterNegTracks.Clear("C");
   fValues[AliReducedVarManager::kEvAverageTPCchi2] = 0.0;
   fValues[AliReducedVarManager::kEvAverageITSchi2] = 0.0;
   fValues[AliReducedVarManager::kEvAveragegoldenchi2] = 0.0;
   fValues[AliReducedVarManager::kEvAverageITSchi2Nonzero] = 0.0;
   
   // loop over the track list(s) and evaluate all the track cuts
   LoopOverTracks(1);      // first array
   LoopOverTracks(2);      // second array (if used)
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::LoopOverTracks(Int_t arrayOption /*=1*/) {
   //
   // Loop over a given track array, apply cuts and add selected tracks to arrays
   //
   AliReducedTrackInfo* track = 0x0;
   TClonesArray* trackList = (arrayOption==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if (!trackList) return;

   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      // do not loop over pure MC truth tracks 
      // NOTE: this can be also handled via AliReducedTrackCut::SetRejectPureMC()
      if(fOptionRunOverMC && track->IsMCTruth()) continue;     
      // reset track variables
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kEMCALmatchedEOverP; ++i) fValues[i]=-9999.;

      AliReducedVarManager::FillTrackInfo(track, fValues);
      fHistosManager->FillHistClass("Track_BeforeCuts", fValues);
      
      if(track->IsA() == AliReducedTrackInfo::Class()) {
         AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(track);
         if(trackInfo) {
            for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
               AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
               fHistosManager->FillHistClass("TrackStatusFlags_BeforeCuts", fValues);
            }
            for(Int_t iLayer=0; iLayer<6; ++iLayer) {
               AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
               fHistosManager->FillHistClass("TrackITSclusterMap_BeforeCuts", fValues);
            }
            for(Int_t iLayer=0; iLayer<8; ++iLayer) {
               AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
               fHistosManager->FillHistClass("TrackTPCclusterMap_BeforeCuts", fValues);
            }
         }
      }
      
      if(IsTrackSelected(track, fValues)) {
         if(track->Charge()>0) fPosTracks.Add(track);
         if(track->Charge()<0) fNegTracks.Add(track);
         
         if(track->IsA() == AliReducedTrackInfo::Class()){
            
         }
      }
      if(IsTrackPrefilterSelected(track, fValues)) {
         if(track->Charge()>0) fPrefilterPosTracks.Add(track);
         if(track->Charge()<0) fPrefilterNegTracks.Add(track);
   //         fValues[AliReducedVarManager::kEvAverageTPCchi2] += ((AliReducedTrackInfo*)track)->TPCchi2();
     //       fValues[AliReducedVarManager::kEvAverageITSchi2] += ((AliReducedTrackInfo*)track)->ITSchi2();
       //     fValues[AliReducedVarManager::kEvAveragegoldenchi2] += ((AliReducedTrackInfo*)track)->Chi2TPCConstrainedVsGlobal();
         //   fValues[AliReducedVarManager::kEvAverageITSchi2Nonzero] = 1; 
      }
   }   // end loop over tracks
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::RunSameEventPairing(TString pairClass /*="PairSE"*/) {
   //
   // Run the same event pairing
   //
   if(fOptionStoreJpsiCandidates) fJpsiCandidates.Clear("C");
   fValues[AliReducedVarManager::kNpairsSelected] = 0;
   
   TIter nextPosTrack(&fPosTracks);
   TIter nextNegTrack(&fNegTracks);
   
   AliReducedTrackInfo* pTrack=0;
   AliReducedTrackInfo* pTrack2=0;
   AliReducedTrackInfo* nTrack=0;
   AliReducedTrackInfo* nTrack2=0;
   for(Int_t ip=0; ip<fPosTracks.GetEntries(); ++ip) {
      pTrack = (AliReducedTrackInfo*)nextPosTrack();
      
      nextNegTrack.Reset();
      for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
         nTrack = (AliReducedTrackInfo*)nextNegTrack();
         
         // verify that the two current tracks have at least 1 common bit
         if(!(pTrack->GetFlags() & nTrack->GetFlags())) continue;
         if(fOptionRunRotation){
            AliReducedTrackInfo pTrackCopy((*pTrack));
            AliReducedTrackInfo nTrackCopy((*nTrack));
            for(int i=0; i<fNRotations; ++i){
               RunTrackRotation(pTrackCopy, nTrackCopy, 1);
             }
         }
         AliReducedVarManager::FillPairInfo(pTrack, nTrack, AliReducedPairInfo::kJpsiToEE, fValues);
         if(IsPairSelected(fValues)) {
            FillPairHistograms(pTrack->GetFlags() & nTrack->GetFlags(), 1, pairClass, (fOptionRunOverMC ? CheckReconstructedLegMCTruth(pTrack, nTrack) : 0));    // 1 is for +- pairs 
            fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
            if(fOptionStoreJpsiCandidates) {
               AliReducedPairInfo* pair = new AliReducedPairInfo();
               pair->SetFlags(pTrack->GetFlags() & nTrack->GetFlags());
               pair->PtPhiEta(fValues[AliReducedVarManager::kPt], fValues[AliReducedVarManager::kPhi], fValues[AliReducedVarManager::kEta]);
               pair->SetMass(fValues[AliReducedVarManager::kMass]);
               pair->CandidateId(AliReducedPairInfo::kJpsiToEE);
               pair->PairType(1);
               pair->SetLegIds(pTrack->TrackId(), nTrack->TrackId());
               fJpsiCandidates.Add(pair);
            }
         }
      }  // end loop over negative tracks
      
      if(fOptionRunLikeSignPairing) {
         for(Int_t ip2=ip+1; ip2<fPosTracks.GetEntries(); ++ip2) {
            pTrack2 = (AliReducedTrackInfo*)fPosTracks.At(ip2);
         
            // verify that the two current tracks have at least 1 common bit
            if(!(pTrack->GetFlags() & pTrack2->GetFlags())) continue;
            if(fOptionRunRotation){
              AliReducedTrackInfo pTrackCopy((*pTrack));
              AliReducedTrackInfo pTrack2Copy((*pTrack2));
              for(int i=0; i<fNRotations; ++i){
                RunTrackRotation(pTrackCopy, pTrack2Copy, 0);
              }
            }
            AliReducedVarManager::FillPairInfo(pTrack, pTrack2, AliReducedPairInfo::kJpsiToEE, fValues);
            if(IsPairSelected(fValues)) {
               FillPairHistograms(pTrack->GetFlags() & pTrack2->GetFlags(), 0, pairClass);       // 0 is for ++ pairs 
               fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
               if(fOptionStoreJpsiCandidates) {
                  AliReducedPairInfo* pair = new AliReducedPairInfo();
                  pair->SetFlags(pTrack->GetFlags() & pTrack2->GetFlags());
                  pair->PtPhiEta(fValues[AliReducedVarManager::kPt], fValues[AliReducedVarManager::kPhi], fValues[AliReducedVarManager::kEta]);
                  pair->SetMass(fValues[AliReducedVarManager::kMass]);
                  pair->CandidateId(AliReducedPairInfo::kJpsiToEE);
                  pair->PairType(0);
                  pair->SetLegIds(pTrack->TrackId(), pTrack2->TrackId());
                  fJpsiCandidates.Add(pair);
               }
            }
         }  // end loop over positive tracks
      }
   }  // end loop over positive tracks
   
   if(fOptionRunLikeSignPairing) {
      nextNegTrack.Reset();
      for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
         nTrack = (AliReducedTrackInfo*)nextNegTrack();
      
         for(Int_t in2=in+1; in2<fNegTracks.GetEntries(); ++in2) {
            nTrack2 = (AliReducedTrackInfo*)fNegTracks.At(in2);
         
            // verify that the two current tracks have at least 1 common bit
            if(!(nTrack->GetFlags() & nTrack2->GetFlags())) continue;
            if(fOptionRunRotation){
              AliReducedTrackInfo nTrack2Copy((*nTrack2));
              AliReducedTrackInfo nTrackCopy((*nTrack));
              for(int i=0; i<fNRotations; ++i){
                RunTrackRotation(nTrack2Copy, nTrackCopy, 2);
              }
            }
            AliReducedVarManager::FillPairInfo(nTrack, nTrack2, AliReducedPairInfo::kJpsiToEE, fValues);
            if(IsPairSelected(fValues)) {
               FillPairHistograms(nTrack->GetFlags() & nTrack2->GetFlags(), 2, pairClass);      // 2 is for -- pairs
               fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
               if(fOptionStoreJpsiCandidates) {
                  AliReducedPairInfo* pair = new AliReducedPairInfo();
                  pair->SetFlags(nTrack->GetFlags() & nTrack2->GetFlags());
                  pair->PtPhiEta(fValues[AliReducedVarManager::kPt], fValues[AliReducedVarManager::kPhi], fValues[AliReducedVarManager::kEta]);
                  pair->SetMass(fValues[AliReducedVarManager::kMass]);
                  pair->CandidateId(AliReducedPairInfo::kJpsiToEE);
                  pair->PairType(2);
                  pair->SetLegIds(nTrack->TrackId(), nTrack2->TrackId());
                  fJpsiCandidates.Add(pair);
               }
            }
         }  // end loop over negative tracks
      }  // end loop over negative tracks
   }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::RunPrefilter() {
   //
   // Run the prefilter selection
   // At this point it is assumed that the track lists are filled
   //
   TIter nextPosTrack(&fPosTracks);
   TIter nextNegTrack(&fNegTracks);
   TIter nextPosPrefilterTrack(&fPrefilterPosTracks);
   TIter nextNegPrefilterTrack(&fPrefilterNegTracks);
   
   // First pair the positive trackes with the prefilter selected tracks
   AliReducedTrackInfo* track=0;
   AliReducedTrackInfo* trackPref=0;
   for(Int_t ip = 0; ip<fPosTracks.GetEntries(); ++ip) {
      track = (AliReducedTrackInfo*)nextPosTrack();
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedTrackInfo*)nextPosPrefilterTrack();
         
         if(track->TrackId()==trackPref->TrackId()) continue;       // avoid self-pairing
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over positive prefilter tracks
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedTrackInfo*)nextNegPrefilterTrack();
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over negative prefilter tracks
   }  // end loop over the positive tracks

   for(Int_t in = 0; in<fNegTracks.GetEntries(); ++in) {
      track = (AliReducedTrackInfo*)nextNegTrack();
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedTrackInfo*)nextPosPrefilterTrack();
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over positive prefilter tracks
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedTrackInfo*)nextNegPrefilterTrack();
         
         if(track->TrackId()==trackPref->TrackId()) continue;       // avoid self-pairing
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over negative prefilter tracks
   }  // end loop over the negative tracks

   // remove tracks
   nextPosTrack.Reset();
   for(Int_t ip = fPosTracks.GetEntries()-1 ; ip >= 0; --ip) {
     track = (AliReducedTrackInfo*)nextPosTrack();
     if(!track->GetFlags()) {
        fPosTracks.Remove(track);
     }
   }
  nextNegTrack.Reset();
  for(Int_t ip = fNegTracks.GetEntries()-1 ; ip >= 0; --ip) {
    track = (AliReducedTrackInfo*)nextNegTrack();
    if(!track->GetFlags()) {
       fNegTracks.Remove(track);
    }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::Finish() {
  //
  // run stuff after the event loop
  //
   if(fOptionRunMixing && !fOptionRunOverMC)
     fMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
}


//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2eeMult::CheckReconstructedLegMCTruth(AliReducedTrackInfo* track) {
   //
   // Check a reconstructed track against all the specified MC truth cuts
   //
   // TODO: In the fLegCandidatesMCcuts one can also add AliSignalMC objects which can then be tested
   //             using the AliReducedTrackInfo::fMCPdg[]
   //
   if(fLegCandidatesMCcuts.GetEntries()==0) return 0;
   
   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fLegCandidatesMCcuts.At(i);
      if(cut->IsSelected(track)) 
         decisionMap |= (UInt_t(1)<<i);
   }
   
   return decisionMap;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2eeMult::CheckReconstructedLegMCTruth(AliReducedTrackInfo* ptrack, AliReducedTrackInfo* ntrack) {
   //
   // check the pair of tracks to see if they match the defined MC cuts and in addition
   // that they have the same mother
   // NOTE: The condition for the 2 tracks to have the same mother requires information on the MC label,
   //             which is available just in the full track information (AliReducedTrackInfo::fMCLabels[]).
   //           The consequence is that for the jpsi2ee analysis, the reconstructed tracks need to be always written as full tracks 
   //
   
   // check that both tracks are full tracks
   if(ptrack->IsA() != AliReducedTrackInfo::Class()) return 0;
   if(ntrack->IsA() != AliReducedTrackInfo::Class()) return 0;
   
   // check that the tracks have the same mother
   if(TMath::Abs(((AliReducedTrackInfo*)ptrack)->MCLabel(1)) != TMath::Abs(((AliReducedTrackInfo*)ntrack)->MCLabel(1))) return 0;
   
   // check the MC requirements on each of the leg and their logical intersection
   if(fLegCandidatesMCcuts.GetEntries()==0) return 0;
   UInt_t pTrackDecisions = CheckReconstructedLegMCTruth(ptrack);
   if(!pTrackDecisions) return 0;
   UInt_t nTrackDecisions = CheckReconstructedLegMCTruth(ntrack);
   return (pTrackDecisions & nTrackDecisions);
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillMCTruthHistograms() {
  //
  // fill histograms with pure signal
  //   
   // loop over the first track array
  LoopOverMCTracks(1);
  // and over the second
  // NOTE: In the current model, handling the MC truth info requires the labels, which are properties of the full track,
  //         so there is no point in looping over the second track array which, if it exists, contains just base tracks
  //LoopOverMCTracks(2);
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::LoopOverMCTracks(Int_t trackArray /*=1*/, Bool_t fillHistos /*= kTRUE*/) {
   //
   // loop over the track array and check the pure MC tracks against the defined MC selections
   //   
   AliReducedTrackInfo* mother=0x0;
   AliReducedTrackInfo* daughter1 = 0x0;
   AliReducedTrackInfo* daughter2 = 0x0;
   
   TClonesArray* trackList = (trackArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if(!trackList) return;
   TIter nextTrack(trackList);
   
   fValues[AliReducedVarManager::kNMCtruthJpsi  ] = 0;
   fValues[AliReducedVarManager::kNMCtruthJpsiLegs ] = 0;
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      mother = (AliReducedTrackInfo*)nextTrack();
      if(!mother->IsMCKineParticle()) continue;
      
      // apply selections on the jpsi mother
      UInt_t motherDecisions = CheckMotherMCTruth(mother);
      if(!motherDecisions) continue;
      fValues[AliReducedVarManager::kNMCtruthJpsi  ] ++;
      // find the jpsi daughters (needed to compute 2-track properties like the polarization, etc.)
      Int_t daughter1Label = 0; Int_t daughter2Label = 0;
      FindJpsiTruthLegs(mother, daughter1Label, daughter2Label);   
      daughter1 = FindMCtruthTrackByLabel(daughter1Label);
      daughter2 = FindMCtruthTrackByLabel(daughter2Label);
      
      if(fillHistos){
        // reset track variables and fill info
        for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kEMCALmatchedEOverP; ++i) fValues[i]=-9999.;
        AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);
        
        // loop over jpsi mother selections and fill histograms before the kine cuts on electrons
        for(Int_t iCut = 0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
          if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
          fHistosManager->FillHistClass(Form("%s_PureMCTruth_BeforeSelection", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);         
        }
      }
      if(!daughter1) continue;
      if(!daughter2) continue;
      
      // apply selections on pure MC daughter electrons (kine cuts)
      UInt_t daughter1Decisions = CheckDaughterMCTruth(daughter1);
      if(!daughter1Decisions) continue;
      UInt_t daughtersDecisions = daughter1Decisions & CheckDaughterMCTruth(daughter2);
      if(!daughtersDecisions) continue;
      
      fValues[AliReducedVarManager::kNMCtruthJpsiLegs  ] ++;
      if(fillHistos){
        for(Int_t iCut = 0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
          if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
          if(!(daughtersDecisions & (UInt_t(1)<<iCut)))  continue;
          fHistosManager->FillHistClass(Form("%s_PureMCTruth_AfterSelection", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);         
        }
      }
   }  // end loop over tracks
   return;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2eeMult::CheckMotherMCTruth(AliReducedTrackInfo* mother) {
   //
   // Check the mother pure MC truth against all defined selections and return a bit map with all decisions
   //
   if(fJpsiMotherMCcuts.GetEntries()==0) return 0;
   
   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fJpsiMotherMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiMotherMCcuts.At(i);
      if(cut->IsSelected(mother)) 
         decisionMap |= (UInt_t(1)<<i);
   }
   
   return decisionMap;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2eeMult::CheckDaughterMCTruth(AliReducedTrackInfo* daughter) {
   //
   // Check the daughter pure MC truth against all defined selections and return a bit map with all decisions
   //
   if(fJpsiElectronMCcuts.GetEntries()==0) return 0;
   
   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fJpsiElectronMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiElectronMCcuts.At(i);
      if(cut->IsSelected(daughter)) 
         decisionMap |= (UInt_t(1)<<i);
   }
   
   return decisionMap;
}

//___________________________________________________________________________
AliReducedTrackInfo* AliReducedAnalysisJpsi2eeMult::FindMCtruthTrackByLabel(Int_t label) {
   //
   // search the track list for pure MC track with label and return the track pointer
   //
   AliReducedTrackInfo* track=0x0;
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCKineParticle()) continue;
      if(track->MCLabel(0)==label) return track;
   }
   return 0x0;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label) {
   //
   // find the jpsi legs in the list of pure MC truth particles
   //
   Int_t mLabel = mother->MCLabel(0);
   AliReducedTrackInfo* track=0x0;
   
   Int_t legsFound = 0;
   
   // loop over the first track array
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      if(legsFound==2) return;
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCKineParticle()) continue;
      if(track->MCLabel(1)==mLabel && TMath::Abs(track->MCPdg(0))==11) {
         legsFound += 1;
         if(legsFound==1) leg1Label = track->MCLabel(0);
         if(legsFound==2) leg2Label = track->MCLabel(0);
      }
   }
   return;
}

void AliReducedAnalysisJpsi2eeMult::RunTrackRotation(AliReducedTrackInfo &pTrack, AliReducedTrackInfo &nTrack, Int_t pairType){

  TString pairClass = "PairTR";
  Double_t phi1 = TMath::TwoPi() * gRandom->Rndm();
  Double_t phi2 = TMath::TwoPi() * gRandom->Rndm();

  if(pTrack.IsCartesian()){
    pTrack.Px( pTrack.Pt() * TMath::Cos(phi1)  );
    pTrack.Py( pTrack.Pt() * TMath::Sin(phi1)  );
    nTrack.Px( nTrack.Pt() * TMath::Cos(phi2)  );
    nTrack.Py( nTrack.Pt() * TMath::Sin(phi2)  );
  }
  else{
    pTrack.Phi( phi1 );
    nTrack.Phi( phi2 );
  }

  AliReducedVarManager::FillPairInfo( (&pTrack), (&nTrack), AliReducedPairInfo::kJpsiToEE, fValues);
  if(IsPairSelected(fValues)) {
    FillPairHistograms(pTrack.GetFlags() & nTrack.GetFlags(), pairType, pairClass );
  }
}


void AliReducedAnalysisJpsi2eeMult::SetRunTrackRotation( Bool_t option){
  fOptionRunRotation = option;

  if(fOptionRunRotation) gRandom->SetSeed();

}
