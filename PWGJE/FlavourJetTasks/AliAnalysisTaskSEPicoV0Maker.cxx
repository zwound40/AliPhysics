/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for RecoDecay object (K0 short, Lambda,
// D mesons ...) filtering
//
// Author: X-M. Zhang, xmzhang@lbl.gov
/////////////////////////////////////////////////////////////

#include <iostream>

#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliVVertex.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliV0vertexer.h"
#include "AliAnalysisUtils.h"
#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliPIDResponse.h"

#include "AliPicoHeaderCJ.h"
#include "AliPicoV0RD.h"
#include "AliPicoV0MC.h"
#include "AliAnalysisTaskSEPicoV0Maker.h"

ClassImp(AliAnalysisTaskSEPicoV0Maker)
//=============================================================================

const Double_t AliAnalysisTaskSEPicoV0Maker::fgkMassPion   = 0.13957;
const Double_t AliAnalysisTaskSEPicoV0Maker::fgkMassKshort = 0.497614;
const Double_t AliAnalysisTaskSEPicoV0Maker::fgkMassProton = 0.938272;
const Double_t AliAnalysisTaskSEPicoV0Maker::fgkMassLambda = 1.11568;

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Maker::AliAnalysisTaskSEPicoV0Maker() :
AliAnalysisTaskSE(),
fEventAOD(0),
fEventESD(0),
fCentInfo(0),
fRespoPID(0),
fAnaUtils(0),
fEventAcptMask(0),
fTriggerMask(),
fCollisionType(0),
fCentEst(""),
fIsAnaInfoMC(kFALSE),
fIsRefitV0sESD(kFALSE),
fIsSkipFastOnly(kFALSE),
fIsDPMjetMC(kFALSE),
fRapidityShift(0.),
fCutMinEventVtxContr(0),
fCutMaxEventVzAbs(0.),
fCutMinV0Pt(0.),
fCutMaxV0Pt(0.),
fCutMinV0Rap(0.),
fCutMaxV0Rap(0.),
fCutMinDauPt(0.),
fCutMinDauEta(0.),
fCutMaxDauEta(0.),
fCutMaxV0Chi2(0.),
fCutMinV0Radius(0.),
fCutMaxV0Radius(0.),
fCutMaxDausDCA(0.),
fCutMinDauDCAtoPV(0.),
fCutMinDauXrowsTPC(0.),
fCutMinDauXrowsOverFindableClusTPC(0.),
fCutMaxKshortSigmaTPC(0.),
fCutMinKshortCosPA(0.),
fCutMaxKshortCtau(0.),
fCutMaxKshortArmFrac(0.),
fCutMinKshortDeltaM(0.),
fCutMaxLambdaSigmaTPC(0.),
fCutMinLambdaCosPA(0.),
fCutMaxLambdaCtau(0.),
fCutMaxLambdaArmFrac(0.),
fCutMinLambdaDeletaM(0.),
fPicoV0sClArr(0),
fOutputListEH(0),
fOutputListMC(0)
{
//
// Default constructor
//

  for (Int_t i=3; i--;) fPrimaryVtx[i] = 0.;
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Maker::AliAnalysisTaskSEPicoV0Maker(const char *name, Bool_t bIsMC) :
AliAnalysisTaskSE(name),
fEventAOD(0),
fEventESD(0),
fCentInfo(0),
fRespoPID(0),
fAnaUtils(0),
fEventAcptMask(0),
fTriggerMask(),
fCollisionType(AliPicoHeaderCJ::kPP),
fCentEst("V0M"),
fIsAnaInfoMC(bIsMC),
fIsRefitV0sESD(kFALSE),
fIsSkipFastOnly(kFALSE),
fIsDPMjetMC(kFALSE),
fRapidityShift(0.),
fCutMinEventVtxContr(2),
fCutMaxEventVzAbs(10.),
fCutMinV0Pt(0.),
fCutMaxV0Pt(100.),
fCutMinV0Rap(-10.),
fCutMaxV0Rap(10.),
fCutMinDauPt(0.),
fCutMinDauEta(-10.),
fCutMaxDauEta(10.),
fCutMaxV0Chi2(33.),
fCutMinV0Radius(0.3),
fCutMaxV0Radius(200.),
fCutMaxDausDCA(1.5),
fCutMinDauDCAtoPV(0.05),
fCutMinDauXrowsTPC(70.),
fCutMinDauXrowsOverFindableClusTPC(0.8),
fCutMaxKshortSigmaTPC(-1.),
fCutMinKshortCosPA(0.95),
fCutMaxKshortCtau(30.),
fCutMaxKshortArmFrac(-1.),
fCutMinKshortDeltaM(0.003),
fCutMaxLambdaSigmaTPC(7.),
fCutMinLambdaCosPA(0.993),
fCutMaxLambdaCtau(40.),
fCutMaxLambdaArmFrac(-1.),
fCutMinLambdaDeletaM(-1.),
fPicoV0sClArr(0),
fOutputListEH(0),
fOutputListMC(0)
{
//
// Constructor
//

  for (Int_t i=3; i--;) fPrimaryVtx[i] = 0.;

  DefineOutput(1, TList::Class());
  if (fIsAnaInfoMC) DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Maker::~AliAnalysisTaskSEPicoV0Maker()
{
//
// Default destructor
//


  if (fEventAOD) { delete fEventAOD; fEventAOD = 0; }
  if (fEventESD) { delete fEventESD; fEventESD = 0; }
  if (fCentInfo) { delete fCentInfo; fCentInfo = 0; }
  if (fRespoPID) { delete fRespoPID; fRespoPID = 0; }
  if (fAnaUtils) { delete fAnaUtils; fAnaUtils = 0; }

  if (fPicoV0sClArr) { delete fPicoV0sClArr; fPicoV0sClArr = 0; }
  if (fOutputListEH) { delete fOutputListEH; fOutputListEH = 0; }
  if (fOutputListMC) { delete fOutputListMC; fOutputListMC = 0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::Init()
{
//
//  AliAnalysisTaskSEPicoV0Maker::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskSEPicoV0Maker::UserCreateOutputObjects
//

  InitAnalysis();
//=============================================================================

  if (fIsAnaInfoMC) {
    fPicoV0sClArr = new TClonesArray("AliPicoV0MC");
    fPicoV0sClArr->SetName("PicoV0s");
  } else {
    fPicoV0sClArr = new TClonesArray("AliPicoV0RD");
    fPicoV0sClArr->SetName("PicoV0s");
  }
//=============================================================================

  fOutputListEH = new TList();
  fOutputListEH->SetOwner();
  CreateHistogramsEH();
  PostData(1, fOutputListEH);

  if (fIsAnaInfoMC) {
    fOutputListMC = new TList();
    fOutputListMC->SetOwner();
    CreateHistogramsMC();
    PostData(2, fOutputListMC);
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::UserExec(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0Maker::UserExec
//

  fPicoV0sClArr->Delete();
  if (IsEventNotAcpt()) return;
  if (!(InputEvent()->FindListObject("PicoV0s"))) InputEvent()->AddObject(fPicoV0sClArr);
//=============================================================================

  FillHistogramsEH();
  if (IsEventNotINEL()) return;
  if (fIsAnaInfoMC) FillHistogramsMC();
//=============================================================================

  if (IsEventNotMBsa()) return;


  FillPicoV0s();
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::Terminate(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0Maker::Terminate
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::NotifyRun()
{
//
//  AliAnalysisTaskSEPicoV0Maker::NotifyRun
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::FillPicoV0s()
{
//
//  AliAnalysisTaskSEPicoV0Maker::FillPicoRecoV0s
//

  Int_t nV0s = 0;
  AliAODv0 *pV0AOD = 0;
  AliESDv0 *pV0ESD = 0;
  if (fEventAOD) nV0s = fEventAOD->GetNumberOfV0s();
  if (fEventESD) nV0s = fEventESD->GetNumberOfV0s();
  if (nV0s<=0) return;
//=============================================================================


  TH2D *hKshortPtInvM = dynamic_cast<TH2D*>(fOutputListEH->FindObject("hKshortPtInvM"));
  TH2D *hLambdaPtInvM = dynamic_cast<TH2D*>(fOutputListEH->FindObject("hLambdaPtInvM"));
  TH2D *hAntiLaPtInvM = dynamic_cast<TH2D*>(fOutputListEH->FindObject("hAntiLaPtInvM"));
//=============================================================================

  AliPicoV0RD *pV0RD = 0;
  AliPicoV0MC *pV0MC = 0;
  Int_t nAt = fPicoV0sClArr->GetEntriesFast();

  for (Int_t iV0=0; iV0<nV0s; iV0++) {

    if (fEventAOD) {
      pV0AOD = fEventAOD->GetV0(iV0); if (!pV0AOD) continue;

      if (fIsAnaInfoMC) {
        pV0MC = SelectV0CandidateMC(pV0AOD); if (!pV0MC) { pV0AOD = 0; continue; }
      } else {
        pV0RD = SelectV0CandidateRD(pV0AOD); if (!pV0RD) { pV0AOD = 0; continue; }
      }
    }

    if (fEventESD) {
      pV0ESD = fEventESD->GetV0(iV0); if (!pV0ESD) continue;

      if (fIsAnaInfoMC) {
        pV0MC = SelectV0CandidateMC(pV0ESD); if (!pV0MC) { pV0ESD = 0; continue; }
      } else {
        pV0RD = SelectV0CandidateRD(pV0ESD); if (!pV0RD) { pV0ESD = 0; continue; }
      }
    }

    if (pV0RD) {
      pV0RD->FillKshortPtInvM(hKshortPtInvM);
      pV0RD->FillLambdaPtInvM(hLambdaPtInvM);
      pV0RD->FillAntiLaPtInvM(hAntiLaPtInvM);
      new ((*fPicoV0sClArr)[nAt++]) AliPicoV0RD(*pV0RD);
      delete pV0RD; pV0RD=0;
    }

    if (pV0MC) {
      pV0MC->FillKshortPtInvM(hKshortPtInvM);
      pV0MC->FillLambdaPtInvM(hLambdaPtInvM);
      pV0MC->FillAntiLaPtInvM(hAntiLaPtInvM);
      new ((*fPicoV0sClArr)[nAt++]) AliPicoV0MC(*pV0MC);
      delete pV0MC; pV0MC=0;
    }

    if (pV0AOD) pV0AOD=0;
    if (pV0ESD) pV0ESD=0;
  }

  return;
}

//_____________________________________________________________________________
AliPicoV0RD* AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD(AliAODv0 const *pV0)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD
//

  if  (pV0->GetOnFlyStatus()) return 0x0;
  if ((pV0->Chi2V0())>fCutMaxV0Chi2) return 0x0;
//=============================================================================

  Double_t dV0Pt  = pV0->Pt(); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return 0x0;
  Double_t dKaRap = pV0->RapK0Short(); if ((dKaRap<fCutMinV0Rap) || (dKaRap>fCutMaxV0Rap)) return 0x0;
  Double_t dLaRap = pV0->RapLambda();  if ((dLaRap<fCutMinV0Rap) || (dLaRap>fCutMaxV0Rap)) return 0x0;
//=============================================================================

  Double_t dV0Vtx[3]; pV0->GetXYZ(dV0Vtx);
  Double_t dV0Radius = TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]);
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return 0x0;

  Double_t dDausDCA = pV0->DcaV0Daughters(); if (dDausDCA>fCutMaxDausDCA) return 0x0;
  Double_t dPosDCAtoPV = pV0->DcaPosToPrimVertex(); if (dPosDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
  Double_t dNegDCAtoPV = pV0->DcaNegToPrimVertex(); if (dNegDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
//=============================================================================

  AliAODTrack *pDauPos = (AliAODTrack*)pV0->GetDaughter(0); if (!pDauPos) return 0x0;
  AliAODTrack *pDauNeg = (AliAODTrack*)pV0->GetDaughter(1); if (!pDauNeg) return 0x0;

  if (!(pDauPos->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;
  if (!(pDauNeg->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;

  if ((pDauPos->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return 0x0;
  if ((pDauNeg->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return 0x0;

  Float_t dPosXrowsTPC = pDauPos->GetTPCClusterInfo(2,1);
  Float_t dNegXrowsTPC = pDauNeg->GetTPCClusterInfo(2,1);
  Float_t dDauXrowsTPC = dPosXrowsTPC; if (dDauXrowsTPC>dNegXrowsTPC) dDauXrowsTPC = dNegXrowsTPC;
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return 0x0;

  UShort_t wPosTPCNClsF = pDauPos->GetTPCNclsF(); if (wPosTPCNClsF<=0) return 0x0;
  UShort_t wNegTPCNClsF = pDauNeg->GetTPCNclsF(); if (wNegTPCNClsF<=0) return 0x0;
  Double_t dPosXrowsOverFindableClusTPC = ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF);
  Double_t dNegXrowsOverFindableClusTPC = ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF);

  Double_t dDauXrowsOverFindableClusTPC = dPosXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC>dNegXrowsOverFindableClusTPC) dDauXrowsOverFindableClusTPC = dNegXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return 0x0;
//=============================================================================

  Short_t nPosCharge = pDauPos->Charge();
  Short_t nNegCharge = pDauNeg->Charge();
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return 0x0;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPos = (AliAODTrack*)pV0->GetDaughter(1);
    pDauNeg = (AliAODTrack*)pV0->GetDaughter(0);

    dPosPxPyPz[0] = pV0->MomNegX(); dPosPxPyPz[1] = pV0->MomNegY(); dPosPxPyPz[2] = pV0->MomNegZ();
    dNegPxPyPz[0] = pV0->MomPosX(); dNegPxPyPz[1] = pV0->MomPosY(); dNegPxPyPz[2] = pV0->MomPosZ();
  } else {
    dPosPxPyPz[0] = pV0->MomPosX(); dPosPxPyPz[1] = pV0->MomPosY(); dPosPxPyPz[2] = pV0->MomPosZ();
    dNegPxPyPz[0] = pV0->MomNegX(); dNegPxPyPz[1] = pV0->MomNegY(); dNegPxPyPz[2] = pV0->MomNegZ();
  }

  TVector3 v3Pos(dPosPxPyPz);
  TVector3 v3Neg(dNegPxPyPz);

  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return 0x0;
  Double_t dPosEta = v3Pos.Eta(); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return 0x0;
  Double_t dNegEta = v3Neg.Eta(); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return 0x0;
//=============================================================================

  Bool_t bIsKshort = kFALSE;
  Bool_t bIsLambda = kFALSE;
  Bool_t bIsAntiLa = kFALSE;

  Float_t dPosPionSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kPion);
  Float_t dNegPionSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kPion);

  Float_t dPosProtonSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kProton);
  Float_t dNegProtonSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kProton);

  if (fCutMaxKshortSigmaTPC>0.) {
    bIsKshort = ((TMath::Abs(dPosPionSigmaTPC)<fCutMaxKshortSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)<fCutMaxKshortSigmaTPC));
  } else {
    bIsKshort = kTRUE;
  }

  if (fCutMaxLambdaSigmaTPC>0.) {
    bIsLambda = ((TMath::Abs(dPosProtonSigmaTPC)<fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)  <fCutMaxLambdaSigmaTPC));

    bIsAntiLa = ((TMath::Abs(dPosPionSigmaTPC)  <fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegProtonSigmaTPC)<fCutMaxLambdaSigmaTPC));
  } else {
    bIsLambda = kTRUE;
    bIsAntiLa = kTRUE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0CosPA = pV0->CosPointingAngle(fPrimaryVtx);

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0DistToPV = 0.;
  for (Int_t i=0; i<3; i++) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  Double_t dV0DistToPVoverP = TMath::Sqrt(dV0DistToPV) / (pV0->P()+1e-10);

  if (bIsKshort) if ((dV0DistToPVoverP*fgkMassKshort)>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*fgkMassLambda)>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0ArmFrac = pV0->PtArmV0() / (TMath::Abs(pV0->AlphaV0())+1e-12);

  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda || bIsAntiLa) && fCutMaxLambdaArmFrac>0.) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, fgkMassPion);
  TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, fgkMassPion);

  TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, fgkMassProton);
  TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, fgkMassProton);

  TLorentzVector vKshort = vPosPion   + vNegPion;
  TLorentzVector vLamvda = vPosProton + vNegPion;
  TLorentzVector vAntiLa = vNegProton + vPosPion;

  Double_t dKshortInvM = vKshort.M();
  Double_t dLambdaInvM = vLamvda.M();
  Double_t dAntiLaInvM = vAntiLa.M();

  if (bIsKshort) {
    Double_t dLower = 0.430006 - 0.0110029*dV0Pt;
    Double_t dUpper = 0.563707 + 0.0114979*dV0Pt;

    if ((dKshortInvM<dLower) || (dKshortInvM>dUpper)) bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) {
    Double_t dLower = 1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt);
    Double_t dUpper = 1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt);

    if (bIsLambda) if ((dLambdaInvM<dLower) || (dLambdaInvM>dUpper)) bIsLambda = kFALSE;
    if (bIsAntiLa) if ((dAntiLaInvM<dLower) || (dAntiLaInvM>dUpper)) bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  if (bIsKshort && (fCutMinKshortDeltaM>0.)) {
    if ((TMath::Abs(dLambdaInvM-fgkMassLambda)<fCutMinKshortDeltaM) ||
        (TMath::Abs(dAntiLaInvM-fgkMassLambda)<fCutMinKshortDeltaM)) bIsKshort = kFALSE;
  }

  if ((bIsLambda || bIsAntiLa) && (fCutMinLambdaDeletaM>0.)) {
    if ((TMath::Abs(dKshortInvM-fgkMassKshort)<fCutMinLambdaDeletaM)) {
      bIsLambda = kFALSE;
      bIsAntiLa = kFALSE;
    }
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  UInt_t wMask = 0;
  if (bIsKshort) wMask |= AliPicoHeaderCJ::kKshort;
  if (bIsLambda) wMask |= AliPicoHeaderCJ::kLambda;
  if (bIsAntiLa) wMask |= AliPicoHeaderCJ::kAntiLambda;

  Bool_t bPosInJC = kFALSE;
  Bool_t bNegInJC = kFALSE;

/*AliAODTrack *pTrkAOD = 0;
  Int_t idPos = pDauPos->GetID();
  Int_t idNeg = pDauNeg->GetID();
  if (fJetContisClArr) for (Int_t i=0; i<fJetContisClArr->GetEntriesFast(); i++) {
    pTrkAOD = static_cast<AliAODTrack*>(fJetContisClArr->At(i)); if (!pTrkAOD) continue;

    Int_t id = pTrkAOD->GetID();
    if (idPos==id) bPosInJC = kTRUE;
    if (idNeg==id) bNegInJC = kTRUE;
    if (bPosInJC && bNegInJC) { pTrkAOD = 0; break; }

    pTrkAOD = 0;
  }*/

  AliPicoV0RD *pPicoV0 = new AliPicoV0RD(wMask,
                                         dV0Radius,
                                         dV0CosPA,
                                         dV0DistToPVoverP,
                                         dDausDCA,
                                         dPosDCAtoPV,
                                         dNegDCAtoPV,
                                         dDauXrowsTPC,
                                         dDauXrowsOverFindableClusTPC,
                                         v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                                         v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                                         bPosInJC, bNegInJC,
                                         dPosPionSigmaTPC, dPosProtonSigmaTPC,
                                         dNegPionSigmaTPC, dNegProtonSigmaTPC);

  return pPicoV0;
}

//_____________________________________________________________________________
AliPicoV0RD* AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD(AliESDv0 const *pV0)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD
//

  if (pV0->GetOnFlyStatus()) return 0x0;
  if (pV0->GetChi2V0()>fCutMaxV0Chi2) return 0x0;
//=============================================================================

  Double_t dV0Pt = pV0->Pt(); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return 0x0;
  Double_t dKaRap = pV0->RapK0Short(); if ((dKaRap<fCutMinV0Rap) || (dKaRap>fCutMaxV0Rap)) return 0x0;
  Double_t dLaRap = pV0->RapLambda();  if ((dLaRap<fCutMinV0Rap) || (dLaRap>fCutMaxV0Rap)) return 0x0;
//=============================================================================

  Double_t dV0Vtx[3];  pV0->GetXYZ(dV0Vtx[0], dV0Vtx[1], dV0Vtx[2]);
  Double_t dV0Radius = TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]);
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return 0x0;

  Double_t dDausDCA = pV0->GetDcaV0Daughters(); if (dDausDCA>fCutMaxDausDCA) return 0x0;
//=============================================================================

  Int_t nPosIndex = pV0->GetPindex(); if (nPosIndex<0) return 0x0;
  Int_t nNegIndex = pV0->GetNindex(); if (nNegIndex<0) return 0x0;

  AliESDtrack *pDauPos = fEventESD->GetTrack(nPosIndex); if (!pDauPos) return 0x0;
  AliESDtrack *pDauNeg = fEventESD->GetTrack(nNegIndex); if (!pDauNeg) return 0x0;

  Double_t dMegField = fEventESD->GetMagneticField();
  Double_t dPosDCAtoPV = TMath::Abs(pDauPos->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)); if (dPosDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
  Double_t dNegDCAtoPV = TMath::Abs(pDauNeg->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)); if (dNegDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
//=============================================================================

  if (!(pDauPos->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;
  if (!(pDauNeg->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;
  if ((pDauPos->GetKinkIndex(0)>0) || (pDauNeg->GetKinkIndex(0)>0)) return 0x0;

  Float_t dPosXrowsTPC = pDauPos->GetTPCClusterInfo(2,1);
  Float_t dNegXrowsTPC = pDauNeg->GetTPCClusterInfo(2,1);
  Float_t dDauXrowsTPC = dPosXrowsTPC; if (dDauXrowsTPC>dNegXrowsTPC) dDauXrowsTPC = dNegXrowsTPC;
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return 0x0;

  UShort_t wPosTPCNClsF = pDauPos->GetTPCNclsF(); if (wPosTPCNClsF<=0) return 0x0;
  UShort_t wNegTPCNClsF = pDauNeg->GetTPCNclsF(); if (wNegTPCNClsF<=0) return 0x0;
  Double_t dPosXrowsOverFindableClusTPC = ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF);
  Double_t dNegXrowsOverFindableClusTPC = ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF);

  Double_t dDauXrowsOverFindableClusTPC = dPosXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC>dNegXrowsOverFindableClusTPC) dDauXrowsOverFindableClusTPC = dNegXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return 0x0;
//=============================================================================

  Short_t nPosCharge = pDauPos->Charge();
  Short_t nNegCharge = pDauNeg->Charge();
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return 0x0;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPos = fEventESD->GetTrack(nNegIndex);
    pDauNeg = fEventESD->GetTrack(nPosIndex);

    pV0->GetNPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0->GetPPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  } else {
    pV0->GetPPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0->GetNPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  }

  TVector3 v3Pos(dPosPxPyPz);
  TVector3 v3Neg(dNegPxPyPz);

  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return 0x0;
  Double_t dPosEta = v3Pos.Eta(); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return 0x0;
  Double_t dNegEta = v3Neg.Eta(); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return 0x0;
//=============================================================================

  Bool_t bIsKshort = kFALSE;
  Bool_t bIsLambda = kFALSE;
  Bool_t bIsAntiLa = kFALSE;

  Float_t dPosPionSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kPion);
  Float_t dNegPionSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kPion);

  Float_t dPosProtonSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kProton);
  Float_t dNegProtonSigmaTPC = fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kProton);

  if (fCutMaxKshortSigmaTPC>0.) {
    bIsKshort = ((TMath::Abs(dPosPionSigmaTPC)<fCutMaxKshortSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)<fCutMaxKshortSigmaTPC));
  } else {
    bIsKshort = kTRUE;
  }

  if (fCutMaxLambdaSigmaTPC>0.) {
    bIsLambda = ((TMath::Abs(dPosProtonSigmaTPC)<fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)  <fCutMaxLambdaSigmaTPC));

    bIsAntiLa = ((TMath::Abs(dPosPionSigmaTPC)  <fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegProtonSigmaTPC)<fCutMaxLambdaSigmaTPC));
  } else {
    bIsLambda = kTRUE;
    bIsAntiLa = kTRUE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0CosPA = pV0->GetV0CosineOfPointingAngle(fPrimaryVtx[0], fPrimaryVtx[1], fPrimaryVtx[2]);

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0DistToPV = 0.;
  for (Int_t i=0; i<3; i++) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  Double_t dV0DistToPVoverP = TMath::Sqrt(dV0DistToPV) / (pV0->P()+1e-10);

  if (bIsKshort) if ((dV0DistToPVoverP*fgkMassKshort)>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*fgkMassLambda)>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0ArmFrac = pV0->PtArmV0() / (TMath::Abs(pV0->AlphaV0())+1e-12);

  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda && bIsAntiLa) && (fCutMaxLambdaArmFrac>0.)) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, fgkMassPion);
  TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, fgkMassPion);

  TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, fgkMassProton);
  TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, fgkMassProton);

  TLorentzVector vKshort = vPosPion   + vNegPion;
  TLorentzVector vLamvda = vPosProton + vNegPion;
  TLorentzVector vAntiLa = vNegProton + vPosPion;

  Double_t dKshortInvM = vKshort.M();
  Double_t dLambdaInvM = vLamvda.M();
  Double_t dAntiLaInvM = vAntiLa.M();

  if (bIsKshort) {
    Double_t dLower = 0.430006 - 0.0110029*dV0Pt;
    Double_t dUpper = 0.563707 + 0.0114979*dV0Pt;

    if ((dKshortInvM<dLower) || (dKshortInvM>dUpper)) bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) {
    Double_t dLower = 1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt);
    Double_t dUpper = 1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt);

    if (bIsLambda) if ((dLambdaInvM<dLower) || (dLambdaInvM>dUpper)) bIsLambda = kFALSE;
    if (bIsAntiLa) if ((dAntiLaInvM<dLower) || (dAntiLaInvM>dUpper)) bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  if (bIsKshort && (fCutMinKshortDeltaM>0.)) {
    if ((TMath::Abs(dLambdaInvM-fgkMassLambda)<fCutMinKshortDeltaM) ||
        (TMath::Abs(dAntiLaInvM-fgkMassLambda)<fCutMinKshortDeltaM)) {
      bIsKshort = kFALSE;
    }
  }

  if ((bIsLambda || bIsAntiLa) && (fCutMinLambdaDeletaM>0.)) {
    if ((TMath::Abs(dKshortInvM-fgkMassKshort)<fCutMinLambdaDeletaM)) {
      bIsLambda = kFALSE;
      bIsAntiLa = kFALSE;
    }
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  UInt_t wMask = 0;
  if (bIsKshort) wMask |= AliPicoHeaderCJ::kKshort;
  if (bIsLambda) wMask |= AliPicoHeaderCJ::kLambda;
  if (bIsAntiLa) wMask |= AliPicoHeaderCJ::kAntiLambda;

  Bool_t bPosInJC = kFALSE;
  Bool_t bNegInJC = kFALSE;
  AliPicoV0RD *pPicoV0 = new AliPicoV0RD(wMask,
                                         dV0Radius,
                                         dV0CosPA,
                                         dV0DistToPVoverP,
                                         dDausDCA,
                                         dPosDCAtoPV,
                                         dNegDCAtoPV,
                                         dDauXrowsTPC,
                                         dDauXrowsOverFindableClusTPC,
                                         v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                                         v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                                         bPosInJC, bNegInJC,
                                         dPosPionSigmaTPC, dPosProtonSigmaTPC,
                                         dNegPionSigmaTPC, dNegProtonSigmaTPC);
  
  return pPicoV0;
}

//_____________________________________________________________________________
AliPicoV0MC* AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC(AliAODv0 const *pV0RD)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC
//

  if (pV0RD->GetOnFlyStatus()) return 0x0;
  if ((pV0RD->Chi2V0())>fCutMaxV0Chi2) return 0x0;

  Double_t dV0Pt = pV0RD->Pt(); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return 0x0;
//=============================================================================

  Double_t dV0Vtx[3]; pV0RD->GetXYZ(dV0Vtx);
  Double_t dV0Radius = TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]);
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return 0x0;

  Double_t dDausDCA = pV0RD->DcaV0Daughters(); if (dDausDCA>fCutMaxDausDCA) return 0x0;
  Double_t dPosDCAtoPV = pV0RD->DcaPosToPrimVertex(); if (dPosDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
  Double_t dNegDCAtoPV = pV0RD->DcaNegToPrimVertex(); if (dNegDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
//=============================================================================

  AliAODTrack *pDauPosRD = (AliAODTrack*)pV0RD->GetDaughter(0); if (!pDauPosRD) return 0x0;
  AliAODTrack *pDauNegRD = (AliAODTrack*)pV0RD->GetDaughter(1); if (!pDauNegRD) return 0x0;

  if (!(pDauPosRD->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;
  if (!(pDauNegRD->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;

  if ((pDauPosRD->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return 0x0;
  if ((pDauNegRD->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return 0x0;

  Float_t dPosXrowsTPC = pDauPosRD->GetTPCClusterInfo(2,1);
  Float_t dNegXrowsTPC = pDauNegRD->GetTPCClusterInfo(2,1);
  Float_t dDauXrowsTPC = dPosXrowsTPC; if (dDauXrowsTPC>dNegXrowsTPC) dDauXrowsTPC = dNegXrowsTPC;
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return 0x0;

  UShort_t wPosTPCNClsF = pDauPosRD->GetTPCNclsF(); if (wPosTPCNClsF<=0) return 0x0;
  UShort_t wNegTPCNClsF = pDauNegRD->GetTPCNclsF(); if (wNegTPCNClsF<=0) return 0x0;
  Double_t dPosXrowsOverFindableClusTPC = ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF);
  Double_t dNegXrowsOverFindableClusTPC = ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF);

  Double_t dDauXrowsOverFindableClusTPC = dPosXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC>dNegXrowsOverFindableClusTPC) dDauXrowsOverFindableClusTPC = dNegXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return 0x0;
//=============================================================================

  Short_t nPosCharge = pDauPosRD->Charge();
  Short_t nNegCharge = pDauNegRD->Charge();
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return 0x0;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPosRD = (AliAODTrack*)pV0RD->GetDaughter(1);
    pDauNegRD = (AliAODTrack*)pV0RD->GetDaughter(0);

    dPosPxPyPz[0] = pV0RD->MomNegX(); dPosPxPyPz[1] = pV0RD->MomNegY(); dPosPxPyPz[2] = pV0RD->MomNegZ();
    dNegPxPyPz[0] = pV0RD->MomPosX(); dNegPxPyPz[1] = pV0RD->MomPosY(); dNegPxPyPz[2] = pV0RD->MomPosZ();
  } else {
    dPosPxPyPz[0] = pV0RD->MomPosX(); dPosPxPyPz[1] = pV0RD->MomPosY(); dPosPxPyPz[2] = pV0RD->MomPosZ();
    dNegPxPyPz[0] = pV0RD->MomNegX(); dNegPxPyPz[1] = pV0RD->MomNegY(); dNegPxPyPz[2] = pV0RD->MomNegZ();
  }

  TVector3 v3Pos(dPosPxPyPz);
  TVector3 v3Neg(dNegPxPyPz);

  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return 0x0;
  Double_t dPosEta = v3Pos.Eta(); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return 0x0;
  Double_t dNegEta = v3Neg.Eta(); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return 0x0;
//=============================================================================

  Int_t inp = TMath::Abs(pDauPosRD->GetLabel()); if (inp<0) return 0x0;
  Int_t inn = TMath::Abs(pDauNegRD->GetLabel()); if (inn<0) return 0x0;
  AliAODMCParticle *pDauPosMC = (AliAODMCParticle*)MCEvent()->GetTrack(inp); if (!pDauPosMC) return 0x0;
  AliAODMCParticle *pDauNegMC = (AliAODMCParticle*)MCEvent()->GetTrack(inn); if (!pDauNegMC) return 0x0;

  Int_t imp = pDauPosMC->GetMother(); if (imp<0) return 0x0;
  Int_t imn = pDauNegMC->GetMother(); if (imn<0) return 0x0;
  if (imp != imn) return 0x0;

  AliAODMCParticle *pV0MC = (AliAODMCParticle*)MCEvent()->GetTrack(imp); if (!pV0MC) return 0x0;
  if (((pV0MC->Y())<fCutMinV0Rap) || ((pV0MC->Y())>fCutMaxV0Rap)) return 0x0;

  Int_t idvMC = pV0MC->GetPdgCode();
  Int_t idp = pDauPosMC->GetPdgCode();
  Int_t idn = pDauNegMC->GetPdgCode();
  Bool_t bIsKshort = ((idp==211)  && (idn==-211)  && (idvMC== 310));
  Bool_t bIsLambda = ((idp==2212) && (idn==-211)  && (idvMC== 3122));
  Bool_t bIsAntiLa = ((idp==211)  && (idn==-2212) && (idvMC==-3122));
  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  UInt_t wsvMC = 0;
  if (pV0MC->IsPrimary())                wsvMC |= AliPicoHeaderCJ::kPrimary;
  if (pV0MC->IsPhysicalPrimary())        wsvMC |= AliPicoHeaderCJ::kPhysicalPrimary;
  if (pV0MC->IsSecondaryFromWeakDecay()) wsvMC |= AliPicoHeaderCJ::kSecondaryFromWeakDecay;
  if (pV0MC->IsSecondaryFromMaterial())  wsvMC |= AliPicoHeaderCJ::kSecondaryFromMaterial;

  Int_t  idmMC = 0;
  UInt_t wsmMC = 0;
  Double_t dMotherPt  = 0.;
  Double_t dMotherEta = 0.;
  Double_t dMotherRap = 0.;
  if (bIsLambda || bIsAntiLa) {
    Int_t imv = pV0MC->GetMother(); if (imv>=0) {
      AliAODMCParticle *pMother = (AliAODMCParticle*)MCEvent()->GetTrack(imv);

      if (pMother) {
        idmMC = pMother->GetPdgCode();
        if ((bIsLambda && ((idmMC== 3312) || (idmMC== 3322))) ||
            (bIsAntiLa && ((idmMC==-3312) || (idmMC==-3322)))) {
          dMotherPt  = pMother->Pt();
          dMotherEta = pMother->Eta();
          dMotherRap = pMother->Y();

          if (pMother->IsPrimary())                wsmMC |= AliPicoHeaderCJ::kPrimary;
          if (pMother->IsPhysicalPrimary())        wsmMC |= AliPicoHeaderCJ::kPhysicalPrimary;
          if (pMother->IsSecondaryFromWeakDecay()) wsmMC |= AliPicoHeaderCJ::kSecondaryFromWeakDecay;
          if (pMother->IsSecondaryFromMaterial())  wsmMC |= AliPicoHeaderCJ::kSecondaryFromMaterial;
        }
      }
    }
  }
//=============================================================================

  Double_t dV0CosPA = pV0RD->CosPointingAngle(fPrimaryVtx);

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0DistToPV = 0.;
  for (Int_t i=0; i<3; i++) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  Double_t dV0DistToPVoverP = TMath::Sqrt(dV0DistToPV) / (pV0RD->P()+1e-10);

  if (bIsKshort) if ((dV0DistToPVoverP*fgkMassKshort)>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*fgkMassLambda)>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0ArmFrac = pV0RD->PtArmV0() / (TMath::Abs(pV0RD->AlphaV0())+1e-12);

  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda || bIsAntiLa) && fCutMaxLambdaArmFrac>0.) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Int_t wMask = 0;
  if (bIsKshort) {
    TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, fgkMassPion);
    TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, fgkMassPion);
    TLorentzVector vKshort = vPosPion + vNegPion;

    Double_t dKshortInvM = vKshort.M();
    Double_t dLower = 0.430006 - 0.0110029*dV0Pt;
    Double_t dUpper = 0.563707 + 0.0114979*dV0Pt;
    if ((dKshortInvM<dLower) || (dKshortInvM>dUpper)) return 0x0;

    if (fCutMinKshortDeltaM>0.) {
      TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, fgkMassProton);
      TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, fgkMassProton);

      TLorentzVector vLamvda = vPosProton + vNegPion;
      TLorentzVector vAntiLa = vNegProton + vPosPion;

      Double_t dLambdaInvM = vLamvda.M();
      Double_t dAntiLaInvM = vAntiLa.M();
      if ((TMath::Abs(dLambdaInvM-fgkMassLambda)<fCutMinKshortDeltaM) ||
          (TMath::Abs(dAntiLaInvM-fgkMassLambda)<fCutMinKshortDeltaM)) return 0x0;
    }

    wMask = AliPicoHeaderCJ::kKshort;
  }

  if (bIsLambda) {
    TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, fgkMassProton);
    TLorentzVector vNegPion;   vNegPion.SetVectM(v3Neg, fgkMassPion);
    TLorentzVector vLamvda = vPosProton + vNegPion;

    Double_t dLambdaInvM = vLamvda.M();
    Double_t dLower = 1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt);
    Double_t dUpper = 1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt);
    if ((dLambdaInvM<dLower) || (dLambdaInvM>dUpper)) return 0x0;

    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, fgkMassPion);
      TLorentzVector vKshort = vPosPion + vNegPion;

      Double_t dKshortInvM = vKshort.M();
      if ((TMath::Abs(dKshortInvM-fgkMassKshort)<fCutMinLambdaDeletaM)) return 0x0;
    }

    wMask = AliPicoHeaderCJ::kLambda;
  }

  if (bIsAntiLa) {
    TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, fgkMassProton);
    TLorentzVector vPosPion;   vPosPion.SetVectM(v3Pos, fgkMassPion);
    TLorentzVector vAntiLa = vNegProton + vPosPion;

    Double_t dAntiLaInvM = vAntiLa.M();
    Double_t dLower = 1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt);
    Double_t dUpper = 1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt);
    if ((dAntiLaInvM<dLower) || (dAntiLaInvM>dUpper)) return 0x0;

    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, fgkMassPion);
      TLorentzVector vKshort = vPosPion + vNegPion;

      Double_t dKshortInvM = vKshort.M();
      if ((TMath::Abs(dKshortInvM-fgkMassKshort)<fCutMinLambdaDeletaM)) return 0x0;
    }

    wMask = AliPicoHeaderCJ::kAntiLambda;
  }
//=============================================================================

  Bool_t bPosInJC = kFALSE;
  Bool_t bNegInJC = kFALSE;
  AliPicoV0MC *pPicoV0 = new AliPicoV0MC(wMask,
                                         dV0Radius,
                                         dV0CosPA,
                                         dV0DistToPVoverP,
                                         dDausDCA,
                                         dPosDCAtoPV,
                                         dNegDCAtoPV,
                                         dDauXrowsTPC,
                                         dDauXrowsOverFindableClusTPC,
                                         v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                                         v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                                         bPosInJC, bNegInJC,
                                         idvMC, wsvMC, pV0MC->Px(), pV0MC->Py(), pV0MC->Pz(), pV0MC->E(),
                                         idmMC, wsmMC, dMotherPt, dMotherEta, dMotherRap);


  return pPicoV0;
}

//_____________________________________________________________________________
AliPicoV0MC* AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC(AliESDv0 const *pV0RD)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC
//

  AliStack *pStack = MCEvent()->Stack(); if (!pStack) return 0x0;
  Int_t nPrimary = pStack->GetNprimary();
//=============================================================================

  if (pV0RD->GetOnFlyStatus()) return 0x0;
  if (pV0RD->GetChi2V0()>fCutMaxV0Chi2) return 0x0;

  Double_t dV0Pt = pV0RD->Pt(); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return 0x0;
//=============================================================================
  
  Double_t dV0Vtx[3];  pV0RD->GetXYZ(dV0Vtx[0], dV0Vtx[1], dV0Vtx[2]);
  Double_t dV0Radius = TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]);
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return 0x0; 

  Double_t dDausDCA = pV0RD->GetDcaV0Daughters(); if (dDausDCA>fCutMaxDausDCA) return 0x0;
//=============================================================================

  Int_t nPosIndex = TMath::Abs(pV0RD->GetPindex()); if (nPosIndex<0) return 0x0;
  Int_t nNegIndex = TMath::Abs(pV0RD->GetNindex()); if (nNegIndex<0) return 0x0;

  AliESDtrack *pDauPosRD = fEventESD->GetTrack(nPosIndex); if (!pDauPosRD) return 0x0;
  AliESDtrack *pDauNegRD = fEventESD->GetTrack(nNegIndex); if (!pDauNegRD) return 0x0;

  Double_t dMegField = fEventESD->GetMagneticField();
  Double_t dPosDCAtoPV = TMath::Abs(pDauPosRD->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)); if (dPosDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
  Double_t dNegDCAtoPV = TMath::Abs(pDauNegRD->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)); if (dNegDCAtoPV<fCutMinDauDCAtoPV) return 0x0;
//=============================================================================

  if (!(pDauPosRD->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;
  if (!(pDauNegRD->GetStatus() & AliESDtrack::kTPCrefit)) return 0x0;
  if ((pDauPosRD->GetKinkIndex(0)>0) || (pDauNegRD->GetKinkIndex(0)>0)) return 0x0;

  Float_t dPosXrowsTPC = pDauPosRD->GetTPCClusterInfo(2,1);
  Float_t dNegXrowsTPC = pDauNegRD->GetTPCClusterInfo(2,1);
  Float_t dDauXrowsTPC = dPosXrowsTPC; if (dDauXrowsTPC>dNegXrowsTPC) dDauXrowsTPC = dNegXrowsTPC;
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return 0x0;

  UShort_t wPosTPCNClsF = pDauPosRD->GetTPCNclsF(); if (wPosTPCNClsF<=0) return 0x0;
  UShort_t wNegTPCNClsF = pDauNegRD->GetTPCNclsF(); if (wNegTPCNClsF<=0) return 0x0;
  Double_t dPosXrowsOverFindableClusTPC = ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF);
  Double_t dNegXrowsOverFindableClusTPC = ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF);

  Double_t dDauXrowsOverFindableClusTPC = dPosXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC>dNegXrowsOverFindableClusTPC) dDauXrowsOverFindableClusTPC = dNegXrowsOverFindableClusTPC;
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return 0x0;
//=============================================================================

  Short_t nPosCharge = pDauPosRD->Charge();
  Short_t nNegCharge = pDauNegRD->Charge();
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return 0x0;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPosRD = fEventESD->GetTrack(nNegIndex);
    pDauNegRD = fEventESD->GetTrack(nPosIndex);

    pV0RD->GetNPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0RD->GetPPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  } else {
    pV0RD->GetPPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0RD->GetNPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  }

  TVector3 v3Pos(dPosPxPyPz);
  TVector3 v3Neg(dNegPxPyPz);

  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return 0x0;
  Double_t dPosEta = v3Pos.Eta(); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return 0x0;
  Double_t dNegEta = v3Neg.Eta(); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return 0x0;
//=============================================================================

  Int_t inp = TMath::Abs(pDauPosRD->GetLabel()); if (inp<0) return 0x0;
  Int_t inn = TMath::Abs(pDauNegRD->GetLabel()); if (inn<0) return 0x0;
  TParticle *pDauPosMC = ((AliMCParticle*)MCEvent()->GetTrack(inp))->Particle(); if (!pDauPosMC) return 0x0;
  TParticle *pDauNegMC = ((AliMCParticle*)MCEvent()->GetTrack(inn))->Particle(); if (!pDauNegMC) return 0x0;

  Int_t imp = pDauPosMC->GetFirstMother(); if (imp<0) return 0x0;
  Int_t imn = pDauNegMC->GetFirstMother(); if (imn<0) return 0x0;
  if (imp != imn) return 0x0;

  TParticle *pV0MC = ((AliMCParticle*)MCEvent()->GetTrack(imp))->Particle(); if (!pV0MC) return 0x0;
  if (((pV0MC->Y())<fCutMinV0Rap) || ((pV0MC->Y())>fCutMaxV0Rap)) return 0x0;

  Int_t idvMC = pV0MC->GetPdgCode();
  Int_t idp = pDauPosMC->GetPdgCode();
  Int_t idn = pDauNegMC->GetPdgCode();
  Bool_t bIsKshort = ((idp==211)  && (idn==-211)  && (idvMC== 310));
  Bool_t bIsLambda = ((idp==2212) && (idn==-211)  && (idvMC== 3122));
  Bool_t bIsAntiLa = ((idp==211)  && (idn==-2212) && (idvMC==-3122));
  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  UInt_t wsvMC = 0;
  if (imp<nPrimary)                          wsvMC |= AliPicoHeaderCJ::kPrimary;
  if (pStack->IsPhysicalPrimary(imp))        wsvMC |= AliPicoHeaderCJ::kPhysicalPrimary;
  if (pStack->IsSecondaryFromWeakDecay(imp)) wsvMC |= AliPicoHeaderCJ::kSecondaryFromWeakDecay;
  if (pStack->IsSecondaryFromMaterial(imp))  wsvMC |= AliPicoHeaderCJ::kSecondaryFromMaterial;

  Int_t  idmMC = 0;
  UInt_t wsmMC = 0;
  Double_t dMotherPt  = 0.;
  Double_t dMotherEta = 0.;
  Double_t dMotherRap = 0.;
  if (bIsLambda || bIsAntiLa) {
    Int_t imv = pV0MC->GetFirstMother(); if (imv>=0) {
      TParticle *pMother = ((AliMCParticle*)MCEvent()->GetTrack(imv))->Particle();

      if (pMother) {
        idmMC = pMother->GetPdgCode();
        if ((bIsLambda && ((idmMC== 3312) || (idmMC== 3322))) ||
            (bIsAntiLa && ((idmMC==-3312) || (idmMC==-3322)))) {
          dMotherPt  = pMother->Pt();
          dMotherEta = pMother->Eta();
          dMotherRap = pMother->Y();

          if (imp<nPrimary)                          wsmMC |= AliPicoHeaderCJ::kPrimary;
          if (pStack->IsPhysicalPrimary(imv))        wsmMC |= AliPicoHeaderCJ::kPhysicalPrimary;
          if (pStack->IsSecondaryFromWeakDecay(imv)) wsmMC |= AliPicoHeaderCJ::kSecondaryFromWeakDecay;
          if (pStack->IsSecondaryFromMaterial(imv))  wsmMC |= AliPicoHeaderCJ::kSecondaryFromMaterial;
        }
      }
    }
  }
//=============================================================================

  Double_t dV0CosPA = pV0RD->GetV0CosineOfPointingAngle(fPrimaryVtx[0], fPrimaryVtx[1], fPrimaryVtx[2]);

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0DistToPV = 0.;
  for (Int_t i=0; i<3; i++) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  Double_t dV0DistToPVoverP = TMath::Sqrt(dV0DistToPV) / (pV0RD->P()+1e-10);

  if (bIsKshort) if ((dV0DistToPVoverP*fgkMassKshort)>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*fgkMassLambda)>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Double_t dV0ArmFrac = pV0RD->PtArmV0() / (TMath::Abs(pV0RD->AlphaV0())+1e-12);

  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda && bIsAntiLa) && (fCutMaxLambdaArmFrac>0.)) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return 0x0;
//=============================================================================

  Int_t wMask = 0;

  if (bIsKshort) {
    TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, fgkMassPion);
    TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, fgkMassPion);
    TLorentzVector vKshort = vPosPion + vNegPion;

    Double_t dKshortInvM = vKshort.M();
    Double_t dLower = 0.430006 - 0.0110029*dV0Pt;
    Double_t dUpper = 0.563707 + 0.0114979*dV0Pt;
    if ((dKshortInvM<dLower) || (dKshortInvM>dUpper)) return 0x0;

    if (fCutMinKshortDeltaM>0.) {
      TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, fgkMassProton);
      TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, fgkMassProton);

      TLorentzVector vLamvda = vPosProton + vNegPion;
      TLorentzVector vAntiLa = vNegProton + vPosPion;

      Double_t dLambdaInvM = vLamvda.M();
      Double_t dAntiLaInvM = vAntiLa.M();
      if ((TMath::Abs(dLambdaInvM-fgkMassLambda)<fCutMinKshortDeltaM) ||
          (TMath::Abs(dAntiLaInvM-fgkMassLambda)<fCutMinKshortDeltaM)) return 0x0;
    }

    wMask = AliPicoHeaderCJ::kKshort;
  }

  if (bIsLambda) {
    TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, fgkMassProton);
    TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, fgkMassPion);
    TLorentzVector vLamvda = vPosProton + vNegPion;

    Double_t dLambdaInvM = vLamvda.M();
    Double_t dLower = 1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt);
    Double_t dUpper = 1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt);
    if ((dLambdaInvM<dLower) || (dLambdaInvM>dUpper)) return 0x0;

    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, fgkMassPion);
      TLorentzVector vKshort = vPosPion + vNegPion;

      Double_t dKshortInvM = vKshort.M();
      if ((TMath::Abs(dKshortInvM-fgkMassKshort)<fCutMinLambdaDeletaM)) return 0x0;
    }

    wMask = AliPicoHeaderCJ::kLambda;
  }

  if (bIsAntiLa) {
    TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, fgkMassProton);
    TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, fgkMassPion);
    TLorentzVector vAntiLa = vNegProton + vPosPion;
  
    Double_t dAntiLaInvM = vAntiLa.M();
    Double_t dLower = 1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt);
    Double_t dUpper = 1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt);
    if ((dAntiLaInvM<dLower) || (dAntiLaInvM>dUpper)) return 0x0;
  
    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, fgkMassPion);
      TLorentzVector vKshort = vPosPion + vNegPion;

      Double_t dKshortInvM = vKshort.M();
      if ((TMath::Abs(dKshortInvM-fgkMassKshort)<fCutMinLambdaDeletaM)) return 0x0;
    }

    wMask = AliPicoHeaderCJ::kAntiLambda;
  }
//=============================================================================

  Bool_t bPosInJC = kFALSE;
  Bool_t bNegInJC = kFALSE;
  AliPicoV0MC *pPicoV0 = new AliPicoV0MC(wMask,
                                         dV0Radius,
                                         dV0CosPA,
                                         dV0DistToPVoverP,
                                         dDausDCA,
                                         dPosDCAtoPV,
                                         dNegDCAtoPV,
                                         dDauXrowsTPC,
                                         dDauXrowsOverFindableClusTPC,
                                         v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                                         v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                                         bPosInJC, bNegInJC,
                                         idvMC, wsvMC, pV0MC->Px(), pV0MC->Py(), pV0MC->Pz(), pV0MC->Energy(),
                                         idmMC, wsmMC, dMotherPt, dMotherEta, dMotherRap);


  return pPicoV0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEPicoV0Maker::IsEventNotAcpt()
{
//
//  AliAnalysisTaskSEPicoV0Maker::IsEventNotAcpt
//

  fEventAcptMask = 0;
  if (!InputEvent())  return (fEventAcptMask==0);
  if (!fInputHandler) return (fEventAcptMask==0);

  if (fCollisionType!=(AliPicoHeaderCJ::kPP)) {
    fCentInfo = InputEvent()->GetCentrality();
    if (!fCentInfo) return (fEventAcptMask==0);
  }

  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fEventESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if ((!fEventAOD) && (!fEventESD)) return (fEventAcptMask==0);

  fRespoPID = fInputHandler->GetPIDResponse();
  if (!fRespoPID) return kTRUE;

  if (fIsAnaInfoMC) {
     if (MCEvent()) {
      if (MCEvent()->GetNumberOfTracks()<=0) return (fEventAcptMask==0);
    } else return (fEventAcptMask==0);

    AliHeader *pHeader = MCEvent()->Header(); if (!pHeader) return (fEventAcptMask==0);

    if (fIsDPMjetMC) {
      AliGenDPMjetEventHeader *pDPMjetH = dynamic_cast<AliGenDPMjetEventHeader*>(pHeader->GenEventHeader());

      if (pDPMjetH) {
        Int_t nd0=0, nd1=0, nd2=0; pDPMjetH->GetNDiffractive(nd1, nd2, nd0);
        if ((nd1+nd2) != (pDPMjetH->ProjectileParticipants() + pDPMjetH->TargetParticipants())) return (fEventAcptMask==0);
      }
    }
  }

  fEventAcptMask |= AliPicoHeaderCJ::kEventAccCheck;
//=============================================================================

  if (fCollisionType==(AliPicoHeaderCJ::kPP)) {
    fEventAcptMask |= AliPicoHeaderCJ::kEventAccMult;
  } else {
    if (fCentInfo->GetQuality()==0)
      fEventAcptMask |= AliPicoHeaderCJ::kEventAccMult;
    else
      return (fEventAcptMask==0);
  }
//=============================================================================

  AliBits wMask = fInputHandler->IsEventSelected();
  if ((wMask & fTriggerMask) != fTriggerMask) return (fEventAcptMask==0);
  if (fIsSkipFastOnly) if ((wMask & AliVEvent::kFastOnly) == AliVEvent::kFastOnly) return (fEventAcptMask==0);

  fEventAcptMask |= AliPicoHeaderCJ::kEventAccTrigger;
//=============================================================================

  const AliVVertex *pVertex = InputEvent()->GetPrimaryVertex(); if (!pVertex) return (fEventAcptMask==0);
  pVertex->GetXYZ(fPrimaryVtx); if (TMath::Abs(fPrimaryVtx[2])>fCutMaxEventVzAbs) return (fEventAcptMask==0);

  if ((fCollisionType==(AliPicoHeaderCJ::kPA)) || (fCollisionType==(AliPicoHeaderCJ::kAP))) {
    if ( fAnaUtils->IsFirstEventInChunk(InputEvent()))    return (fEventAcptMask==0);
    if (!fAnaUtils->IsVertexSelected2013pA(InputEvent())) return (fEventAcptMask==0);

/*  if (fEventAOD) {
      const AliAODVertex *pVtxSPD = fEventAOD->GetPrimaryVertexSPD();
      const AliAODVertex *pVtxTrk = fEventAOD->GetPrimaryVertex();
      if ((!pVtxSPD) && (!pVtxTrk)) return (fEventAcptMask==0);
    }*/

/*  if (fEventESD) {
      Bool_t fHasVertex = kFALSE;
      const AliESDVertex *pVtxESD = fEventESD->GetPrimaryVertexTracks();
      if (pVtxESD->GetNContributors()<1) {
        pVtxESD = fEventESD->GetPrimaryVertexSPD();
        if (pVtxESD->GetNContributors()<1) fHasVertex = kFALSE;
        else fHasVertex = kTRUE;

        Double_t cov[6] = { 0., 0., 0., 0., 0., 0. };
        pVtxESD->GetCovarianceMatrix(cov);
        Double_t zRes = TMath::Sqrt(cov[5]);
        if (pVtxESD->IsFromVertexerZ() && (zRes>0.25)) fHasVertex = kFALSE;
      } else fHasVertex = kTRUE;

      if (!fHasVertex) return (fEventAcptMask==0);
    }*/

  } else {
    if (fEventAOD) {
      const AliAODVertex *pVtxSPD = fEventAOD->GetPrimaryVertexSPD(); if (!pVtxSPD) return (fEventAcptMask==0);
      const AliAODVertex *pVtxTrk = fEventAOD->GetPrimaryVertex();    if (!pVtxTrk) return (fEventAcptMask==0);
    }

    if (fEventESD) {
      const AliESDVertex *pVtxPri = fEventESD->GetPrimaryVertex();       if (!pVtxPri) return (fEventAcptMask==0);
      const AliESDVertex *pVtxSPD = fEventESD->GetPrimaryVertexSPD();    if (!pVtxSPD) return (fEventAcptMask==0);
      const AliESDVertex *pVtxTrk = fEventESD->GetPrimaryVertexTracks(); if (!pVtxTrk) return (fEventAcptMask==0);
      if ((!(pVtxPri->GetStatus())) && (!(pVtxSPD->GetStatus())) && (!(pVtxTrk->GetStatus()))) return (fEventAcptMask==0);
    }
  }

  fEventAcptMask |= AliPicoHeaderCJ::kEventAccVertex;
//=============================================================================

  if ((fCollisionType==AliPicoHeaderCJ::kPP) ||
      (fCollisionType==AliPicoHeaderCJ::kPA) ||
      (fCollisionType==AliPicoHeaderCJ::kAP)) {
    if (fAnaUtils->IsPileUpEvent(InputEvent())) return (fEventAcptMask==0);
  }

  fEventAcptMask |= AliPicoHeaderCJ::kEventAccPileup;
//=============================================================================

  if (fIsRefitV0sESD && fEventESD) {
    Double_t dCuts[7] = { fCutMaxV0Chi2,
                          fCutMinDauDCAtoPV,
                          fCutMinDauDCAtoPV,
                          fCutMaxDausDCA,
                          fCutMinKshortCosPA,
                          fCutMinV0Radius,
                          fCutMaxV0Radius };

    fEventESD->ResetV0s();
    AliV0vertexer aV0vtxer;
    aV0vtxer.SetDefaultCuts(dCuts);
    aV0vtxer.Tracks2V0vertices(fEventESD);
  }

  return (fEventAcptMask==0);
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEPicoV0Maker::IsEventNotINEL()
{
//
//  AliAnalysisTaskSEPicoV0Maker::IsEventNotINEL
//

  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccCheck) != AliPicoHeaderCJ::kEventAccCheck) return kTRUE;
  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccMult)  != AliPicoHeaderCJ::kEventAccMult)  return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEPicoV0Maker::IsEventNotMBsa()
{
//
//  AliAnalysisTaskSEPicoV0Maker::IsEventNotMBsa
//

  if (IsEventNotINEL()) return kTRUE;

  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccTrigger) != AliPicoHeaderCJ::kEventAccTrigger) return kTRUE;
  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccVertex)  != AliPicoHeaderCJ::kEventAccVertex)  return kTRUE;
  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccPileup)  != AliPicoHeaderCJ::kEventAccPileup)  return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::FillHistogramsEH()
{
//
//  AliAnalysisTaskSEPicoV0Maker::FillHistogramsEH
//

  Float_t dV0M = fCentInfo->GetCentralityPercentile("V0M");
  Float_t dV0A = fCentInfo->GetCentralityPercentile("V0A");
  Float_t dCL1 = fCentInfo->GetCentralityPercentile("CL1");
  Float_t dZNA = fCentInfo->GetCentralityPercentile("ZNA");

  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccCheck) == AliPicoHeaderCJ::kEventAccCheck) {
    ((TH1D*)fOutputListEH->FindObject("hEventAccCheck_V0M"))->Fill(dV0M);
    ((TH1D*)fOutputListEH->FindObject("hEventAccCheck_V0A"))->Fill(dV0A);
    ((TH1D*)fOutputListEH->FindObject("hEventAccCheck_CL1"))->Fill(dCL1);
    ((TH1D*)fOutputListEH->FindObject("hEventAccCheck_ZNA"))->Fill(dZNA);
  }

  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccMult) == AliPicoHeaderCJ::kEventAccMult) {
    ((TH1D*)fOutputListEH->FindObject("hEventAccMult_V0M"))->Fill(dV0M);
    ((TH1D*)fOutputListEH->FindObject("hEventAccMult_V0A"))->Fill(dV0A);
    ((TH1D*)fOutputListEH->FindObject("hEventAccMult_CL1"))->Fill(dCL1);
    ((TH1D*)fOutputListEH->FindObject("hEventAccMult_ZNA"))->Fill(dZNA);
  }

  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccTrigger) == AliPicoHeaderCJ::kEventAccTrigger) {
    ((TH1D*)fOutputListEH->FindObject("hEventAccTrigger_V0M"))->Fill(dV0M);
    ((TH1D*)fOutputListEH->FindObject("hEventAccTrigger_V0A"))->Fill(dV0A);
    ((TH1D*)fOutputListEH->FindObject("hEventAccTrigger_CL1"))->Fill(dCL1);
    ((TH1D*)fOutputListEH->FindObject("hEventAccTrigger_ZNA"))->Fill(dZNA);
  }

  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccVertex) == AliPicoHeaderCJ::kEventAccVertex) {
    ((TH1D*)fOutputListEH->FindObject("hEventAccVertex_V0M"))->Fill(dV0M);
    ((TH1D*)fOutputListEH->FindObject("hEventAccVertex_V0A"))->Fill(dV0A);
    ((TH1D*)fOutputListEH->FindObject("hEventAccVertex_CL1"))->Fill(dCL1);
    ((TH1D*)fOutputListEH->FindObject("hEventAccVertex_ZNA"))->Fill(dZNA);
  }

  if ((fEventAcptMask & AliPicoHeaderCJ::kEventAccPileup) == AliPicoHeaderCJ::kEventAccPileup) {
    ((TH1D*)fOutputListEH->FindObject("hEventAccPileup_V0M"))->Fill(dV0M);
    ((TH1D*)fOutputListEH->FindObject("hEventAccPileup_V0A"))->Fill(dV0A);
    ((TH1D*)fOutputListEH->FindObject("hEventAccPileup_CL1"))->Fill(dCL1);
    ((TH1D*)fOutputListEH->FindObject("hEventAccPileup_ZNA"))->Fill(dZNA);
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::FillHistogramsMC()
{
//
//  AliAnalysisTaskSEPicoV0Maker::FillHistogramsMC
//

  Int_t   nPrimary = 0;
  AliStack *pStack = 0;

  if (fEventESD) {
    pStack   = MCEvent()->Stack(); if (!pStack) return;
    nPrimary = pStack->GetNprimary();
  }
//=============================================================================

  Double_t dEvType = -0.5;
  if (IsEventNotMBsa()) dEvType = 0.5;

  Double_t dV0M = fCentInfo->GetCentralityPercentile("V0M");
  Double_t dV0A = fCentInfo->GetCentralityPercentile("V0A");
  Double_t dCL1 = fCentInfo->GetCentralityPercentile("CL1");
  Double_t dZNA = fCentInfo->GetCentralityPercentile("ZNA");

  THnSparseD *hsV0 = dynamic_cast<THnSparseD*>(fOutputListMC->FindObject("hsV0"));
  THnSparseD *hsXi = dynamic_cast<THnSparseD*>(fOutputListMC->FindObject("hsXi"));

  if (hsV0 == 0 || hsXi == 0) { // Keep Coverity happy
    AliFatal("Cannot find hsV0 or hsXi; should not happen");
    return;
  }  // Keep Coverity happy
//=============================================================================

  TParticle        *pESD = 0;
  AliAODMCParticle *pAOD = 0;
  for (Int_t i=0; i<MCEvent()->GetNumberOfTracks(); i++) {
    if (fEventAOD) { pAOD = (AliAODMCParticle*)MCEvent()->GetTrack(i);              if (!pAOD) continue; }
    if (fEventESD) { pESD =   ((AliMCParticle*)MCEvent()->GetTrack(i))->Particle(); if (!pESD) continue; }

    Bool_t bPri = kFALSE;
    if (pAOD) bPri = pAOD->IsPrimary();
    if (pESD) bPri = (i<nPrimary);

    Bool_t bPhy = kFALSE;
    if (pAOD) bPhy =   pAOD->IsPhysicalPrimary();
    if (pESD) bPhy = pStack->IsPhysicalPrimary(i);
    if ((!bPri) && (!bPhy)) { pAOD=0; pESD=0; continue; }

    Int_t id = 0;
    if (pAOD) id = pAOD->GetPdgCode();
    if (pESD) id = pESD->GetPdgCode();

    Bool_t bXi = (bPri && ((id==3312) || (id==-3312)));
    Bool_t bV0 = (bPhy && ((id==3122) || (id==-3122) || (id==310)));
    if (!(bXi || bV0)) { pAOD=0; pESD=0; continue; }

    Double_t  dEta = 0.;
    if (pAOD) dEta = pAOD->Eta();
    if (pESD) dEta = pESD->Eta();
    if ((dEta<-5.) || (dEta>=5.)) { pAOD=0; pESD=0; continue; }

    Double_t  dRapLab = 0.;
    if (pAOD) dRapLab = pAOD->Y();
    if (pESD) dRapLab = pESD->Y();
    if ((dRapLab<-5.) || (dRapLab>=5.)) { pAOD=0; pESD=0; continue; }

    Double_t dRapCMS = dRapLab + fRapidityShift;
    if ((dRapCMS<-5.) || (dRapCMS>=5.)) { pAOD=0; pESD=0; continue; }

    Double_t dVar[10];
    if (pAOD) dVar[9] = pAOD->Pt();
    if (pESD) dVar[9] = pESD->Pt();

    dVar[1] = dEvType;
    dVar[2] = dV0M;
    dVar[3] = dV0A;
    dVar[4] = dCL1;
    dVar[5] = dZNA;
    dVar[6] = dEta;
    dVar[7] = dRapLab;
    dVar[8] = dRapCMS;

    if (bXi) {
      if (id== 3312) dVar[0] = -0.5;
      if (id==-3312) dVar[0] =  0.5;
      hsXi->Fill(dVar);
    }

    if (bV0) {
      if (id== 310 ) dVar[0] = 0.;
      if (id== 3122) dVar[0] = 1.;
      if (id==-3122) dVar[0] = 2.;
      hsV0->Fill(dVar);
    }

    pAOD = 0;
    pESD = 0;
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::CreateHistogramsEH()
{
//
//  AliAnalysisTaskSEPicoV0Maker::CreateHistogramsEH
//

  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
//=============================================================================

  TH1D *h1 = 0;
  h1 = new TH1D("hEventAccCheck_V0M", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccCheck_V0A", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccCheck_CL1", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccCheck_ZNA", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);

  h1 = new TH1D("hEventAccMult_V0M", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccMult_V0A", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccMult_CL1", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccMult_ZNA", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);

  h1 = new TH1D("hEventAccTrigger_V0M", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccTrigger_V0A", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccTrigger_CL1", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccTrigger_ZNA", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);

  h1 = new TH1D("hEventAccVertex_V0M", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccVertex_V0A", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccVertex_CL1", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccVertex_ZNA", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);

  h1 = new TH1D("hEventAccPileup_V0M", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccPileup_V0A", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccPileup_CL1", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);
  h1 = new TH1D("hEventAccPileup_ZNA", "", 210, -10., 200.); h1->Sumw2(); fOutputListEH->Add(h1);


  TH2D *h2 = 0;
  h2 = new TH2D("hKshortPtInvM", "", 1000, 0., 100., 300, fgkMassKshort-0.15, fgkMassKshort+0.15);
  h2->Sumw2(); fOutputListEH->Add(h2); h2=0;

  h2 = new TH2D("hLambdaPtInvM", "", 1000, 0., 100., 200, fgkMassLambda-0.10, fgkMassLambda+0.10);
  h2->Sumw2(); fOutputListEH->Add(h2); h2=0;

  h2 = new TH2D("hAntiLaPtInvM", "", 1000, 0., 100., 200, fgkMassLambda-0.10, fgkMassLambda+0.10);
  h2->Sumw2(); fOutputListEH->Add(h2); h2=0;

  TH1::AddDirectory(bStatusTmpH);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::CreateHistogramsMC()
{
//
//  AliAnalysisTaskSEPicoV0Maker::CreateHistogramsMC
//

// TODO: add the multiplicity bins

  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nV0 = 10; // 0: particle type
                        //    ==0, Kshort
                        //    ==1, Lambda
                        //    ==2, AntiLa
                        // 1: Event type
                        //    ==-0.5, INEL
                        //    == 0.5, MB
                        // 2: V0M
                        // 3: V0A
                        // 4: CL1
                        // 5: ZNA
                        // 6: eta
                        // 7: rap in Lab
                        // 8: rap in CMS
                        // 9: Pt
  const Int_t    nV0Bin[nV0] = {  3,    2,  210,  210,  210,  210, 100, 100, 100, 1000  };
  const Double_t dV0Min[nV0] = { -0.5, -1., -10., -10., -10., -10., -5., -5., -5.,   0. };
  const Double_t dV0Max[nV0] = {  2.5,  1., 200., 200., 200., 200.,  5.,  5.,  5., 100. };
  THnSparseD *hsV0 = new THnSparseD("hsV0", "", nV0, nV0Bin, dV0Min, dV0Max); fOutputListMC->Add(hsV0);
//=============================================================================

  const Int_t nXi = 10; // 0: particle type
                        //    ==-0.5, XiNeg
                        //    == 0.5, XiPos
                        // 1: Event type
                        //    ==-0.5, INEL
                        //    == 0.5, MB
                        // 2: V0M
                        // 3: V0A
                        // 4: CL1
                        // 5: ZNA
                        // 6: eta
                        // 7: rap in Lab
                        // 8: rap in CMS
                        // 9: Pt
  const Int_t    nXiBin[nV0] = {  2,   2,  210,  210,  210,  210, 100, 100, 100, 1000  };
  const Double_t dXiMin[nV0] = { -1., -1., -10., -10., -10., -10., -5., -5., -5.,   0. };
  const Double_t dXiMax[nV0] = {  1.,  1., 200., 200., 200., 200.,  5.,  5.,  5., 100. };
  THnSparseD *hsXi = new THnSparseD("hsXi", "", nXi, nXiBin, dXiMin, dXiMax); fOutputListMC->Add(hsXi);

  TH1::AddDirectory(bStatusTmpH);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitAnalysis()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitAnalysis
//

  if (fCollisionType==(AliPicoHeaderCJ::kPP)) InitParamsPP();
  if (fCollisionType==(AliPicoHeaderCJ::kPA)) InitParamsPA();
  if (fCollisionType==(AliPicoHeaderCJ::kAP)) InitParamsAP();
  if (fCollisionType==(AliPicoHeaderCJ::kAA)) InitParamsAA();

  fAnaUtils = new AliAnalysisUtils();
  fAnaUtils->SetMinVtxContr(fCutMinEventVtxContr);
  fAnaUtils->SetMaxVtxZ(fCutMaxEventVzAbs);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsPP()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersPP
//

  fRapidityShift = 0.;

  fCutMaxEventVzAbs = 10.;
//=============================================================================

  fCutMaxV0Chi2   = 33.;
  fCutMinV0Radius = 0.3;   // default: >0.5; uncertainty: >0.3, 0.4, 0.6, 0.7;
  fCutMaxV0Radius = 200.;  // default: not applied in pp
//=============================================================================

  fCutMaxDausDCA                     = 1.5;  // default: <1;    uncertainty: <0.5, 0.75, 1.25, 1.5
  fCutMinDauDCAtoPV                  = 0.05; // default: >0.06; uncertainty: >0.05, 0.055, 0.07, 0.08
  fCutMinDauXrowsTPC                 = 70.;  // default: >70;   uncertainty: >75, 80
  fCutMinDauXrowsOverFindableClusTPC = 0.8;  // default: >0.8;  uncertainty: >0.95
//=============================================================================

  fCutMaxKshortSigmaTPC = -1.;   // default: <5;     uncertainty: w/o cut
  fCutMinKshortCosPA    = 0.95;  // default: >0.97;  uncertainty: >0.95, 0.96, 0.98, 0.99
  fCutMaxKshortCtau     = 30.;   // default: <20;    uncertainty: <12, 30
  fCutMaxKshortArmFrac  = -1.;   // default: not applied in pp
  fCutMinKshortDeltaM   = 0.003; // default: >0.005; uncertainty: >0.003, 0.006
//=============================================================================

  fCutMaxLambdaSigmaTPC  = 7.;    // default: <5;     uncertainty: 4, 6, 7
  fCutMinLambdaCosPA     = 0.993; // default: >0.995; uncertainty: >0.993, 0.994, 0.996, 0.997
  fCutMaxLambdaCtau      = 40.;   // default: <30;    uncertainty: <20, 40
  fCutMaxLambdaArmFrac   = -1.;   // default: not applied in pp
  fCutMinLambdaDeletaM   = -1.;   // default: >0.01;  uncertainty: w/o rejection

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsPA()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersPA
//

  InitParamsPP();

  fRapidityShift = 0.465;

  fCutMaxLambdaSigmaTPC = 6;  // default: <5; uncertaity: <4, 6

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsAP()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersAP
//

  InitParamsPP();

  fRapidityShift = -0.465;

  fCutMaxLambdaSigmaTPC = 6;  // default: <5; uncertaity: <4, 6

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsAA()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersAA
//

  InitParamsPP();

  fCutMinV0Radius = 0.9;  // default: 5; uncertainty: varying around
  fCutMaxV0Radius = 100.;

//fCutMinDauPt      = 0.16;  // default: 0.16
  fCutMinDauDCAtoPV = 0.08;  // default: 0.1; uncertainty: 0.08, 0.12

  fCutMaxKshortSigmaTPC = 6;      // default: <5;     uncertainty: <4, 6;
  fCutMinKshortCosPA    = 0.997;  // default: >0.998; uncertainty: 0.997, 0.999
  fCutMaxKshortArmFrac  = 0.2;    // default: <0.2

  fCutMaxLambdaSigmaTPC = 6;      // default: <5;     uncertaity: <4, 6
  fCutMinLambdaCosPA    = 0.997;  // default: >0.998; uncertainty: 0.997, 0.999
  fCutMinLambdaDeletaM  = 0.008;  // default: >0.01;  uncertainty: 0.008, 0.012

  return;
}
