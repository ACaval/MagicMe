#ifndef MAGICME_H
#define MAGICME_H

#include <TH2I.h>
#include <TMath.h>
#include <TString.h>
#include <Riostream.h>
#include <TRandom3.h>
#include <time.h>
#include <TGraph.h>

class MagicMe {
 public:
  MagicMe();
  void MakeMeFeel(Int_t step,Int_t iter);
  Double_t MaxTrack(Double_t prob,Int_t iter);

 protected:
  void ApplyProbability(Double_t prob);
  void CopyArray(Int_t **matrix, Int_t **m);
  Int_t FastCount(Int_t **matrix, Int_t value);
  TH2I * MakeHisto(Int_t **matrix, TString title);
  //Double_t MaxTrack(Double_t prob,Int_t iter);
  void SearchAndModify(Int_t x, Int_t y, Int_t **matrix);

  Int_t fCounter;
  const Int_t fkDimensions;
  Int_t **fLines;
  Int_t **fMatrix;
  Int_t fMatrixEntries;
  TH2I *fMatrixHisto;
  const Int_t fkMatrixSize;
  Int_t **fProbMatrix;
  Int_t fProbMatrixEntries;
  TH2I *fProbMatrixHisto;

  static const Int_t fgkDimensions;
};

#endif
