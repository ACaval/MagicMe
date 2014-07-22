#include <TCanvas.h>
#include <TFile.h>
#include "MagicMe.h"


const Int_t MagicMe::fgkDimensions = 4;

//________________________________________________________
MagicMe::MagicMe():
  fCounter(0),
  fkDimensions(fgkDimensions),
  fkMatrixSize(TMath::Power(4,fgkDimensions))
{
  //Creazione degli array iniziali
  fLines = new Int_t*[fkMatrixSize];
  fMatrix = new Int_t*[fkMatrixSize];
  fProbMatrix = new Int_t*[fkMatrixSize];
  for(Int_t i=0; i<fkMatrixSize;++i) {
    fLines[i] = new Int_t[fkDimensions];
    fMatrix[i] = new Int_t[fkMatrixSize];
    fProbMatrix[i] = new Int_t[fkMatrixSize];
    for(Int_t k=0; k<fkDimensions;++k) {
      fLines[i][k]=0;
    }
    for(Int_t j=0; j<fkMatrixSize;++j) {
      fMatrix[i][j]=0;
      fProbMatrix[i][j]=0;
      }
  }
  
  //Riempimento dell'array di righe e colonne
  Int_t count[fkDimensions];
  for(Int_t i=0; i<fkDimensions; ++i) count[i]=0;
  
  for(Int_t k=0; k<fkMatrixSize; ++k) {
    for(Int_t j=fkDimensions-1;j>=0;--j) {
      if(count[j]==4) {
	if(j!=0) count[j-1]++;
	count[j]=0;
      }
    }
    for(Int_t i=0;i<fkDimensions;++i) {
      fLines[k][i]=count[i];
      //cout << count[i];
    }
    //cout << endl;
    ++count[fkDimensions-1];
  }

  //Riempimento delle matrici
  for(Int_t r=0; r<fkMatrixSize;++r) {
    for(Int_t c=fkMatrixSize-1; c>=r;--c) {
      Int_t equality=0;
      for(Int_t n=0; n<fkDimensions;++n) {
	if(fLines[r][n]==fLines[c][n]) ++equality;
      }
      if(equality==fkDimensions-1) {
	fMatrix[r][c]=1;
	fMatrix[c][r]=1;
	fProbMatrix[r][c]=1;
	fProbMatrix[c][r]=1;
      }
    }
  }
  fMatrixHisto=MakeHisto(fMatrix,"Link statici");
  //fMatrixHisto->Draw("colz");
  fMatrixEntries=FastCount(fMatrix,1);
}

//________________________________________________________
void MagicMe::MakeMeFeel(Int_t step,Int_t iter) {
  Double_t prob[step];
  Double_t stat[step];
  for(Int_t i=0;i<step;i++) {
    prob[i]=0;
    stat[i]=0;
  }

  for(Int_t k=0;k<step;k++) { 
    prob[k]=((k+1)/(Double_t)step);
    printf("Step #%i, probability %f. Now iterating...\n",k,prob[k]);
    stat[k]+=((Double_t)MaxTrack(prob[k],iter)/(Double_t)fMatrixEntries);
    //cout << prob[k] << "\t" << stat[k] << "\t" << MaxTrack(prob[k],iter) << endl;
  }
  TGraph * Res = new TGraph(step,prob,stat);
  TCanvas * c1 = new TCanvas(Form("Dim%i_Step%i_Iter_%i",fkDimensions,step,iter),Form("Dim%i_Step%i_Iter_%i",fkDimensions,step,iter));
  c1->cd();
  Res->SetMarkerStyle(22);
  Res->SetMarkerSize(1.0);
  //Res->Fit("pol2");
  Res->Draw("AP");

  TFile resultFile("result.root","update");
  c1->Write();
  resultFile.Close();
  return;
}

//________________________________________________________
Double_t MagicMe::MaxTrack(Double_t prob,Int_t iter) {
  Double_t retValue=0;
  for(Int_t i=0; i<iter;++i) {
    //printf("\tIteration #%i\n",i+1);
    ApplyProbability(prob);
    fCounter=0;
    //fProbMatrixHisto=MakeHisto(fProbMatrix,Form("Variazione %f",prob));
    Int_t max=0;
    for(Int_t r=0; r<fkMatrixSize; ++r) {
      for(Int_t c=0; c<fkMatrixSize; ++c) {
	if(fProbMatrix[r][c]==1) {
	  ++fCounter;
	  Int_t **m;
	  m=new Int_t*[fkMatrixSize];
	  for(Int_t count=0;count<fkMatrixSize;++count) m[count]=new Int_t[fkMatrixSize];
	  CopyArray(fProbMatrix,m);
	  SearchAndModify(r,c,m);
	  Int_t val = FastCount(m,2);
	  //printf("(%i,%i) -> %i\n",r,c,fCounter);
	  for(Int_t count=0;count<fkMatrixSize;++count) delete [] m[count];
	  delete [] m;
	  if(max<val) max=val;
	}
      }
    }
    retValue+=max/(Double_t)iter;
  }
  //printf("%f->%i\n",prob,fCounter);
  return retValue;
  //cout << prob << "\t" << max << "\t" << max/(Double_t)fMatrixEntries<<endl;
}

//________________________________________________________
void MagicMe::ApplyProbability(Double_t prob) {
  TRandom3 gen(time(NULL));
  CopyArray(fMatrix,fProbMatrix);
  for(Int_t r=0; r<fkMatrixSize; ++r) {
    for(Int_t c=fkMatrixSize-1; c>=r; --c) {
      if(gen.Uniform(1)>prob) {
	fProbMatrix[r][c]=0;
	fProbMatrix[c][r]=0;
      } 
    } 
  }
  fProbMatrixEntries=FastCount(fProbMatrix,1);
  //printf("Prob. %f) Matrix Entries %i\n",prob,fProbMatrixEntries);
}

//________________________________________________________
void MagicMe::CopyArray(Int_t **matrix, Int_t **m) {
  for(Int_t r=0; r<fkMatrixSize;++r) {
    for(Int_t c=0; c<fkMatrixSize;++c) {
      m[r][c]=matrix[r][c];
    }
  }
}

//________________________________________________________
TH2I * MagicMe::MakeHisto(Int_t ** matrix,TString title) {
  TH2I * histo = new TH2I(title,title,fkMatrixSize,0,fkMatrixSize,fkMatrixSize,0,fkMatrixSize);
  for(Int_t r=0; r<fkMatrixSize; ++r) {
    for(Int_t c=0; c<fkMatrixSize; ++c) {
      if(matrix[r][c]!=0) histo->Fill(r,c,matrix[r][c]);
    }
  }
  return histo;
}

//________________________________________________________
void  MagicMe::SearchAndModify(Int_t x,Int_t y,Int_t **matrix) {
  matrix[x][y]=2;
  matrix[y][x]=2;
  for(Int_t onX=0;onX<fkMatrixSize;++onX) {
    if(matrix[onX][y]!=1) continue; 
    SearchAndModify(onX,y,matrix);
    return;
  }

  for(Int_t onY=0;onY<fkMatrixSize;++onY) {
    if(matrix[x][onY]!=1) continue; 
    SearchAndModify(x,onY,matrix);
    return;
  }
}

//________________________________________________________
Int_t MagicMe::FastCount(Int_t **matrix,Int_t value) {
  Int_t count=0;
  for(Int_t r=0; r<fkMatrixSize; ++r) {
    for(Int_t c=0; c<fkMatrixSize; ++c) {
      if(matrix[r][c]==value) ++count;
    }
  }
  return count;
}
