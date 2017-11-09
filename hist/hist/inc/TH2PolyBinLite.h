#ifndef ROOT_TH2PolyBinLite
#define ROOT_TH2PolyBinLite

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TH2Poly                                                              //
//                                                                      //
// 2-Dim histogram with polygon bins                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TH2
#include "TH2.h"
#endif

/** \class TH2PolyBinLite
\ingroup Hist
Helper class to represent a bin in the TH2PolyLite histogram
*/

class TH2PolyBinLite {
	// LITE VERSION OF THE CLASS THAT CAN HANDLE TRIANGLES AND RECTANGLES ONLY
	// BINS ARE ANONYMOUS (DO NOT HAVE TITLE, NAME USW.) JUST ITS ID (INT)
	TH2PolyBinLite(Int_t nVerts, Double_t* x, Double32_t* y, Int_t bin_number);
	virtual ~TH2PolyBinLite();

	void      ClearContent() { fContent = 0; }
	void      Fill(Double_t w) { fContent = fContent + w; SetChanged(true); }

	Double_t  GetArea();

	Double_t  GetContent() const { return fContent; }
	Bool_t    GetChanged() const { return fChanged; }
	Int_t     GetBinNumber() const { return fNumber; }

	Double_t  GetXMax();
	Double_t  GetXMin();
	Double_t  GetYMax();
	Double_t  GetYMin();

	const Double_t* GetXList() const { return x; }
	const Double_t* GetYList() const { return y; }
	Int_t     GetNVerts() const { return nVerts; }

    Int_t     BuildFullPolyDescription(Double_t *x, Double_t *y);

	Bool_t    IsInside(Double_t x, Double_t y) const;
	void      SetChanged(Bool_t flag) { fChanged = flag; }
	void      SetContent(Double_t content) { fContent = content; SetChanged(true); }

protected:
	Bool_t    fChanged;   //For the 3D Painter
	Int_t     fNumber;    //Bin number of the bin in TH2Poly
						  //TObject  *fPoly;      //Object holding the polygon definition
	Double_t* x;
	Double_t* y;
	Int_t nVerts;

	Double_t  fArea;      //Bin area
	Double_t  fContent;   //Bin content
	Double_t  fXmin;      //X minimum value
	Double_t  fYmin;      //Y minimum value
	Double_t  fXmax;      //X maximum value
	Double_t  fYmax;      //Y maximum value

private:
	Double_t GetDet();

	ClassDef(TH2PolyBinLite, 1)  //2-Dim polygon bins
};
