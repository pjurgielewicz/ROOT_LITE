#include "TH2PolyBinLite.h"

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

//TH2PolyBinLite::TH2PolyBinLite()
//{
//	fPoly = 0;
//	fContent = 0.;
//	fNumber = 0;
//	fXmax = -1111;
//	fXmin = -1111;
//	fYmax = -1111;
//	fYmin = -1111;
//	fArea = 0;
//	SetChanged(kTRUE);
//}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TH2PolyBinLite::TH2PolyBinLite(Int_t nVerts, const Double_t* x, const Double32_t* y, Int_t bin_number) : nVerts(nVerts), x(x), y(y)
{
	fContent = 0.;
	fNumber = bin_number;
	fArea = 0.;
	fXmax = -1111;
	fXmin = -1111;
	fYmax = -1111;
	fYmin = -1111;
	SetChanged(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2PolyBinLite::~TH2PolyBinLite()
{
	// ???
	if (x) delete[] x;
	if (y) delete[] y;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the area of the bin.

Double_t TH2PolyBinLite::GetDet()
{
	if (nVerts != 3) return 0;

	Double_t v1x = x[1] - x[0];
	Double_t v2x = x[2] - x[0];
	Double_t v1y = y[1] - y[0];
	Double_t v2y = y[2] - y[0];

    return v1x * v2y - v1y * v2x;
}

Double_t TH2PolyBinLite::GetArea()
{
	if (nVerts == 2)
	{
        return std::fabs((x[1] - x[0]) * (y[1] - y[0]));
	}
	else if (nVerts == 3)
	{
        return std::fabs(GetDet()) * 0.5;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum value for the x coordinates of the bin.

#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)

Double_t TH2PolyBinLite::GetXMax()
{
	if (fXmax != -1111) return fXmax;

	fXmax = (nVerts == 2 ? (MAX(x[0], x[1])) : (MAX(MAX(x[0], x[1]), x[2])));

	return fXmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the minimum value for the x coordinates of the bin.

Double_t TH2PolyBinLite::GetXMin()
{
	if (fXmin != -1111) return fXmin;

	fXmax = (nVerts == 2 ? (MIN(x[0], x[1])) : (MIN(MIN(x[0], x[1]), x[2])));

	return fXmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum value for the y coordinates of the bin.

Double_t TH2PolyBinLite::GetYMax()
{
	if (fYmax != -1111) return fYmax;

	fYmax = (nVerts == 2 ? (MAX(y[0], y[1])) : (MAX(MAX(y[0], y[1]), y[2])));

	return fYmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the minimum value for the y coordinates of the bin.

Double_t TH2PolyBinLite::GetYMin()
{
	if (fYmin != -1111) return fYmin;

	fYmin = (nVerts == 2 ? (MIN(y[0], y[1])) : (MIN(MIN(y[0], y[1]), y[2])));

	return fYmin;
}

#undef MIN
#undef MAX

////////////////////////////////////////////////////////////////////////////////
/// Return the set of vertices to build TGraph object

Int_t TH2PolyBinLite::BuildFullPolyDescription(Double_t *x, Double_t *y)
{
    if (fNumber == 2) { // RECTANGLE CASE
        x[0] = this->x[0]; y[0] = this->y[0];
        x[1] = this->x[0]; y[1] = this->y[1];
        x[2] = this->x[1]; y[2] = this->y[1];
        x[3] = this->x[1]; y[3] = this->y[0];
        x[4] = this->x[0]; y[4] = this->y[0];

        return fNumber;
    }
    else if (fNumber == 3) { // TRIANGLE CASE
        x[0] = this->x[0]; y[0] = this->y[0];
        x[1] = this->x[1]; y[1] = this->y[1];
        x[2] = this->x[2]; y[2] = this->y[2];
        x[4] = this->x[0]; y[4] = this->y[0];

        return fNumber;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return "true" if the point (x,y) is inside the bin.

Bool_t TH2PolyBinLite::IsInside(Double_t x, Double_t y)
{
	Int_t in = 0;

	if (nVerts == 2)
	{
		if (GetXMin() <= x && x <= GetXMax() &&
			GetYMin() <= y && y <= GetYMax())
			in = 1;
	}
	else if (nVerts == 3)
	{
		Double_t detInv = 1. / GetDet();

        Double_t l1 = ((this->y[1] - this->y[2]) * (x - this->x[2]) +
                       (this->x[2] - this->x[1]) * (y - this->y[2])) * detInv;

        Double_t l2 = ((this->y[2] - this->y[0]) * (x - this->x[2]) +
                       (this->x[0] - this->x[2]) * (y - this->y[2])) * detInv;

		in = (l1 > 0.0 && l2 > 0.0 && l1 + l2 < 1.0) ? 1 : 0;
	}

	return in;
}
