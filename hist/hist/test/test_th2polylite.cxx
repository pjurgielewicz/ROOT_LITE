#include "TClass.h"
#include "TList.h"
#include "TRandom3.h"
#include "TH2Poly.h"
#include "TH2PolyLite.h"

//#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>

using namespace std;

#define ASSERT_NEAR(x1, x2, delta) (assert(std::fabs(x1 - x2) < delta))

Float_t delta = 0.00000000001;

void FillForTest(TH2Poly* th2poly, TH2PolyLite* th2polyLite, TRandom& ran) {
	Double_t value, weight;
	Double_t px, py;

	for (Int_t i = 0; i < 1000000; i++) {
		px = ran.Gaus(5, 2);
		py = ran.Gaus(4, 2);
		value = ran.Gaus(20, 5);
		//value = ran.Uniform(0, 20);
		weight = ran.Gaus(17, 20);

		th2poly->Fill(px, py, value);
		th2polyLite->Fill(px, py, value);
	}
}

void BinContentCompare(TH2Poly* th2poly, TH2PolyLite* th2polyLite) {
	Double_t cont1, cont2;
	for (Double_t y = 0.5; y<10; y += 2.0) {
		for (Double_t x = 0.5; x<10; x += 2.0) {
			cont1 = th2poly->GetBinContent(th2poly->FindBin(x, y));
			cont2 = th2polyLite->GetBinContent(th2polyLite->FindBin(x, y));
			ASSERT_NEAR(cont1, cont2, delta);
		}
	}
	// test overflow
	cont1 = th2poly->GetBinContent(th2poly->FindBin(11, 11));
	cont2 = th2polyLite->GetBinContent(th2polyLite->FindBin(11, 11));
	ASSERT_NEAR(cont1, cont2, delta);
}

void CreateGrid(TH2Poly* th2poly, TH2PolyLite* th2polyLite)
{
	double minx = -10; double maxx = 10;
	double miny = -10; double maxy = 10;
	double binsz = 2;

	for (double i = minx; i < maxx; i += binsz) {
		for (double j = miny; j < maxy; j += binsz) {
			th2poly->AddBin(i, j, i + binsz, j + binsz);
			th2polyLite->AddBin(i, j, i + binsz, j + binsz);
		}
	}
}

void test_globalStats() {

	TH2Poly* th2poly = new TH2Poly("TH2Poly", "TH2Poly", 10, -10, 10, 10, -10, 10);
	TH2PolyLite* th2polyLite = new TH2PolyLite("TH2PolyLite", "TH2PolyLite", 10, -10, 10, 10, -10, 10);

	TRandom3 ran(1);

	// ABSOLUTE BASICS

	CreateGrid(th2poly, th2polyLite);
	FillForTest(th2poly, th2polyLite);
	BinContentCompare(th2poly, th2polyLite);

	// ADD/MERGE
	// ...

	delete th2poly;
	delete th2polyLite;
}

// ------------ TEST CALLS ------------

//TEST(TProfile2Poly, GlobalCompare) {
//	test_globalStats();
//}

void testMain()
{
	test_globalStats();
}