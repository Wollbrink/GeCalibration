/*
 * new.cpp
 *
 *  Created on: Jun 16, 2017
 *      Author: chrisw94
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TGraphErrors.h"

using namespace std;

vector<Double_t> fListOfEnergies;

class DetectorConfig672 {
public:
	double fHighVoltage, fCoarseGain, fFineGain, fShapeTime;

	bool mode; //if true, in triangle mode and if false, Gauss mode
	bool BLRrate; // if true, auto and if false, high
};

class RadioactiveSource {
public:
	string fName; //name of the source
	vector<double> fGammas; //collection of gamma energies from the source
	bool fVerbose;
	void AddGamma(double energy) {
		fGammas.push_back(energy);
		fListOfEnergies.push_back(energy);
		if (fVerbose) {
			cout << "The number of gammas: " << fGammas.size() << endl;
		}
	}
	void PrintStatus() {
		cout << "Source name is: " << fName << endl;

		vector<double>::iterator it;
		int counter = 1;
		for (it = fGammas.begin(); it != fGammas.end(); it++) {
			cout << "Energy of Gamma " << counter << " = " << *it << endl;
			counter++;
		}
	}
};

class HistogramAndADCPeaks {
public:
	vector<Double_t> fMeasuredADCError;

	TH1D* ConvertMaestroToHistogram(TString file, TString name) {
		ifstream fin(file.Data());
		char buffer[1024];
		if (!fin) {
			cout << " Oh dear could not open file " << endl;
		}
		int lineVals[7];
		int col0 = -1;
		int line = 0;
		vector<int> data;
		vector<int> channel;
		while (fin.getline(buffer, 1024)) {
			line++;
			if (line > 4) {
				for (int i = 0; i < 7; i++) {
					lineVals[i] = -1;
					sscanf(buffer, "%d:%d%d%d%d%d%d%d", &col0, &lineVals[0],
							&lineVals[1], &lineVals[2], &lineVals[3],
							&lineVals[4], &lineVals[5], &lineVals[6]);
					for (int i = 0; i < 7; i++) {
						channel.push_back(col0 + i);
						data.push_back(lineVals[i]);
					}
				}
			}
		}
		fin.close();
		int chanMin = channel[0];
		int chanMax = channel[channel.size() - 1];
		TH1D* hOut = new TH1D(name.Data(), "ADC Histogram", (chanMax - chanMin),
				chanMin, chanMax);

		for (unsigned int i = 0; i < channel.size(); i++) {
			hOut->Fill(hOut->FindBin(channel[i]), data[i]);
		}

		hOut->Sumw2(0); //Get rid of error bars for now
		return hOut;
	}

	TF1* CalibrationFunction() {
		TH1D* hOut = ConvertMaestroToHistogram(
				"/Users/chrisw94/Desktop/SubtractedBackgroundFromSources.Txt",
				"Run1Data");
		TSpectrum* s = new TSpectrum(11);
		Int_t nFound = s->Search(hOut, 2, "nobackground", 0.1);
		Double_t *fMeasuredADC = s->GetPositionX();

		//sorting fMeasuredADC

		cout << "The number of peaks that have been found: " << nFound << endl;
		cout
				<< "If less than expected, adjust 'a' in Search(-,-,-,a), where a is 0<a<1"
				<< endl;

		//defining crystalball

		//fit Gaussian and extract parameters
		TCanvas *can2 = new TCanvas("can2", "Histograms for Fits", 2100, 900);
		can2->Divide(2, 1);
		can2->cd(1);

		//create guess for standard deviations at peaks
		Double_t sigmaGuess[nFound];
		for (int i = 0; i < nFound; i++) {
			TF1 *a = new TF1("a", "gaus(0)", fMeasuredADC[i] - 15,
					fMeasuredADC[i] + 15);
			hOut->Fit("a", "Q", "", fMeasuredADC[i] - 15, fMeasuredADC[i] + 15);
			a = hOut->GetFunction("a");
			sigmaGuess[i] = a->GetParameter(2);
			delete a;
		}

		//perform fits with crystal ball function
		for (int i = 0; i < nFound; i++) {
			TF1 *b =
					new TF1("b",
							"(((x-[2])/[3])<=(-[0]))*([4])*(([1]/abs([0]))^[1])*(exp((-[0]^2)/2))*(([1]/abs([0]))-abs([0])-((x-[2])/[3]))^(-[1]) + \(((x-[2])/[3])>(-[0]))*([4])*(exp(-(x-[2])^2/(2*[3]^2)))",
							fMeasuredADC[i] - 20, fMeasuredADC[i] + 15);
			b->SetParameter(0, 1.5);
			b->SetParLimits(0,0,10);
			b->SetParameter(1, 2);
			b->SetParameter(2, fMeasuredADC[i]);
			b->SetParameter(3, sigmaGuess[i]);
			b->SetParameter(4, hOut->GetBinContent(fMeasuredADC[i]));
			hOut->Fit("b", "R+", "", fMeasuredADC[i] - 20,
					fMeasuredADC[i] + 15);
			cout << "NDF: " << b->GetNDF() << endl;
			b = hOut->GetFunction("b");
			b->DrawCopy("SAME");
			//need to extract parameters
			fMeasuredADC[i] = b->GetParameter(2);
			fMeasuredADCError.push_back(b->GetParameter(3));
			//delete b;

		}

		//sort the ADC's from low to high
		Double_t temp1;
		Double_t temp2;
		for (int i = 0; i < nFound; i++) {
			for (int j = i + 1; j < nFound; j++) {
				if (fMeasuredADC[i] > fMeasuredADC[j]) {
					temp1 = fMeasuredADC[i];
					temp2 = fMeasuredADCError.at(i);
					fMeasuredADC[i] = fMeasuredADC[j];
					fMeasuredADCError.at(i) = fMeasuredADCError.at(j);
					fMeasuredADC[j] = temp1;
					fMeasuredADCError.at(j) = temp2;
				}
			}
		}

		//sorting the energies from low to high
		Double_t temp3;
		for (int i = 0; i < fListOfEnergies.size() - 1; i++) {
			for (int j = i + 1; j < fListOfEnergies.size(); j++) {
				if (fListOfEnergies.at(i) > fListOfEnergies.at(j)) {
					temp3 = fListOfEnergies.at(i);
					fListOfEnergies.at(i) = fListOfEnergies.at(j);
					fListOfEnergies.at(j) = temp3;
				}

			}
		}
		can2->cd(2);
		double* x = fMeasuredADC;
		double* dx = &fMeasuredADCError[0];
		double* y = &fListOfEnergies[0];
		TGraphErrors* grE = new TGraphErrors(nFound, x, y, dx, 0);
		grE->Draw("AP");
		TF1* myFit1 = new TF1("myFit1", "pol1", 0, 6200);
		myFit1->SetParameters(0, 5);
		grE->Fit("myFit1", "R");
		//TF1* myFit2 = new TF1("myFit2", "pol2", 0, 6200);
		//myFit2->SetParameters(0, 5);
		//grE->Fit("myFit2", "R");
		//TF1* myFit3 = new TF1("myFit3", "pol3", 0, 6200);
		//myFit3->SetParameters(0, 5);
		//grE->Fit("myFit3", "R");
		can2->SaveAs("fits.pdf");
		return myFit1;	//change depending on the degree of the polynomial fit
	}
};

int main() {

	RadioactiveSource source1; //Declaring source object
	source1.fName = "Cesium137"; //Naming the object
	source1.AddGamma(661.7); /*Adding the Energies associated with the source.
	 *If the source emits a low energy then the search function will not be able to detect
	 *the peak, so make sure the peaks are above the Compton edges
	 */

	RadioactiveSource source2;
	source2.fName = "Cobalt60";
	source2.AddGamma(1173.2);
	source2.AddGamma(1332.5);

	RadioactiveSource source3;
	source3.fName = "Barium133";
	source3.AddGamma(81.0);
	source3.AddGamma(276.4);
	source3.AddGamma(302.9);
	source3.AddGamma(356.0);
	source3.AddGamma(383.8);

	RadioactiveSource source4;
	source4.fName = "Sodium22";
	source4.AddGamma(511.0);

	HistogramAndADCPeaks hOut;

	TCanvas* can = new TCanvas("can", "ADC vs Counts", 2100, 900); //creating a canvas for ADC plot
	TH1D* adcHisto = hOut.ConvertMaestroToHistogram(
			"/Users/chrisw94/Desktop/SubtractedBackgroundFromSources.Txt",
			"ADC");
	can->cd(1);
	adcHisto->Draw();
	hOut.CalibrationFunction(); //making the calibration function, where p0 is the intercept
	//p1 is the slope of the line
	can->SaveAs("plot.pdf");

	return 0;
}
