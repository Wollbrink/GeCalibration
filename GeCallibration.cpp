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

	TF1* CalibrationFunction(TString file, TString name) {
		TH1D* hOut = ConvertMaestroToHistogram(file, name);
		TSpectrum* s = new TSpectrum(9);
		Int_t nFound = s->Search(hOut, 2, "nobackground", 0.05);
		Double_t *fMeasuredADC = s->GetPositionX();

		//sorting fMeasuredADC

		cout << "The number of peaks that have been found: " << nFound << endl;
		cout
				<< "If less than expected, adjust 'a' in Search(-,-,-,a), where a is 0<a<1"
				<< endl;

		//fit Gaussian and extract parameters
		TCanvas *can2 = new TCanvas("can2", "Histograms for Fits", 2100, 900);
		can2->Divide(2, 1);
		can2->cd(1);
		//create guess for standard deviations at peaks
		Double_t sigmaGuess[nFound];
		Double_t halfPeakPosition[nFound];
		Double_t peakHeight[nFound];
		for (int i = 0; i < nFound; i++) {
			int cont = fMeasuredADC[i] - 50; // Continuum location. Adjust according to peak spacing
			peakHeight[i] = (hOut->GetBinContent(fMeasuredADC[i]))
					- (hOut->GetBinContent(cont));
			for (int j = 1; j < 50; j--) {
				halfPeakPosition[i] = fMeasuredADC[i] - j;
				if ((peakHeight[i] / 2)
						> (hOut->GetBinContent(halfPeakPosition[i])
								- hOut->GetBinContent(cont))) {
					halfPeakPosition[i] += 1;
					break;
				}
			}
			sigmaGuess[i] = (fMeasuredADC[i] - halfPeakPosition[i])
					/ (sqrt(2 * log(2)));
			TF1 *a = new TF1("a", "gaus(0)",
					fMeasuredADC[i] - (3 * sigmaGuess[i]),
					fMeasuredADC[i] + (3 * sigmaGuess[i]));
			a->SetParameter(2, sigmaGuess[i]);
			hOut->Fit("a", "RQ");
			a = hOut->GetFunction("a");
			//a->DrawCopy("SAME");
			//cout << "Chi2/NDF:" << (a->GetChisquare()) / (a->GetNDF()) << endl;
			sigmaGuess[i] = a->GetParameter(2);
			delete a;
		}
		//Find guesses for linear fit on left of peak
		Double_t fLeftSlope[nFound];
		Double_t fLeftIntercept[nFound];
		for (int i = 0; i < nFound; i++) {
			TF1* left = new TF1("left", "pol1",
					fMeasuredADC[i] - 6 * sigmaGuess[i],
					fMeasuredADC[i] - 4 * sigmaGuess[i]);
			hOut->Fit("left", "RQ");
			fLeftIntercept[i] = left->GetParameter(0);
			fLeftSlope[i] = left->GetParameter(1);
			delete left;
		}

		//Find guesses for linear fit on right of peak
		Double_t fRightSlope[nFound];
		Double_t fRightIntercept[nFound];
		for (int i = 0; i < nFound; i++) {
			TF1* right = new TF1("right", "pol1",
					fMeasuredADC[i] + 4 * sigmaGuess[i],
					fMeasuredADC[i] + 6 * sigmaGuess[i]);
			hOut->Fit("right", "RQ");
			fRightIntercept[i] = right->GetParameter(0);
			fRightSlope[i] = right->GetParameter(1);
			delete right;
		}

		//perform fits with gauss-linear function

		for (int i = 0; i < nFound; i++) {
			TF1 *b = new TF1("b", "pol1(0) + gaus(2) + pol1(5)",
					fMeasuredADC[i] - 6 * sigmaGuess[i],
					fMeasuredADC[i] + 6 * sigmaGuess[i]);
			b->SetParameter(0, fLeftIntercept[i]);
			b->SetParameter(1, fLeftSlope[i]);
			b->SetParameter(2, hOut->GetBinContent(fMeasuredADC[i]));
			b->SetParameter(3, fMeasuredADC[i]);
			b->SetParameter(4, sigmaGuess[i]);
			b->SetParameter(5, fRightIntercept[i]);
			b->SetParameter(6, fRightSlope[i]);
			hOut->Fit("b", "R");
			cout << "Chi2/NDF: " << (b->GetChisquare()) / (b->GetNDF()) << endl;
			b = hOut->GetFunction("b");
			b->DrawCopy("SAME");
			//need to extract parameters
			fMeasuredADC[i] = b->GetParameter(3);
			fMeasuredADCError.push_back(b->GetParameter(4));
			delete b;

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
		Double_t* x = fMeasuredADC;
		Double_t* dx = &fMeasuredADCError[0];
		Double_t* y = &fListOfEnergies[0];
		TGraphErrors* grE = new TGraphErrors(nFound, x, y, dx, 0);
		grE->Draw("AP");
		//TF1* myFit1 = new TF1("myFit1", "pol1", 0, 6200);
		//myFit1->SetParameters(0, 5);
		//grE->Fit("myFit1", "R");
		//cout << "Chi2/NDF: " << (myFit1->GetChisquare())/(myFit1->GetNDF()) << endl;
		//TF1* myFit2 = new TF1("myFit2", "pol2", 0, 6200);
		//myFit2->SetParameters(0, 5);
		//grE->Fit("myFit2", "R");
		//cout << "Chi2/NDF: " << (myFit2->GetChisquare())/(myFit2->GetNDF()) << endl;
		TF1* myFit3 = new TF1("myFit3", "pol1", 0, 6200);
		myFit3->SetParameters(0, 5);
		grE->Fit("myFit3", "R");
		cout << "Chi2/NDF: " << (myFit3->GetChisquare()) / (myFit3->GetNDF())
				<< endl;
		can2->SaveAs("fits.pdf");
		return myFit3;	//change depending on the degree of the polynomial fit
	}

	TH1D* energyCalibratedSpectrum(TF1* fun, TH1D* adcHisto) {
		Int_t nbins = adcHisto->GetXaxis()->GetNbins();
		Double_t adcMax = adcHisto->GetXaxis()->GetXmax();
		Double_t adcMin = adcHisto->GetXaxis()->GetXmin();

		Double_t eMin = fun->Eval(adcMin);
		Double_t eMax = fun->Eval(adcMax);

		TString name = Form("Energy_%s", adcHisto->GetName());
		TString title = Form("Energy_%s", adcHisto->GetTitle());

		TH1D* eSpectrum = new TH1D(name, title, nbins, eMin, eMax);

		for (int i = 0; i < nbins; i++) {
			Double_t adc = adcHisto->GetBinCenter(i);
			Double_t e = fun->Eval(adc);
			Double_t counts = adcHisto->GetBinContent(i);
			eSpectrum->Fill(e, counts);
			eSpectrum->Sumw2(0);
		}

		return eSpectrum;
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
	can->Divide(2, 1);
	TH1D* adcHisto = hOut.ConvertMaestroToHistogram(
			"/Users/chrisw94/Desktop/20170728_BackgroundSubtractedSources.Txt",
			"ADC");
	can->cd(1);
	adcHisto->Draw();
	TF1* fun = (TF1*) hOut.CalibrationFunction(
			"/Users/chrisw94/Desktop/20170728_BackgroundSubtractedSources.Txt",
			"Calibration"); //making the calibration function, where p0 is the intercept
	//p1 is the slope of the line
	can->cd(2);
	TH1D* eHisto = hOut.energyCalibratedSpectrum(fun, adcHisto);
	eHisto->Draw();
	can->SaveAs("plot.pdf");

	return 0;
}
