// simple example of using ROOT libraries in a C++ program with graphics
// and use of TASImage class

#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"


#include "assert.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <random>
#include <csignal>
#include <cstdio>
#include <cstdlib>

using namespace std;

//random number generator device
mt19937 gen(random_device{}());
int randInt(int a, int b){
	uniform_int_distribution <> d(a,b);
	return d(gen);
}


double randDouble(double a, double b){
	uniform_real_distribution <> d(a,b);
	return d(gen);
}

//function to caluclate the relative colormetric distance between two pixels
double colometricDistance(UInt_t pixel1, UInt_t pixel2){
	
	double alpha1 = pixel1 & 0xff000000;
	double red1 = pixel1 & 0x00ff0000;
	double green1 = pixel1 & 0x0000ff00;
	double blue1 = pixel1 & 0x000000ff;

	double alpha2 = pixel2 & 0xff000000;
	double red2 = pixel2 & 0x00ff0000;
	double green2 = pixel2 & 0x0000ff00;
	double blue2 = pixel2 & 0x000000ff;

	double delAlpha = alpha2 - alpha1;
	double delRed = red2-red1;
	double delGreen = green2-green1; 
	double delBlue = blue2-blue1;


	double r = sqrt( (delAlpha*delAlpha) + (delRed*delRed) + (delGreen*delGreen) +  (delBlue * delBlue) );

	return r;
}

//function to calculate total concatenated distance of each pixel pair
double totalColometricDistance(UInt_t * tgtPix, UInt_t * srcPix, Long_t numPix){
	double totalDistance = 0.0;
	for (int i = 0; i<numPix; i++){
		double pairDistance = colometricDistance(srcPix[i], tgtPix[i]);
		totalDistance += pairDistance; 
	}

	return totalDistance;
}

void melt(UInt_t * tgtPix, UInt_t * srcPix, Long_t numPix, double T0, double meltingIterations){
	double oldDist = totalColometricDistance(tgtPix, srcPix, numPix);
	double newDist = 0.0;
	for (int i = 0; i<numPix; i++){
		int randPixel1 = randInt(0, numPix-1);
		int randPixel2 = randInt(0, numPix-1);
		while (randPixel2 == randPixel1){
			randPixel2 = randInt(0, numPix-1);
		}

		double initialPairDistancePixel1 = colometricDistance(srcPix[randPixel1], tgtPix[randPixel1]);
		double initialPairDistancePixel2 = colometricDistance(srcPix[randPixel2], tgtPix[randPixel2]);
		double totalInitialPairDistance = initialPairDistancePixel1 + initialPairDistancePixel2;
		
		swap(srcPix[randPixel1], srcPix[randPixel2]);

		double finalPairDistancePixel1 = colometricDistance(srcPix[randPixel1], tgtPix[randPixel1]);
		double finalPairDistancePixel2 = colometricDistance(srcPix[randPixel2], tgtPix[randPixel2]);
		double totalFinalPairDistance = finalPairDistancePixel1 + finalPairDistancePixel2;
	
		newDist = oldDist-totalInitialPairDistance+totalFinalPairDistance;
		
		double deltaDist = newDist - oldDist;

		if (deltaDist < 0){
			oldDist = newDist;
			continue;
		}
		else {
			double p = exp(-deltaDist/T0);
			double r = randDouble(0.0, 1.0);

			if (r<p){
				oldDist = newDist;
				continue;
			}
			swap(srcPix[randPixel1], srcPix[randPixel2]);
		}		

	}
	printf("Finished melting\n"); 

}

void simulatedAnnealingPixelSwap(UInt_t * tgtPix, UInt_t * srcPix, Long_t numPix, double T0, double iterationsPerTemperature){
	printf("Running pixel swap formula\n"); 
	double T = T0;
	double oldDist = totalColometricDistance(tgtPix, srcPix, numPix);
	double newDist = 0.0;

	while(T>0){
		for (int i = 0; i < iterationsPerTemperature; i++){
			int randPixel1 = randInt(0, numPix);
			int randPixel2 = randInt(0, numPix);
			double initialPairDistancePixel1 = colometricDistance(srcPix[randPixel1], tgtPix[randPixel1]);
			double initialPairDistancePixel2 = colometricDistance(srcPix[randPixel2], tgtPix[randPixel2]);
			double totalInitialPairDistance = initialPairDistancePixel1 + initialPairDistancePixel2;
			
			swap(srcPix[randPixel1], srcPix[randPixel2]);
			
			double finalPairDistancePixel1 = colometricDistance(srcPix[randPixel1], tgtPix[randPixel1]);
			double finalPairDistancePixel2 = colometricDistance(srcPix[randPixel2], tgtPix[randPixel2]);
			double totalFinalPairDistance = finalPairDistancePixel1 + finalPairDistancePixel2;
		
			newDist = oldDist-totalInitialPairDistance+totalFinalPairDistance;

			double deltaDist = newDist - oldDist;
			
			if (deltaDist < 0){
				oldDist = newDist;
				continue;
			}
			else {
				double p = exp(-deltaDist/T);
				double r = randDouble(0.0,1.0);

				if (r < p){
					oldDist = newDist;
					continue;
				}
				swap(srcPix[randPixel1], srcPix[randPixel2]);
			}
		}
		T-=.1; 	
	}
}


int main(int argc, char **argv){

  clock_t tStart = clock();
  if (argc<3) {
    cout << "Usage: simapix_start image1 image2 <output=out.png>" << endl;
    return 0; 
  }
  TString fsrc=argv[1];
  TString ftgt=argv[2];
  TString fout;
  argc>3 ? fout = argv[3] : fout="out.png";
  cout << "Reading images: source= " << fsrc << " target= " << ftgt << endl;
  cout << "Output= " << fout << endl;

  TApplication theApp("App", &argc, argv);

  // create image objects
  TASImage *src = new TASImage(fsrc.Data());
  TASImage *tgt = new TASImage(ftgt.Data());
  TASImage *out = new TASImage(*src); // start with copy of source

  // Test image geometry, exit if they are not the same dimensions
  assert ( src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight() );
  cout << "Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;
  Long_t numPix=src->GetWidth()*src->GetHeight();

  // *** The work happens here
  // access the pixels for the output image 
  // each pixel is a 32-bit word, 1 byte each for (alpha,red,green,blue)
  // don't touch alpha (bits 31:28)
  UInt_t *tgtPix = tgt->GetArgbArray();  
  UInt_t *srcPix = out->GetArgbArray();
  
  double T0 = 10000;
  double meltingIterations = 100000;
  double iterationsPerTemperature = 1000;

  melt(tgtPix, srcPix, numPix, T0, meltingIterations);
  simulatedAnnealingPixelSwap(tgtPix, srcPix, numPix, T0, iterationsPerTemperature);

  printf("Solution execution time: %.2fs\n", (double) (clock()-tStart)/CLOCKS_PER_SEC);
  // *************************


  // print the results
  TCanvas *c1 = new TCanvas("c1", "images", 640, 480);
  c1->Divide(2,2);

  c1->cd(1);
  c1->Draw();
  src->Draw("X");
  c1->cd(2);
  tgt->Draw("X");
  c1->cd(3);
  out->Draw("X");
  c1->Print("scottToRotunda.png");
  
  // save the new image
  out->WriteImage(fout.Data());

  // coment out the lines for running in batch mode
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();

}
