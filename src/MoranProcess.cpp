#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
//using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]

int chooseCell(NumericVector nCells, double r){
// nCells have two values, nCells[0] is number of WT cells, and nCells[1] is number of MT cells
// r > 1 if we want to choose division cell
// r = 1 if we want to select cell to die
// WT cell is chosen if 0 returned
// MT cell is chosen if 1 returned

	//Get total fitness of WT and MT
	double fWT = nCells[0];
	double fMT = r*nCells[1];
	double fCut= fWT/(fWT + fMT);

	//random number generator
	NumericVector ru = runif(1);
	if(ru[0] <= fCut)
		return(0);
	else
		return(1);
}

// [[Rcpp::depends(RcppArmadillo)]]

int mutateCell(double u){
// u is mutation rate
// Cell will be mutated if 1 is returned
// Cell will not be mutated if 0 is returned

	//random number generator
	NumericVector ru = runif(1);
	if(ru[0] <= u)
		return(1);
	else
		return(0);
}

// [[Rcpp::depends(RcppArmadillo)]]

double mutantDection(NumericVector nCells, double tau, double alpha){
// nCells have two values, nCells[0] is number of WT cells, and nCells[1] is number of MT cells
// N*q = alpha = detection probablity

	//get detection probablity
	double Nc=  (nCells[0] + nCells[1]);
	double dP = nCells[1]*(alpha/Nc)*(tau/Nc);
	
	//random number generator
	NumericVector ru = runif(1);
	if(dP > ru[0])
		return(1);
	else
		return(0);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MoranPorcess(double N, double Age, double tau, double r, double u, double alpha){

	// initialized nCells
	NumericVector nCells(2);
	nCells[0] = N;
	nCells[1] = 0;
	int dP    = 0;
	double nDivision = Age*365/tau;

	for(int i=0; i<nDivision; i++){

		for(int j=0; j<N; j++){

			// choose one cell to divide
			int flagBirth = chooseCell(nCells, r);
			int flagMutant = 0;

			if(flagBirth == 1)
				nCells[1] = nCells[1] + 1;
			else {
				flagMutant = mutateCell(u);
				if(flagMutant == 0)
					nCells[0] = nCells[0] + 1;
				else
					nCells[1] = nCells[1] + 1;
			}

			// choose one cell to die
			int flagDeath = chooseCell(nCells, 1);
			if(flagDeath == 0)
				nCells[0] = nCells[0] - 1;
			else
				nCells[1] = nCells[1] - 1;

			// detection
			dP = dP + mutantDection(nCells, tau, alpha);
		}
	}

	// return list
	return List::create(nCells, dP);
}
