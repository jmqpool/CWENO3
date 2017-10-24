#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 
#include <string.h>

#include "driver.h"
#include "debug.h"
#include "globals.h"
#include "timings.h"
#include "monitor.h"
#include "typeNumbers.h"
#include "arrayWrite.h"


using namespace std;
static const double pi = 4.0*atan(1.0);



/*******************************************************************************
 Title: Shu's Linear Test
 Description: Gives the value at location x for the piecewise function used in Shu's Linear Test
 */

long double shutest(long double x){
 long  double value;

 long double a = 0.5L;
 long double z = 0.3L;
 long double delta = 0.005L;
 long double alpha = 10.0L;
 long double beta = log( 2.0L)/(36.0L*delta*delta);

  if( 0.2 <= x && x <= 0.4){
    value = 1.0L/6.0L * (exp(-beta*((x)-(z+delta))*((x)-(z+delta)))+
			 exp(-beta*((x)-(z-delta))*((x)-(z-delta)))+
			 4.0L*exp(-beta*((x)-(z))*((x)-(z))));
    }
  else if (0.6 <= x && x <= 0.8){
    value = 1L ;}
  else if (1.0 <= x && x <= 1.2){
      value = 1-fabs(10.0L*(x-1.1L));}
  else if( 1.4 <= x && x <= 1.6){
    value = 1.0L/6.0L * (sqrt(max(1-alpha*alpha *((x-1)-(a-delta))*((x-1)-(a-delta)),0.0L))+
			  sqrt(max(1-alpha*alpha *((x-1)-(a+delta))*((x-1)-(a+delta)),0.0L))+
			  4.0L*sqrt(max(1-alpha*alpha *((x-1)-a)*((x-1)-a),0.0L)));
    }
    else value=0L;

    return value;
}


/******************************************************************************
Title: Initial Condition
Description: Computes the average value of the initial condition on the element ix.
*/
double InitialCondition(int ix){

	double value = 0.0;
	if(testCase == 1 || testCase == 4){ //Sine wave
		value =2*sin(pi*dx/2)*sin(pi*(grid(ix)+dx/2))/(pi*dx);
	}
	if(testCase ==2){ //Line
		value = (grid(ix)+dx/2.0);
	}
	if(testCase == 5){ //Shu's Linear Tesst

			   value = 0.0L;
			    long double quadPt[5];
			    long double quadWt[5];
			    quadPt[0] =  0.0L;
			    quadPt[1] =  1.0L/3.0L * sqrt( 5.0L - 2.0L * sqrt( 10.0L/7.0L ) );
			    quadPt[2] = -1.0L/3.0L * sqrt( 5.0L - 2.0L * sqrt( 10.0L/7.0L ) );
			    quadPt[3] =  1.0L/3.0L * sqrt( 5.0L + 2.0L * sqrt( 10.0L/7.0L ) );
			    quadPt[4] = -1.0L/3.0L * sqrt( 5.0L + 2.0L * sqrt( 10.0L/7.0L ) );

			    quadWt[0] = 128.0L/225.0L;
			    quadWt[1] = ( 322.0L + 13.0L*sqrt(70.0L) ) / 900.0L;
			    quadWt[2] = ( 322.0L + 13.0L*sqrt(70.0L) ) / 900.0L;
			    quadWt[3] = ( 322.0L - 13.0L*sqrt(70.0L) ) / 900.0L;
			    quadWt[4] = ( 322.0L - 13.0L*sqrt(70.0L) ) / 900.0L;
			    for(int i=0; i<5; i++){
			      value +=1.0L/2.0L*quadWt[i]*shutest(dx/2.0L*quadPt[i]+grid(ix)+dx/2.0L); }

		}
	if(testCase == 9 ){ // Pure Shock for Buckley Leverett Flux
		if( grid(ix)<0.5 ){
			value = 1.0/sqrt(2);
		}
        else{value = 0.0;}
	}
    
    
	return value;

}

/******************************************************************************
Title: Analytical Solution
Description: Computes the average value of the  analytical solution on the element ix at a given time.
 */
double AnalyticalSolution(double time, int timeStep,int ix){

	double value = 0.0;
	if(testCase == 1){
		value = sin(pi*(grid(ix)+dx/2+timeStep%2*dx/2.0-time));
	}
    else{cout<<"NO ANALYTICAL SOLUTION COMPUTED!"<<endl;}

  return value;
}



/******************************************************************************
 Title: flux
 Description: Gives the flux function. We have u_t+f(u)_x=0, gives f(u)
*/
double flux(double uloc,  double tloc ){
  if(testCase == 1) return uloc;
  if(testCase == 2) return uloc;
  if(testCase == 4) return uloc *uloc/2.0;
  if(testCase == 5) return uloc;
  if(testCase == 9) return (uloc*uloc)/(uloc*uloc+(1-uloc)*(1-uloc));
  return uloc;
  // return 0.0L;
}


/******************************************************************************
 Title: Boundary Condition (bc)
 Description: The boundary type = 0 is periodic boundary condition.  Boundary type  = 1 is reproduce the inerior edge.   This portion of the code tells us what cell ii is taking into acount the boundary type.
 */
/*******************************************************************************/
int bc(int ii){
    if (boundaryType == 0) return (ii+nx)%nx;
    if(boundaryType ==1){
        if(ii<0) return 0;
        if(ii>=nx) return nx-1;
    }
    return ii;
}


/******************************************************************************
Title: WENOReconstructAvg
Description: Computes a 3rd order WENO reconstruction of the average value on the left half of the cell (ix).
*/
static double WENOReconstructAvg(int ix){

	 double epsilon=1e-6; //Used to prevent division by 0.
	 int p=2;
	 int p_ord = (order+1)/2; //order (or degree + 1) of lower dimensional polynomial

	Array1 <double> concLoc(order);
	for(int i =0; i<order; i++){concLoc(i) = concOld(bc(ix -(p_ord-1)+i));}

	Array1<double> IS(p_ord);
	Array1<double> interpolant(p_ord);
	Array1<double> linWeight(p_ord);
	Array1<double> weight(p_ord);


	//Get Smoothness indicators and Polynomial interpolants, and linear weighting
		if(weightType == 0 ){    for(int k =0; k<p_ord; k++)  IS(k)= 1.0L; }
		else{
			IS(0) = (concLoc(0)-concLoc(1))*(concLoc(0)-concLoc(1));
			IS(1) = (concLoc(2)-concLoc(1))*(concLoc(2)-concLoc(1));
		}

		interpolant(0)=(double)1.0/4.0*(concLoc(0)+3.0L*concLoc(1)) ;
		interpolant(1)=(double) 1.0/4.0*(5.0L*concLoc(1)-concLoc(2));

		linWeight(0) = 1.0/2.0;
		linWeight(1) = 1.0/2.0;

	//Compute the nonlinear weights
	for(int i=0; i<p_ord; i++) weight(i) = linWeight(i)/(pow(epsilon+IS(i),p));
    
	//Normalize the weights
	double sum =0.0;
	for(int i =0; i<p_ord; i++) sum +=  weight(i);
	for(int i=0; i<p_ord; i++) weight(i)=weight(i)/sum;
	
    //Compute the reconstruction of the average value.
	double reconstruct = 0.0L;
	for(int i=0; i<p_ord; i++) reconstruct += weight(i)*interpolant(i);


	return reconstruct;
}

/******************************************************************************
Title: redistributeMass
Description: project mass onto a grid shifted by dx/2.  gridIndicator 0: backward shift (time step 1 -> time step 2), and gridIndicator = 1: forward shift(time step 0 -> time step 1).
*/

static void redistributeMass(int gridIndicator){
	switch(gridIndicator) {
	case 0: {
		double conc0 = 2*concOld(nx-1)-WENOReconstructAvg(-1);
		double conc1 = WENOReconstructAvg(0);
		for(int ix = 0; ix<nx; ix++){
			conc(ix)=0.5L*(conc0+conc1);
			conc0 = 2*concOld(ix)-conc1;
			conc1 = WENOReconstructAvg(ix+1);
		}
		return;
	}//end grid indicator 0

	case 1: {
		double conc0= 2*concOld(0)-WENOReconstructAvg(0);
		double conc1 =WENOReconstructAvg(1);
		for(int ix = 0; ix<nx; ix++){
			conc(ix)=0.5L*(conc0+conc1);
			conc0 = 2*concOld(bc(ix+1))-conc1;
			conc1 = WENOReconstructAvg(ix+2);

		}
		return;
	}//end grid indicator 1

	default: {
		cerr << "ERROR IN CODE: BAD GRIDINDICATOR" << endl;
	}

	}//end grid indicator

	return;
}// end redistribute mass

/******************************************************************************
Title: WENOmd
Description: Calculates the mass at the midpoint of cell ix using the LPR method.  ConcOption = 0 use concOld (cell averages from previous time step).  ConcOption = 1 use conc (cell averages from current time step).
 */

static double WENOmd(int ix, int concOption){
//Ususal Constants for WENO Scheme
	 double epsilon=1e-6;
	 int p=2;

//Get cell averages from surrounding cells
	Array1 <double> concLoc(order);
	if(concOption ==0){
	for(int i =0; i<order; i++){ concLoc(i) = concOld(bc(ix -1+i));}}
	else {for(int i =0; i<order; i++){ concLoc(i) = conc(bc(ix -1+i));}}

//Get Smoothness indicators and Polynomial interpolants, and linear weighting
	Array1< double> IS(3);
	Array1<double>  interpolant(3);
	Array1<double> linweight(3);
	Array1<double> weight(3);
	linweight(0) =1.0/4.0;
	linweight(1) =1.0/2.0;
	linweight(2) =1.0/4.0;
	interpolant(0) = concLoc(1);
	interpolant(1) = concLoc(1)-1.0/12.0*( concLoc(2)-2.0* concLoc(1)+ concLoc(0));
	interpolant(2) = concLoc(1);
    if(weightType == 0 ){ for(int k =0; k<p_ord; k++)  IS(k)= 1.0L; }
    else{
        IS(0) = pow( concLoc(0)- concLoc(1),2);
        IS(2)= pow( concLoc(1)- concLoc(2),2);
        IS(1)= 13.0/3.0*pow( concLoc(2)-2.0* concLoc(1)+ concLoc(0),2)+1.0/4.0*pow( concLoc(0)- concLoc(2),2);
        }

//Compute the nonlinear weights
	for(int i=0; i<3; i++) weight(i) = linweight(i)/(pow(epsilon+IS(i),p));
    
//Normalize the weights
	double sum =0.0;
	for(int i =0; i<3; i++)sum +=  weight(i);
	for(int i=0; i<3; i++) {weight(i)=weight(i)/sum;
	}

//Compute the value at the midpoint
	double reconstruct = 0.0L;
	for(int i=0; i<3; i++) reconstruct += weight(i)* interpolant(i);

    return reconstruct;

}//END WENOmd
/******************************************************************************
Title: WENOptDeriv
Description: Given 3 point values u0,u1,u2 a distance dx between each other.  This is a WENO Scheme to compute the derivative at the midpoint.
 */
static double WENOptDeriv(double u0, double u1, double u2){

//Get Constants, Smoothness indicators and Polynomial interpolants, and linear weighting
    double epsilon=1e-6;
    int p=2;
	Array1< double> IS(2);
	Array1<double> interpolant(2);
	Array1<double> linWeight(2);
	Array1<double> weight(2);

	linWeight(0) =1.0/2.0;
	linWeight(1) =1.0/2.0;
	interpolant(0) = (u1-u0)/dx;
	interpolant(1) = (u2-u1)/dx;
	if(weightType == 0 ){for(int k =0; k<2; k++)  IS(k)= 1.0L; }
	else{

		IS(0) = pow(u0-u1,2);
			IS(1)= pow(u1-u2,2);
	}

	//Compute the nonlinear weights
	for(int i=0; i<2; i++) weight(i) = linWeight(i)/(pow(epsilon+IS(i),p));

    //Normalize the weights
	double sum =0.0;
	for(int i =0; i<2; i++)sum +=  weight(i);
	for(int i=0; i<2; i++) {weight(i)=weight(i)/sum;
	}

    //Compute the derivative at the midpoint
	double reconstruct = 0.0L;
	for(int i=0; i<2; i++) reconstruct += weight(i)*interpolant(i);

	return reconstruct;
}//END WENOmd

/******************************************************************************
Title: rungeKutta
Description: Runs the forwad time step to get concentration values at the time quadrature points for every grid point.
 */
int rungeKutta( Array1<double> concMid, Array2<double>&
		concQuad, double tLoc){
	
    //Runge Kutta Coefficients
	Array1<double> concLoc(nx);
	Array1<double> fluxLoc(nx);
	Array2<double> gLoc(nx,2);

	//// RK Stage 0
	for (int ix = 0; ix < nx; ix++) concLoc(ix) = concMid(ix);
	for (int ix = 0; ix < nx; ix++) fluxLoc(ix) = flux(concLoc(ix), tLoc);
	for (int ix = 0; ix < nx; ix++) {gLoc(ix,0) = - WENOptDeriv( fluxLoc(bc(ix-1)),fluxLoc(ix), fluxLoc(bc(ix+1))); }

	////RK Stage 1
	tLoc += dt ;
	for (int ix = 0; ix < nx; ix++) concLoc(ix) = concMid(ix) + dt * gLoc(ix,0);
	for (int ix = 0; ix < nx; ix++) fluxLoc(ix) = flux(concLoc(ix),  tLoc);
	for (int ix = 0; ix < nx; ix++) {gLoc(ix,1) =  - WENOptDeriv( fluxLoc(bc(ix-1)),fluxLoc(ix),fluxLoc(bc(ix+1))); }

	//theta gives the gauss quadrature points on the interval [0,1]
	double theta[3] = {0.0, 1.0 / 2.0 , 1.0 };

    //Runge Kutta Coefficients
    double b0[3];
	double b1[3];

	for (int k = 0; k<3; k++) {
		b0[k]=  (-0.5) * theta[k] * theta[k] + theta[k];
		b1[k] = 1.0 / 2.0 * theta[k] * theta[k];
	}

	//Get the concentrations at the time quadrature points
	for(int ix =0; ix < nx; ix++){
	concQuad(ix,0) = concMid(ix) + dt * (b0[0] * gLoc(ix,0) + b1[0] * gLoc(ix,1));
	concQuad(ix,1) = concMid(ix) + dt * (b0[1] * gLoc(ix,0) + b1[1] * gLoc(ix,1));
	concQuad(ix,2) = concMid(ix) + dt * (b0[2] * gLoc(ix,0) + b1[2] * gLoc(ix,1));
	}

	return 0;
}
/******************************************************************************
Title: cweno3
Description: Runs the cweno3 portion of the code for one time step
*/
int cweno3(int timeStep){
    
    double mass0 = 0.0;
    double mass1 = 0.0;
  
    //Save concentration from previous time step and total mass
    for(int ix = 0; ix<nx; ix++)  mass0 += conc(ix)*dx;
    concOld = conc;
  
    //Move mass forward onto staggered grid
    redistributeMass(timeStep%2);
  
    //Get concentrations at the midpoints
    Array1<double> concMid(nx);
    for(int ix = 0; ix<nx; ix++) concMid(ix)=WENOmd(ix,0);

    //Use Runge Kutta to get concentration at Quadrature Points (Simpson's Rule)
    Array2<double> concQuad(nx,3);
    
    //Get cell concentration at quadrature point values using RK
    double tloc = (timeStep-1)*dt;
    rungeKutta(concMid, concQuad, tloc);
    
    //Compute the flux integral
    Array1<double> fluxIntegral(nx);
    for(int ix = 0; ix<nx; ix++) fluxIntegral(ix) = 0.0;
    for(int ix=0; ix<nx; ix++){
	//Quadrature times and weights
        Array1<double> tquad(3);
        tquad(0) = tloc;
        tquad(1) = tloc +1.0/2.0 * dt;
        tquad(2) = tloc + dt;
        Array1<double>wquad(3);
        wquad(0) = 1.0;
        wquad(1) = 4.0;
        wquad(2) = 1.0;
	  for(int k =0; k<3; k++){
		  fluxIntegral(ix) += 1.0/6.0* wquad(k) * flux(concQuad(ix,k),  tquad(k)) ;
	  }
    }

  //Update Concentrations
  for(int ix = 0; ix<nx; ix++){
    if(timeStep%2 == 1){
      conc(ix) += dt/(dx)*(fluxIntegral(ix) - fluxIntegral(bc(ix+1)) ) ;}
    else
      conc(ix) += dt/(dx)*( fluxIntegral(bc(ix-1))-fluxIntegral(ix));
  }
//Adjust for inflow boundary conditions
    if(boundaryType == 1) conc(0)=InitialCondition(0);

    // CHECK MASS CONSERVATION
    for(int ix=0;ix<nx;ix++) mass1 += conc(ix)*dx;
    monitor(4, "Mass loss during time step is = ", (mass0-mass1));
  
  return 0;
}

/******************************************************************************
Title: Driver
Description: Control the simulation
*/
int driver(clock_t* times){

  int error = 0;
  monitor(1,"COMPUTE SOLUTION"); 
  
  // GENERAL SETUP
  monitor(2,"General setup"); 
  
  // Output file
  ofstream os;
  
  if(outOption != OUTFORMAT_NONE) {
    strcpy(fname,dirName);
    switch(outOption) {
    case OUTFORMAT_RAW: {
      strcat(fname,"conc");
      break;
    }
    case OUTFORMAT_TECPLOT: {
      strcat(fname,"conc.plt");
      break;
    }
    }
    os.open(fname,ios::out);
    if(!os) cerr<<"WARNING: failed to open output file "<<fname<<endl;
  }
  

  //Refinement Loop
    double L1ErrorOld;
    double LinfErrorOld;

    //Compute the solutions over nref number of grid refinements
   for (int nref= 0; nref<Nrefine +1; nref++){
   // INITIAL CONDITION for grid size nx
   conc.realloc(nx);
   for(int ix =0; ix<nx; ix++)conc(ix)=InitialCondition(ix);
  // cout<<nx<<endl;
  concOld.realloc(nx);
  concOld = conc;

  // TIME LOOP SETUP
  double time = 0.0; 
  int timeStep = 0;
  
  switch(outOption) {
  case OUTFORMAT_RAW: {
    os << "Concentration at step " << timeStep << ", time " << time << endl;
    os << conc;
    break;
  }
  case OUTFORMAT_TECPLOT: {
    write_tecplot_pt(conc, os, "conc",1, &grid(), 0);
    break;
  }
  }

  //TIME LOOP
  while(time <( timeFinal- dt/2.0)) {
    time += dt;
    timeStep++;
    monitor(3,"Time Step ",timeStep,", to time = ",time); 
    
    //RUN CWENO3 for 1 time step 
    error =  cweno3(timeStep);
    
    //END TIME LOOP
  }/*end time step*/
  
  // OUTPUT
   if(outOption != OUTFORMAT_NONE) monitor(4,"Writing output");
   switch(outOption) {
   case OUTFORMAT_RAW: {
     os << "\nConcentration at step " << timeStep << ", time " << time << endl;
     os << conc;
     break;
   }
   case OUTFORMAT_TECPLOT: {
     write_tecplot_pt(conc, os, "conc",1, &grid(), 2);
     break;
   }
   }

  // Error ANALYSIS, compares L1 and Linf norms of computational solution with given analytical solution
    if (testCase == 1 ){
      double L1Error = 0.0;
      double LinfError = 0.0;
      for(int ix=0; ix<nx;ix++){
      	double SolnPt;
          	SolnPt = WENOmd( ix,1);
        double solnDiff = fabs(SolnPt - AnalyticalSolution(time,timeStep,ix));
        L1Error += solnDiff*dx;
        LinfError = max(LinfError, solnDiff);
      }/*end loop over elements*/
        
        //Print the Error Results
        cout<<" "<<endl;
        cout<<" Refinement: "<<nref<<", nx: "<<nx<<endl;
        cout<<" L1 Error: "<< L1Error;
        if(nref>0)cout<<",  Rate: "<< log(L1ErrorOld/L1Error)/log(2) ;
        cout<<endl;
        cout<<" Linf Error: "<< LinfError;
        if(nref>0)cout<<",  Rate: "<< log(LinfErrorOld/LinfError)/log(2) ;
        cout<<endl;

        //Save Errors for next refinement
        L1ErrorOld = L1Error;
        LinfErrorOld = LinfError;
    }

  //Refine grid
    nx= nx*2;
    grid.realloc(nx+1);
    grid.setUniform(0,2);
    dx = grid.d(0);
    dt = dt/2;
   }//End refinement loop
  return error;
}/*end function driver*/
