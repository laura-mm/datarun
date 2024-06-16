// This program plots 3d phase diagrams for a combination of 2nd and 3rd order
// the x-y plane is mu2 and mu3, the height is sigmac and sigmad
// these figures come in pairs, one for sigma3 = 0, and for sigma2 = 0
// these can be plotted for various vaues of gamma.

// codes for plots many lines
// A gam2 = 0, gam3 = 0
// B gam2 = -1, gam3 = -0.5
// C gam2 = 1, gam3 = 1
// D gam2 = -1, gam3 = 1
// E gam2 = 1, gam3 = -0.5

// codes for plots plotly surfaces
// codes for plots plotly surfaces
// F gam2 = 0, gam3 = 0
// G gam2 = -1, gam3 = -0.5
// H gam2 = 1, gam3 = 1
// I gam2 = 0.5, gam3 = 0.5
// J gam2 = 0.4, gam3 = 0.4
// K gam2 = -1, gam3 = 1
// L gam2 = 1, gam3 = -0.5
// M gam2 = -0.8, gam3 = -0.4



#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <iomanip>
#include <Eigen/Dense>
//#include <eigen3/Eigen/Dense>
using namespace Eigen;
using namespace std;

IOFormat pp(20, DontAlignCols, ",", ",", "", "", "", "");
IOFormat cc(6, DontAlignCols, ",", ",", "", "", "", "");

// these are global variables
double k;
ArrayXd mu(2);
ArrayXd gama(2);
double z2c;
double z3c;


double phi(double z1)
{
	return (1.0 - erf(z1/sqrt(2.0)))/2.0;
}
// this becomes zero when z1 = infinity, which means sigma = 0

double second(double z1)
{
	return (((z1*z1) + 1.0)*(1.0 - erf(z1/sqrt(2.0)))/2.0) - (z1*exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
	// this is the second order integral, this is never zero, im not sure if this is true!
}

double X(double z1, double help)
{
	return phi(z1)/help;
}
// this is zero when phi = 0, so when sigma = 0, also for some gamma that probably is not in the range?

double first(double z1)
{
	return (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) + (exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
}
// this is the first order integral, this can not be zero

double mutot(double M)
{
	double m = 0.0;
	for (int i = 0; i < 2; i++)
	{
		m += mu(i)*pow(M, i + 1);
	}
	return m;
}

double M(double z1, double help)
{	
	if (mu(1) == 0.0) return -k/((z1*help/first(z1)) + mu(0));
	//else if (pow((z1*help/first(z1)) + mu(0), 2.0) == (4.0*mu(1)*k)) return sqrt(k/mu(1)); // divergent condition
	else return ((z1*help/first(z1)) + mu(0) + sqrt(pow((z1*help/first(z1)) + mu(0), 2.0) - (4.0*mu(1)*k)))/(-2.0*mu(1));
}

double sigtot(double z1, double help)
{
	return help*M(z1, help)/first(z1);
}

double q(double z1, double help)
{
	return second(z1)*pow(M(z1, help)/first(z1), 2.0);
}


ArrayXd sigmap(double z1) // this gives the two values if sig3 = 0 and sig2 = 0
{
	ArrayXd sigma(2);
	for (int i = 0; i < 2; i ++)
	{
		double p = (double)(i + 2);
		double help = second(z1)/(second(z1) + (gama(i)*(p - 1.0)*phi(z1)));
		sigma(i) = help*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1, help), p - 2.0);
	}
	return sigma;
}

ArrayXd sigmap(double z1, double gamm0, double gamm1) // for a different specified gamma
{
	ArrayXd gamm(2); gamm(0) = gamm0; gamm(1) = gamm1;
	ArrayXd sigma(2);
	for (int i = 0; i < 2; i ++)
	{
		double p = (double)(i + 2);
		double help = second(z1)/(second(z1) + (gamm(i)*(p - 1.0)*phi(z1)));
		sigma(i) = help*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1, help), p - 2.0);
	}
	return sigma;
}


double crit_z1(double p) // depends on p only
{
	double z1 = -10.0;
	double small = pow(10.0, -6.0);
	double y = second(z1) - (phi(z1)*(p - 1.0));
	double grad;
	//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = ((second(z1 + small) - (phi(z1 + small)*(p - 1.0))) - (second(z1 - small) - (phi(z1 - small)*(p - 1.0))))/(2.0*small);
		if (grad == 0.0) break;
		//if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = second(znew) - (phi(znew)*(p - 1.0));
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}

double maxsig_z1(double zstart) // returns the z1 where sigma3 reaches a stationary point, depends on mu and gamma, for sig2 = 0
{
	int grid = 1000000;
	double small = pow(10.0, -6.0);
	
	/*if (abs(mu(1)) < 0.001)
	{
	ofstream ab; ab.open("abs.txt");
	cout << "now plotting function for maxsig for mu = " << mu(1) << endl;
	
	for (int i = 0; i <= grid; i++)
	{
		double z = 100.0*(double)i/(double)grid - 90.0;
		double sigrad = (sigmap(z + small, help, 3.0) - sigmap(z - small, help, 3.0))/(2.0*small);
		ab << z << "," << sigmap(z, help, 3.0) << "," << sigrad;
		if (i != grid) {ab << ",";}
	}
	ab.close();
	cout << "finished plotting z" << endl;
	cin >> carryon;
	}*/
	
	
	double z1 = zstart;
	//double small = pow(10.0, -6.0);
	double y = (sigmap(z1 + small)(1) - sigmap(z1 - small)(1))/(2.0*small);
	double grad;
	//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (sigmap(z1 + small)(1) + sigmap(z1 - small)(1) - (2.0*sigmap(z1)(1)))/(small*small);
		if (grad == 0.0) break;
		//if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = (sigmap(znew + small)(1) - sigmap(znew - small)(1))/(2.0*small);
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		//if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}


double noneg(double s) // if negitive then returns infinity
{
	if (s >= 0.0) return s;
	else return 1.0/0.0;
}
double nobig(double s, double big) // if bigger than big returns infinity
{
	if (s <= big) return s;
	else return 1.0/0.0;
}


void plotlyy(int grid1, int grid2, int grid3)
{
	gama(0) = 0.0;
	gama(1) = 0.0;
	
	ofstream crit3; crit3.open("Pcrit233T.txt");
	ofstream div2; div2.open("Pdiv232T.txt");
	ofstream div3; div3.open("Pdiv233T.txt");
	ofstream crit2; crit2.open("Pcrit232T.txt");
	ofstream extra3; extra3.open ("Pextra3T.txt"); // to plot an extra surface to stop it looking squarey
	
	// this part plots the instability line and sigmax line for sig2 = 0
	
	for (int i = 0; i <= grid1; i++)
	{
		//mu(0) = 10.0*(double)i/(double)grid1 - 6.0;
		mu(0) = 10.0*(double)i/(double)grid1 - 4.0;
		double zmax; // = 0.0;
		if (i < 80) zmax = 0.1;
		else zmax = 0.2;
		double zstop = 100.0; // z1 at the highest value of mu3 where sigmax can be found, where to end the divergence lines
		
		int lmstop = -1;
		int lcstop = -1;
		
		for (int l = 0; l <= grid3; l++)
		{
			mu(1) = 11.0*(double)l/(double)grid3 - 4.0;
			//mu(1) = 11.0*(double)l/(double)grid3 - 6.0;
			if (i != 0 || l != 0) crit3 << ",";
			if (sigmap(zmax)(1) > 0.0) zmax = maxsig_z1(zmax);
			//if (i == 89) cout << mu(1) << ", " << zmax << ", " << sigmap(zmax)(1) << endl;
			if (sigmap(zmax)(1) > 0.0)
			{
				//if (i == 0) cout << zmax << endl;
				zstop = zmax;
				lmstop = l;
				if (zmax < z3c) crit3 << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c)(0)) << "," << 1.0/0.0 << "," << noneg(sigmap(zmax)(1)); // dont observe sigma3c
				else // observe sigma3c
				{
					crit3 << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c)(0)) << "," << noneg(sigmap(z3c)(1)) << "," << noneg(sigmap(zmax)(1));
					if (sigmap(z3c)(1) > 0.0) lcstop = l;
					//cout << mu(0) << ", " << mu(1) << ", " << sigmap(z3c)(1) << endl;
				}
			}
			else if (l == lmstop + 1) // so it meets the blue surface
			{
				double help = second(zstop)/(second(zstop) + (gama(1)*2.0*phi(zstop)));
				mu(1) = pow((zstop*help/first(zstop)) + mu(0), 2.0)/(4.0*k); // this is divergent condition
				//cout << i << ", " << zstop << ", " << mu(0) << ". " << mu(1) << ". " << sigmap(zstop)(1) << endl;
				if (zmax < z3c) crit3 << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c)(0)) << "," << 1.0/0.0 << "," << noneg(sigmap(zstop)(1)); // dont observe sigma3c
				else // observe sigma3c
				{
					crit3 << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c)(0)) << "," << noneg(sigmap(z3c)(1)) << "," << noneg(sigmap(zstop)(1));
					if (sigmap(z3c)(1) > 0.0) lcstop = l;
				}
			}
			else
			{
				if (zmax < z3c) crit3 << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c)(0)) << "," << 1.0/0.0 << "," << noneg(sigmap(zmax)(1)); // dont observe sigma3c
				else // observe sigma3c
				{
					crit3 << mu(0) << "," << mu(1) << "," << noneg(sigmap(z2c)(0)) << "," << noneg(sigmap(z3c)(1)) << "," << noneg(sigmap(zmax)(1));
					if (sigmap(z3c)(1) > 0.0) lcstop = l;
					//cout << mu(0) << "; " << mu(1) << "; " << sigmap(z3c)(1) << endl;
				}
			}

		}
		
		cout << i << ", " << lmstop << endl;
		
		if (i != 0) extra3 << ",";
		if (lmstop == grid3) {extra3 << mu(0) << "," << lcstop << "," << lmstop; cout << "a" << endl;}
		else {extra3 << mu(0) << "," << lcstop << "," << lmstop + 1; cout << "b" << endl;}
		
		/*if (i == 0)
		{
			cout << "zstop " << zstop << endl;
			double p = 3.0;
			double help = second(zstop)/(second(zstop) + (gama(0)*(p - 1.0)*phi(zstop)));
			mu(1) = pow((zstop*help/first(zstop)) + mu(0), 2.0)/(4.0*k); // this is divergent condition
			cout << "mu3 " << mu(1) << endl;
			cout << "sigma " << sigmap(zstop)(1) << endl;
			zmax = maxsig_z1(zstop);
			cout << "zmax " << zmax << endl;
			cout << "sigma " << sigmap(zmax)(1) << endl;
			help = second(zmax)/(second(zmax) + (gama(0)*(p - 1.0)*phi(zmax)));
			mu(1) = pow((zmax*help/first(zmax)) + mu(0), 2.0)/(4.0*k); // this is divergent condition
			cout << "mu3 " << mu(1) << endl;
			cout << "sigma " << sigmap(zmax)(1) << endl;
		}*/
			
		
			
		
		// this part plots the divergence line for sig3 = 0
	
		double intersection2;
		double p = 2.0;
		bool fold = false; // div graph folds over, has it past the fold point
		double minmu = 10.0; // to find the fold point
		// the divergent line curves back on itself, but when crossed twice the M is negitive, so instead of curving back we only plot if M is positive, and artifical straight line
		
		for (int j = 0; j <= grid2; j++)
		{
			double z1 = 50.0*(double)j/(double)grid2 - 20.0;
			//double z1 = -exp(-1.0*(20.0*(double)j/(double)grid2 - 4.3)) + 20.0;
			double help = second(z1)/(second(z1) + (gama(0)*(p - 1.0)*phi(z1)));		
			mu(1) = pow((z1*help/first(z1)) + mu(0), 2.0)/(4.0*k); // this is divergent condition
			if (mu(1) <= minmu) minmu = mu(1);
			else if (z1 >= z2c) fold = true; // this condition is because there were some issues where it was setting to zero prematurely, needed to be commented out for E
			if ((i != 0 || j != 1) && j != 0) div2 << ",";
			if (sigmap(z1)(0) >= 0.0) // && mu(1) < 7.0 && mu(1) > -5.0)
			{
				if (sigmap(z1)(1) < 0.0 || fold == true) // artificial straight line
				{
					mu(1) = 0.0;
					if (j != 0) div2 << mu(0) << "," << mu(1) << "," << nobig(sigmap(z1)(0), 10.0);
					//if (i == 52) cout << z1 << ", " << mu(0) << ", " << mu(1) << ", " << sigmap(z1)(0) << endl;
					//if (abs(z2c - z1) < pow(10.0, -15.0)) {intersection2 = mu(1)}//; if (i == 0) cout << "intersection " << z1 << ", " << mu(1) << endl;}
				}
				else
				{
					if ((mu(1) > 7.0 || mu(1) < -5.0) && j != 0) div2 << mu(0) << "," << nobig(mu(1), 7.0) << "," << 1.0/0.0;
					//if ((mu(1) > 5.0 || mu(1) < -7.0) && j != 0) div2 << mu(0) << "," << mu(1) << "," << 1.0/0.0;
					else if (j != 0) div2 << mu(0) << "," << nobig(mu(1), 7.0) << "," << nobig(sigmap(z1)(0), 10.0);
					//if (i == 52) cout << z1 << ". " << mu(0) << ". " << mu(1) << ". " << sigmap(z1)(0) << endl;
					//if (abs(z2c - z1) < pow(10.0, -15.0)) {intersection2 = mu(1);}// if (i == 0) cout << "intersection " << z1 << ", " << mu(1) << endl;}
				}
				if (abs(z2c - z1) < pow(10.0, -15.0)) intersection2 = mu(1);
			}
			else if (j != 0) div2 << 1.0/0.0 << "," << 1.0/0.0 << "," << 1.0/0.0;
				
		}
	
		//if (i == 0) cout << endl << endl;
	
		// this part plots the divergence line for sig2 = 0

		p = 3.0;
		
		int jdstart = -1;
		int jdstop = -1;
		
		for (int j = 0; j <= grid2; j++)
		{
			//double z1 = -1.0/((1000.0)*(double)j/(double)grid2); // for antisymmetirc
			//double z1 = (70.0)*(double)j/(double)grid2 - 50.0;
			double z1 = (zstop + 50.0)*(double)j/(double)grid2 - 50.0;
			//double z1 = -exp(-1.0*(20.0*(double)j/(double)grid2 - 4.0));
			double help = second(z1)/(second(z1) + (gama(1)*(p - 1.0)*phi(z1)));
			mu(1) = pow((z1*help/first(z1)) + mu(0), 2.0)/(4.0*k); // this is divergent condition
			if (i != 0 || j != 0) div3 << ",";
			if (sigmap(z1)(0) >= 0.0 && sigmap(z1)(1) >= 0.0 && sigmap(z1)(1) <= 10.0 && mu(1) <= 7.0)
			//if (sigmap(z1)(0) >= 0.0 && sigmap(z1)(1) >= 0.0 && sigmap(z1)(1) <= 10.0 && mu(1) <= 5.0)
			{
				if (jdstart == -1) jdstart = j;
				jdstop = j;
				div3 << mu(0) << "," << nobig(mu(1), 7.0) << "," << nobig(sigmap(z1)(1), 10.0);
				//div3 << mu(0) << "," << nobig(mu(1), 5.0) << "," << nobig(sigmap(z1)(1), 10.0);
				//if (j == grid2) cout << "j " << j << ", z " << z1 << ", mu3 " << mu(1) << ", sig " << sigmap(z1)(1) << endl;
			}
			else {div3 << mu(0) << "," << 1.0/0.0 << "," << 1.0/0.0;}// if (i == 71) cout << "j " << j << ", z " << z1 << ", mu3 " << mu(1) << ", sig " << sigmap(z1)(1) << endl;}
		}
		
		extra3 << "," << jdstart << "," << jdstop;
		
		// this part plots the instability line for sig3 = 0
		mu(1) = intersection2;
		if (i != 0) crit2 << ",";
		crit2 << mu(0) << "," << mu(1) << "," << sigmap(z2c)(0);
		
			
			
	}	
		crit3.close();
		div2.close();
		div3.close();
		crit2.close();
		extra3.close();
}








int main()
{
	gama = ArrayXd::Zero(2);
	mu = ArrayXd::Zero(2);
	k = 1.0;
	z2c = crit_z1(2.0);
	z3c = crit_z1(3.0);
	//cout << z3c << endl;
	
	//phase23(100, 14000, 10000);
	
	//plotlyy(100, 14000, 1500);
	//plotlyy(4, 4, 4);
	/*
	int grid2 = 10000;
	
	gama(0) = 0.0;
	gama(1) = 0.0;
	mu(0) = 0.0;
	mu(1) = -1.0;
	double z1 = 1.086;
	double p = 2.0;
	double help = second(z1)/(second(z1) + (gama(0)*(p - 1.0)*phi(z1)));
	cout << sigmap(z1)(0) << ", " << M(z1, help) << ", " << q(z1, help) << ", " << sigtot(z1, help) << ", " << k + mutot(M(z1, help)) << endl;
	*/
	
	/*for (int j = 0; j <= grid2; j++)
	{
		double z1 = 70.0*(double)j/(double)grid2 - 50.0;
		cout << z1 << ", " << sigmap(z1)(0) << endl;
	}*/
	
	
	
	/*
	double z1 = -0.0005;
	cout << phi(z1) << "," << first(z1) << "," << second(z1) << endl;
	double p = 2.0;
	double help = second(z1)/(second(z1) + (gama(0)*(p - 1.0)*phi(z1)));
	cout << "help " << help << endl;
	mu(1) = pow((z1*help/first(z1)) + mu(0), 2.0)/(4.0*k); 
	cout << "mu3 " << mu(1) << endl;
	cout << "M " << M(z1, help) << endl;
	cout << "MM " << sqrt(k/mu(1)) << endl;
		//else return ((z1*help/first(z1)) + mu(0) + sqrt(pow((z1*help/first(z1)) + mu(0), 2.0) - (4.0*mu(1)*k)))/(-2.0*mu(1));
	cout << z1*help/first(z1) << endl;
	cout << "sigma2 " << help*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/sqrt(k/mu(1)), p - 2.0) << endl;
	cout << "sigma2c " << sigmap(z2c)(0) << endl;
	p = 3.0;
	help = second(z1)/(second(z1) + (gama(1)*(p - 1.0)*phi(z1)));
	cout << "help " << help << endl;
	mu(1) = pow((z1*help/first(z1)) + mu(0), 2.0)/(4.0*k); // this is divergent condition
	cout << "mu3 " << mu(1) << endl;
	cout << "M " << M(z1, help) << endl;
	cout << "MM " << sqrt(k/mu(1)) << endl;
	cout << "sigma3 " << help*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/sqrt(k/mu(1)), p - 2.0) << endl;
	cout << "sigma3c " << sigmap(z3c)(1) << endl;
	*/
	
	
	gama(0) = -0.8;
	gama(1) = -0.4;
	mu(0) = 1.25;
	mu(1) = 0.4;
	
	double z1 = z2c;
	
	cout << sigmap(z1)(0) << ", " << sigmap(z1)(1) << endl;
	
	
	return 0;
}
