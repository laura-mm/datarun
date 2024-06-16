
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

string carryon;
bool info = false;
double p = 2.0;
double gama;
double mu;
double k;

vector<double> muval;
vector<double> gammaval;
vector<double> sigmaval;

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

double X(double z1)
{
	double pphi = phi(z1);
	return (gama*(p - 1.0)*pphi*pphi/second(z1)) + pphi;
}
// this is zero when phi = 0, so when sigma = 0, also for some gamma that probably is not in the range?

double help(double z1)
{
	double sec = second(z1);
	return sec/(sec + (gama*(p - 1.0)*phi(z1))); // may be better because doesnt rely on other parameters
}
// this is undefined for sigma = 0

double first(double z1)
{
	return (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) + (exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
}
// this is the first order integral, this can not be zero


double M(double z1)
{
	return -k/((z1*help(z1)/first(z1)) + mu); // this could give -infinity
}

double sigtot(double z1)
{
	return help(z1)*M(z1)/first(z1);
}

double q(double z1)
{
	return second(z1)*pow(M(z1)/first(z1), 2.0);
}

double sigma(double z1)
{
	return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1), p - 2.0);
}

double mud2(double z1) // for p = 2 this is the value of mu where we find divergence
{
	return -z1*help(z1)/first(z1);
}

///////////////////////////////// end of equation type stuff

double crit_z1() // depends on p only
{	
	double z1 = -10.0;
	double small = pow(10.0, -6.0);
	double y = second(z1) - (phi(z1)*(p - 1.0));
	double grad;
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = ((second(z1 + small) - (phi(z1 + small)*(p - 1.0))) - (second(z1 - small) - (phi(z1 - small)*(p - 1.0))))/(2.0*small);
		if (grad == 0.0) break;
		znew = z1 - (y/grad);
		ynew = second(znew) - (phi(znew)*(p - 1.0));
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}

double div2_z1(double zstart) // for p = 2, for the set value of mu and gama, finds the value of zd
{
	double z1 = zstart;
	double small = pow(10.0, -6.0);
	double y = z1*help(z1)/first(z1) + mu;
	double grad;
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (((z1 + small)*help(z1 + small)/first(z1 + small)) - ((z1 - small)*help(z1 - small)/first(z1 - small)))/(2.0*small);
		if (grad == 0) {if (abs(y) < pow(10.0, -5.0)) return z1; else return -1.0/0.0;}
		znew = z1 - (y/grad);
		ynew = znew*help(znew)/first(znew) + mu;
		if (info == true) {cout << "znew " << znew << ", ynew " << ynew << endl;}
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	if (gama == -0.5 && mu >= 1.5) return -1.0/0.0;
	else return z1;
}

double sig2_z1(double zstart, double sigg2) // find z1 for given sig2 and set gamma and mu
{
	if (gama >= 0.0 && mu > 1.0) return 1.0/0.0;
	double zdhigh = div2_z1(0.0);
	if (info == true) cout << "zdhigh " << zdhigh << ", " << sigma(zdhigh) << endl;
	double zdlow = -1000.0;
	if (gama < 0.0 && mu >= 1.0) zdlow = div2_z1(-1.0);
	if (info == true) cout << "zdlow " << zdlow << ", " << sigma(zdlow) << endl;
	if (!(abs(sigma(zdlow)) >= 0.0) || sigma(zdlow) < 0.0) return 1.0/0.0;
	if (sigg2 < sigma(zdlow) || sigg2 > sigma(zdhigh))
	{
		if (info == true)
		{
			cout << "sigg2 " << sigg2 << " out of range" << endl;
			//cin >> carryon;
		}
		return 1.0/0.0;
	}
	
	int grid = 1000000;
	double small = pow(10.0, -6.0);
	
	if (info == true)
	{
	ofstream ab; ab.open("abss.txt");
		
	for (int i = 0; i <= grid; i++)
	{
		double z = 100.0*(double)i/(double)grid - 90.0;
		double y = sigma(z) - sigg2;
		ab << z << "," << y;
		if (i != grid) {ab << ",";}
	}
	ab.close();
	cin >> carryon;
	}
	
	double z1 = zstart;
	if (z1 < zdlow) z1 =  zdlow;
	if (z1 > zdhigh) z1 = zdhigh - small;
	double y = sigma(z1) - sigg2;
	double grad;
	if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (sigma(z1 + small) - sigma(z1 - small))/(2.0*small);
		if (grad == 0.0) break;
		if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		if (znew < zdlow) znew =  zdlow;
		if (znew > zdhigh) znew = zdhigh - small;
		ynew = sigma(znew) - sigg2;
		if (abs(ynew) >= abs(y) && abs(y) < pow(10.0, -12.0)) break;
		y = ynew;
		z1 = znew;
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -12.0)) small /= 10.0;
	}
	return z1;	
}

void fiveplot2(int grid, int mi) // plots order parameters for varying sigma
{
	mu = muval[mi];
	
	//cout << endl << mu << endl << endl;

	double zcrit = crit_z1();
	
	for (int gi = 0; gi < 6; gi++)
	{
		cout << mi << "," << gi << endl;
		
		string filename = "measfp2_" + to_string(mi) + "_" + to_string(gi) + ".txt";
		ofstream file; file.open(filename);

		gama = gammaval[gi];
		
		//cout << endl << gama << endl << endl;
	
		for (int i = 0; i <= grid; i++)
		{
			double z1 = ((double)i*100.0/(double)grid) - 50.0;
			if (sigma(z1) < 0.0) continue;
			if (i != 0) file << ",";
			//if (i == grid) cout << sigma(z1) << endl;
			if (M(z1) >= 0.0) file << sigma(z1) << "," << z1 << "," << phi(z1) << "," << M(z1) << "," << q(z1) << "," << help(z1);
			else file << sigma(z1) << "," << z1 << "," << phi(z1) << "," << 1.0/0.0 << "," << q(z1) << "," << help(z1);
		}
		file.close();
	}
}

void bunin2(int grid, int gi) // plots phase diagram and finds sigma_c and sigma_d for fiveplots
{
	gama = gammaval[gi];
	double critz1 = crit_z1();
	double critstop = mud2(critz1);
	
	cout << gi << endl;
	
	ofstream fil; fil.open("mu2_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		mu = ((min(critstop, 2.0) + 5.0)*(double)i/(double)grid) - 5.0;
		if (i != 0) fil << ",";
		fil << mu << "," << sigma(critz1); // critical line
	}
	fil.close();
	
	ofstream file; file.open("bunin2_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		double z1 = 100.0*(double)i/(double)grid - 90.0;
		mu = mud2(z1);
		if (!(sigma(z1) >= 0.0)) break;
		if (i != 0) file << ",";
		file << mu << "," << sigma(z1);
	}
	file.close();
	
	ofstream points; points.open("points2_" + to_string(gi) + ".txt");
	for (int mi = 0; mi < 25; mi ++)
	{
		mu = muval[mi];
		if (mu <= critstop) points << sigma(critz1) << "," << sigma(div2_z1(0.0)) << ",";
		else points << 1.0/0.0 << "," << sigma(div2_z1(0.0)) << ",";
		if (mu >= 1.0 && gama < 0.0) points << sigma(div2_z1(-1.0));
		else points << 1.0/0.0;
		if (mi != 24) points << ",";
	}
	points.close();
}

void histogram(int mi, int gi) // for paper
{
	cout << mi << ", " << gi << endl;
	mu = muval[mi];
	gama = gammaval[gi];
	
	ofstream hfile; hfile.open("hist2fp_" + to_string(mi) + "_" + to_string(gi) + ".txt");
	
	double z1 = -10.0;
	for (int si = 0; si < 21; si ++)
	{
		double sigg2 = sigmaval[si];
		cout << "sigg2 " << sigg2 << endl;
		z1 = sig2_z1(z1, sigg2);
		if (abs(sigma(z1) - sigg2) < pow(10.0, -10.0))
		{
			if (info == true)
			{
				cout << sigma(z1) << ", " << sigg2 << endl;
				cin >> carryon;
			}
			
			double sd = M(z1)/first(z1);
			double mean = -z1*sd;
			
			hfile << sd << "," << mean << "," << z1;
			if (info == true) cout << sd << ", " << mean << ", " << z1 << endl;
		}
		else
		{
			hfile << 1.0/0.0 << "," << 1.0/0.0 << "," << z1;
			z1 = -100.0;
		}
		if (si != 20) hfile << ",";
	}
}			
			


int main()
{
	k = 1.0;
	
	muval = vector<double>(25);
	for (int i = 0; i <= 24; i++) muval[i] = (0.25*(double)i) - 4.0;
	gammaval = vector<double>{-1.0, -0.8, -0.5, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);
	
	//info = true;
	
	//for (int mi = 0; mi < 25; mi++) fiveplot2(100000, mi);
	//for (int gi = 0; gi < 6; gi++) bunin2(100000, gi);
	
	for (int mi = 0; mi < 25; mi ++)
	{
		for (int gi = 0; gi < 6; gi ++)
		{
			histogram(mi, gi);
		}
	}
	
	/*
	gama = -1.0;
	mu =-1.0;
	
	double z1 = 0.0000000001;
	cout << "phi " << phi(z1) << endl;
	cout << "S " << second(z1) << endl;
	cout << "X " << X(z1) << endl;
	cout << "H " << help(z1) << endl;
	cout << "F " << first(z1) << endl;
	cout << "M " << M(z1) << endl;
	cout << "q " << q(z1) << endl;
	cout << "sig " << sigma(z1) << endl;
	cout << "mud " << mud2(z1) << endl;
	cout << "erf " << (erf(z1/sqrt(2.0)) - 1.0) << endl;
	cout << "erf/2 " << (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) << endl;
	cout << "z1H " << z1*help(z1) << endl;
	cout << "sigsum " << sigtot(z1) << endl;
	*/
	
	return 0;
}





