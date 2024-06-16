
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <iomanip>
#include <Eigen/Dense>
#include <complex>
//#include <eigen3/Eigen/Dense>
using namespace Eigen;
using namespace std;

IOFormat pp(20, DontAlignCols, ",", ",", "", "", "", "");
IOFormat cc(6, DontAlignCols, ",", ",", "", "", "", "");

string carryon;
bool info = false;
double p = 4.0;
double gama;
double mu;
double k;
bool broke;
bool joined; // from mu = 0 we are looking for the last value of z1 where sigma is defined, then as mu increases this gap shrinks until it joins to another branch, after this happens we use newton raphson again instead of lastsig        

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
	//return phi(z1)/X(z1);
	double sec = second(z1);
	return sec/(sec + (gama*(p - 1.0)*phi(z1))); // may be better because doesnt rely on other parameters
}
// this is undefined for sigma = 0

double first(double z1)
{
	return (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) + (exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
}
// this is the first order integral, this can not be zero

complex<double> crt(complex<double> w)
{
	if (w.imag() == 0.0 && w.real() < 0.0) return -pow(-w, 1.0/3.0);
	else return pow(w, 1.0/3.0);
}
double M(double z1)
{	
	if (mu == 0.0) return -k*first(z1)/(z1*help(z1));
	else
	{
		double Q = z1*help(z1)/(3.0*mu*first(z1)); //cout << "Q " << Q << endl;
		double R = -k/(2.0*mu); //cout << "R " << R << endl;
		double D = pow(Q, 3.0) + pow(R, 2.0); //cout << "D " << D << endl;
		complex<double> root;
		if (D >= 0.0) root = pow(D, 0.5);
		else root = pow(-D, 0.5)*1i;
		complex<double> S = crt(R + root); //cout << "S " << S << endl;
		complex<double> T = crt(R - root); //cout << "T " << T << endl;
		complex<double> Trotate = exp(2i*M_PI/3.0);
		complex<double> M1 = S + T; //cout << "M1 " << M1 << endl;
		complex<double> M2 = -0.5*(S + T) + 0.5i*pow(3, 0.5)*(S - T); //cout << "M2 " << M2 << endl;
		complex<double> M3 = -0.5*(S + T) - 0.5i*pow(3, 0.5)*(S - T); //cout << "M3 " << M3 << endl;
		if (mu < 0.0) return M1.real();
		else if (mu > 0.0) //return M3.real();
		{
			if (M3.imag() == 0.0) return M3.real();
			else return pow(-1.0, 0.5);
		}
	}
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
	//if (info == true) cout << help(z1) << ", " << second(z1) << ", " << first(z1) << ", " << M(z1) << endl;
	return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1), p - 2.0);
}

double mud4(double z1) // for p = 4, value of mu where you get 3 real solutions
{		
	return pow(2.0/k, 2.0)*pow((-z1*help(z1))/(3.0*first(z1)), 3.0);
}


///////////////////////////////// end of equation type stuff


double crit_z1() // depends on p only
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

double maxsig_z1(double zstart) // returns the z1 where sigma reaches a stationary point, depends on mu and gamma
{
	double small = pow(10.0, -6.0);
	if (info == true)
	{
	ofstream ab; ab.open("abs.txt");
	int grid = 1000000;

	cout << "now plotting function for maxsig for mu = " << mu << " gamma = " << gama << endl;
	
	for (int i = 0; i <= grid; i++)
	{
		double z = 100.0*(double)i/(double)grid - 90.0;
		double sigrad = (sigma(z + small) - sigma(z - small))/(2.0*small);
		ab << z << "," << sigma(z) << "," << sigrad;
		//cout << z << "," << sigma(z) << "," << sigrad << endl;
		if (i != grid) {ab << ",";}
	}
	ab.close();
	cout << "finished plotting z" << endl;
	}
	double z1 = zstart;
	//double small = pow(10.0, -6.0);
	if (info == true) cout << "z1 = " << z1 << ", bigsig = " << sigma(z1 + small) << ", smallsig = " << sigma(z1 - small) << endl;
	double y = (sigma(z1 + small) - sigma(z1 - small))/(2.0*small);
	double grad;
	if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0 && y < 1.0/0.0)
	{
		grad = (sigma(z1 + small) + sigma(z1 - small) - (2.0*sigma(z1)))/(small*small);
		if (grad == 0.0) break;
		if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = (sigma(znew + small) - sigma(znew - small))/(2.0*small);
		if (info == true) {cout << "znew " << znew << ", ynew " << ynew << endl;}
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}

void fiveplot4(int grid, int mi) // plots order parameters for varying sigma
{
	mu = muval[mi]; //
	//mu = -50.0;
	
	double zcrit = crit_z1();
	//cout << "crit " << zcrit << endl;
	
	for (int gi = 0; gi < 6; gi++)
	{
	
		//int gi = 1;
		
		cout << mi << "," << gi << endl;
		
		string filename = "measfpb4_" + to_string(mi) + "_" + to_string(gi)+ ".txt";
		ofstream file; file.open(filename);

		gama = gammaval[gi];
	
		double sig_prev = -1.0; // next sigma has to be bigger than previous
		int flag = 0;

		for (int i = 0; i <= grid; i++)
		{
			//cin >> carryon;
			double z1 = ((double)i*100.0/(double)grid) - 99.0;
			if (sigma(z1) < 0.0) continue; //{cout << "sigma < 0 " << z1 << endl; continue;}
			if (M(z1) < 0.0) continue;
			if (!(abs(sigma(z1)) >= 0.0)) continue; //{cout << "sigma nan " << z1 << endl; continue;}
			//if (sigma(z1) < sig_prev) break; //{cout << "sigma increase " << z1 << endl; break;} // this one!
			sig_prev = sigma(z1);
			if (flag == 1) file << ",";
			//cout << "i " << i << ", z1 " << z1 << ", phi " << phi(z1) << ", M " << M(z1) << ", sig " << sigma(z1) << ", q " << q(z1) << endl;
			file << sigma(z1) << "," << z1 << "," << phi(z1) << "," << M(z1) << "," << q(z1) << "," << help(z1);
			flag = 1;
		}
		file.close();
	}
}

double last_z1(double zstart)
{
	if (info == true)
	{
	ofstream ab; ab.open("abs.txt");
	int grid = 1000000;

	cout << "now plotting function for maxsig for mu = " << mu << " gamma = " << gama << endl;
	double small = pow(10.0, -6.0);
	for (int i = 0; i <= grid; i++)
	{
		double z = 100.0*(double)i/(double)grid - 90.0;
		double sigrad = (sigma(z + small) - sigma(z - small))/(2.0*small);
		ab << z << "," << sigma(z) << "," << sigrad;
		cout << z << "," << sigma(z) << "," << sigrad << endl;
		if (i != grid) {ab << ",";}
	}
	ab.close();
	cout << "finished plotting z" << endl;
	}
	double z1 = zstart;
	while (abs(sigma(z1) >= 0.0))
	{
		z1 += pow(10.0, -6.0);
		if (z1 >= 0.0)
		{
			joined = true;
			return maxsig_z1(zstart);
		}
	} // (fixed) this is clumsy, better to use a switch if we should now use maxsig instead of last sig
	z1 -= pow(10.0, -6.0);
	return z1;
}
		


void bunin4(int grid, int gi)
{
	gama = gammaval[gi];
	double critz1 = crit_z1();
	cout << critz1 << endl;
	double zstart;
	if (gi == 0) zstart = 0.0;
	else if (gama < 0.0) zstart = -0.05;
	else zstart = critz1;
	double zstop;
	
	ofstream points; points.open("points4_" + to_string(gi) + ".txt");
	
	ofstream fil; fil.open("mu4_" + to_string(gi) + ".txt");
	info = false;
	joined = false;
	for (int i = 0; i <= grid; i++)
	{
		mu = (10.0*(double)i/(double)grid) - 7.0;
		
		if (i != 0) fil << ",";
		
		zstart = maxsig_z1(zstart);
		
		if (info == true) cin >> carryon;
		
		
		if (zstart < critz1) fil << mu << "," << 1.0/0.0 << "," << sigma(zstart); // no crit and maxsig line
		else fil << mu << "," << sigma(critz1) << "," << sigma(zstart); // critical line and maxsig line
		
		
		
		//cout << "mu = " << mu << ", zstart = " << zstart << ", maxsig = " << sigma(zstart) << ", sigc = " << sigma(z1) << endl;
		for (int mi = 0; mi < 25; mi ++)
		{
			if (abs(mu - muval[mi]) < pow(10.0, -15.0)) cout << mu << ", " << abs(mu - muval[mi]) << endl;
			if (abs(mu - muval[mi]) < pow(10.0, -15.0))
			{
				//cout << mu << endl;
				if (zstart < critz1) points << 1.0/0.0 << "," << sigma(zstart) << ",";
				else points << sigma(critz1) << "," << sigma(zstart) << ",";
				double zlow;
				for (int i = 0; i <= grid; i++)
				{
					double z1 = 100.0*(double)i/(double)grid - 50.0;
					if (sigma(z1) > 0.0 && sigma(z1) < 1.0/0.0) {zlow = z1; break;}
				}
				if (abs(mud4(zlow) - mu) < 0.001) points << sigma(zlow);
				else points << 1.0/0.0;
				if (mi != 24) points << ",";
			}
		}
				
		if (sigma(zstart) > 0.0) zstop = zstart;
		//cout << "sig(zstop) = " << sigma(zstop) << endl;
	}
	fil.close();
	points.close();
	//cout << "zstop = " << zstop << endl;
	ofstream file; file.open("bunin4_" + to_string(gi) + ".txt"); 
	for (int i = 0; i <= grid; i++)
	{
		double z1 = (50.0 + zstop)*(double)i/(double)grid - 50.0;
		//double z1 = 110.0*(double)i/double(grid) - 100.0;
		mu = mud4(z1);
		//double MM = pow(4.0*k/mu, 1.0/3.0);
		double MM = pow(0.5*k/mu, 1.0/3.0);
		double sigmaa = help(z1)*sqrt(1.0/(2.0*pow(second(z1), 3.0)))*pow(first(z1)/MM, 2.0);
		if (!(sigmaa >= 0.0)) break; // this one!
		if (i != 0) file << ",";
		file << mu << "," << sigmaa;
		//cout << z1 << ", " << mu << "," << sigma(crit_z1()) << "," << sigma(z1) << endl;
	}
	file.close();
	
}


void z1_plot(int gi, int mi, int grid) // plots how things vary with z1
{
	mu = muval[mi];
	gama = gammaval[gi];
	bool flag = false;
	ofstream z1plot; z1plot.open("z1plot4_" + to_string(gi) + "_" + to_string(mi) + ".txt");
		
	for (int i = 0; i <= grid; i++)
	{
		double z1 = ((double)i*12.0/(double)grid) - 10.0;
		//if (!(abs(M(z1)) > 0.0) || !(abs(M2(z1)))) continue;
		if (flag == true) z1plot << ",";
		flag = true;
		
		complex<double> M1;
		complex<double> M2;
		complex<double> M3;
		
		if (mu == 0.0) {M1 = -k*first(z1)/(z1*help(z1)); M2 = -k*first(z1)/(z1*help(z1)); M3 = -k*first(z1)/(z1*help(z1));}
		
		else
		{
			double Q = z1*help(z1)/(3.0*mu*first(z1));
			double R = -k/(2.0*mu);
			double D = pow(Q, 3.0) + pow(R, 2.0);
			complex<double> root;
			if (D >= 0.0) root = pow(D, 0.5);
			else root = pow(-D, 0.5)*1i;
			complex<double> S = crt(R + root);
			complex<double> T = crt(R - root);
			complex<double> Trotate = exp(2i*M_PI/3.0);
			M1 = S + T; //cout << "M1 " << M1 << endl;
			M2 = -0.5*(S + T) + 0.5i*pow(3, 0.5)*(S - T); //cout << "M2 " << M2 << endl;
			M3 = -0.5*(S + T) - 0.5i*pow(3, 0.5)*(S - T); //cout << "M3 " << M3 << endl;
		}	
		
		complex<double> q1 = second(z1)*pow(M1/first(z1), 2.0);
		complex<double> q2 = second(z1)*pow(M2/first(z1), 2.0);
		complex<double> q3 = second(z1)*pow(M3/first(z1), 2.0);
		
		complex<double> sig1 = help(z1)*sqrt(2.0/(4.0*pow(second(z1), 3.0)))*pow(first(z1)/M1, 2.0);
		complex<double> sig2 = help(z1)*sqrt(2.0/(4.0*pow(second(z1), 3.0)))*pow(first(z1)/M2, 2.0);
		complex<double> sig3 = help(z1)*sqrt(2.0/(4.0*pow(second(z1), 3.0)))*pow(first(z1)/M3, 2.0);
	
			
		z1plot << z1 << "," << phi(z1) << ",";
		z1plot << M1.real() << "," << M1.imag() << "," << q1.real() << "," << q1.imag() << "," << sig1.real() << "," << sig1.imag() << ",";
		z1plot << M2.real() << "," << M2.imag() << "," << q2.real() << "," << q2.imag() << "," << sig2.real() << "," << sig2.imag() << ",";
		z1plot << M3.real() << "," << M3.imag() << "," << q3.real() << "," << q3.imag() << "," << sig3.real() << "," << sig3.imag();	
		
	}
	z1plot.close();
}

void mu_plot( int gi, int grid) // plot all solutions of sigma for varying z1
{
	gama = gammaval[gi];
	ofstream mudplot; mudplot.open("mudplot4_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		double z1 = ((double)i*14.0/(double)grid) - 13.9;
		mu = mud4(z1);
		double MM = pow(0.5*k/mu, 1.0/3.0);
		double sigmaa = help(z1)*sqrt(1.0/(2.0*pow(second(z1), 3.0)))*pow(first(z1)/MM, 2.0);
		mudplot << z1 << "," << mu << "," << sigmaa;
		if (i != grid) mudplot << ",";
	}
	mudplot.close();
	double z1 = crit_z1();
	ofstream mucplot; mucplot.open("mucplot4_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		mu = ((double)i*5.0/(double)grid) - 3.0;

		complex<double> M1;
		complex<double> M2;
		complex<double> M3;
		
		if (mu == 0.0) {M1 = -k*first(z1)/(z1*help(z1)); M2 = 1.0/0.0; M3 = 1.0/0.0;}
		
		else
		{
			double Q = z1*help(z1)/(3.0*mu*first(z1));
			double R = -k/(2.0*mu);
			double D = pow(Q, 3.0) + pow(R, 2.0);
			complex<double> root;
			if (D >= 0.0) root = pow(D, 0.5);
			else root = pow(-D, 0.5)*1i;
			complex<double> S = crt(R + root);
			complex<double> T = crt(R - root);
			M1 = S + T; //cout << "M1 " << M1 << endl;
			M2 = -0.5*(S + T) + 0.5i*pow(3, 0.5)*(S - T); //cout << "M2 " << M2 << endl;
			M3 = -0.5*(S + T) - 0.5i*pow(3, 0.5)*(S - T); //cout << "M3 " << M3 << endl;
		}	


		complex<double> sig1 = help(z1)*sqrt(2.0/(4.0*pow(second(z1), 3.0)))*pow(first(z1)/M1, 2.0);
		complex<double> sig2 = help(z1)*sqrt(2.0/(4.0*pow(second(z1), 3.0)))*pow(first(z1)/M2, 2.0);
		complex<double> sig3 = help(z1)*sqrt(2.0/(4.0*pow(second(z1), 3.0)))*pow(first(z1)/M3, 2.0);
		
		/*
		mucplot << mu << "," << sig1.real() << ",";
		if (sig2.imag() == 0.0) mucplot << sig2.real() << ",";
		else mucplot << 1.0/0.0 << ",";
		if (sig3.imag() == 0.0) {mucplot << sig3.real(); cout << sig1.real() << ", " << sig2.real() << ", " << sig2.real() << endl;}
		else mucplot << 1.0/0.0;
		*/
		
		mucplot << mu << "," << sig1.real() << "," << sig1.imag() << "," << sig2.real() << "," << sig2.imag() << "," << sig3.real() << "," << sig3.imag();
		
		
		if (i != grid) mucplot << ",";
	}
	mucplot.close();
}




int main()
{

	k = 1.0;
	p = 4.0;
	
	muval = vector<double>(25);
	for (int i = 0; i <= 24; i++) muval[i] = (0.25*(double)i) - 4.0;
	//muval = vector<double>(31);
	//for (int i = 0; i <= 30; i++) muval[i] = (0.1*(double)i) - 1.0;
	//for (int i = 0; i <= 30; i++) muval[i] = (1.0*(double)i) - 30.0;
	gammaval = vector<double>{-1.0/3.0, -0.3, -0.2, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);
	
	
	//z1_plot(1, 18, 100000);
	
	//mu_plot(5, 100000);
	
	
	//for (int mi = 0; mi < 25; mi++) fiveplot4(1000000, mi);
	//for (int gi = 0; gi < 6; gi++) bunin4(1000000, gi);
	fiveplot4(100000, 16);
	
	
	

	
	//cout << crit_z1() << ", " << phi(crit_z1()) << endl;
	
	
	
	
	return 0;
}
/*
int main()
{
	//info = true;
	p = 4.0;
	mu = 0.0001;
	k = 1.0;
	gama = 0.0;
	//double z1 = -1.3259;
	double z1 = -1000.0;
	
	//cout << "sig " << sigma(-0.8449) <<", " << sigma(-1.3259) << endl;
	//double z1 = -1000.0;
	/*cout << "phi " << phi(z1) << endl;
	cout << "S " << second(z1) << endl;
	cout << "X " << X(z1) << endl;
	cout << "H " << help(z1) << endl;
	cout << "F " << first(z1) << endl;
	cout << "M " << M(z1) << endl;
	cout << "q " << q(z1) << endl;
	cout << "sig " << sigma(z1) << endl;
	cout << "mud " << mud4(z1) << endl;
	
	ofstream file; file.open("M4sig0.txt");
	
	int j = 0;
	
	for (int i = 0; i <= 1001; i++)
	{
		if (i == 747) {mu = 0.0; j = 1;}
		else mu = 0.6*(double)(i-j)/1000.0 - 0.3;
		cout << i << endl;
		double Q = z1*help(z1)/(3.0*mu*first(z1));
		double R = -k/(2.0*mu);
		double D = pow(Q, 3.0) + pow(R, 2.0);
		complex<double> root;
		if (D >= 0.0) root = pow(D, 0.5);
		else root = pow(-D, 0.5)*1i;
		complex<double> S = crt(R + root);
		complex<double> T = crt(R - root);
		complex<double> Trotate = exp(2i*M_PI/3.0);
		complex<double> M1 = S + T;
		complex<double> M2 = -0.5*(S + T) + 0.5i*pow(3, 0.5)*(S - T);
		complex<double> M3 = -0.5*(S + T) - 0.5i*pow(3, 0.5)*(S - T);
		
		cout << M2.real() << endl;
		
		file << mu << "," << M1.real() << "," << M1.imag() << "," << M2.real() << "," << M2.imag() << "," << M3.real() << "," << M3.imag() << endl;
		
	}
	
	
	return 0;
}


*/


