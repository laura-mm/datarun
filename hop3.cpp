// numerical solutions for p = 3
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
double p = 3.0;
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
	return sec/(sec + (gama*(p - 1.0)*phi(z1)));
}
// this is undefined for sigma = 0

double first(double z1)
{
	return (z1*(erf(z1/sqrt(2.0)) - 1.0)/2.0) + (exp(-z1*z1/2.0)/sqrt(2.0*M_PI));
}
// this is the first order integral, this can not be zero

double M(double z1) //this is a solution to a p-order polynomial
{	
	if (mu == 0.0) return -k*first(z1)/(z1*help(z1));
	else
	{
		return ((z1*help(z1)/first(z1)) + sqrt(pow(z1*help(z1)/first(z1), 2.0) - (4.0*mu*k)))/(-2.0*mu);
	}

}

double M2(double z1) //other solution
{	
	if (mu == 0.0) return -k*first(z1)/(z1*help(z1));
	else
	{
		return ((z1*help(z1)/first(z1)) - sqrt(pow(z1*help(z1)/first(z1), 2.0) - (4.0*mu*k)))/(-2.0*mu);
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

double q2(double z1)
{
	return second(z1)*pow(M2(z1)/first(z1), 2.0);
}

double sigma(double z1)
{
	return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M(z1), p - 2.0);
}

double sigma2(double z1)
{
	return help(z1)*sqrt(2.0/(p*pow(second(z1), p - 1.0)))*pow(first(z1)/M2(z1), p - 2.0);
}

double mud3(double z1) // for p = 3, value of mu where no solutions for M, fixed gama
{
	return pow(z1*help(z1)/first(z1), 2.0)/(4.0*k);
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


double div3_z1(double zstart) // for p = 3, for the set value of mu and gama, finds the value of zd
{
	double z1 = zstart;
	double small = pow(10.0, -15.0);
	double y = mud3(z1) - mu; if (info == true) cout << "y " << y << endl;
	double grad;
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (mud3(z1 + small) - mud3(z1 - small))/(2.0*small);
		if (info == true) cout << "grad " << grad << endl;
		if (grad == 0) {if (abs(y) < pow(10.0, -10.0)) return z1; else return -1.0/0.0;}
		znew = z1 - (y/grad);
		ynew = znew*help(znew)/first(znew) + mu;
		if (info == true) {cout << "znew " << znew << ", ynew " << ynew << endl;}
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -10.0))) break;
		y = ynew;
		z1 = znew;
		if (small > pow(10.0, -20.0)) small /= 10.0;
	}
	cout << "div3 found " << z1 << endl;
	return z1;
}


		
void z1_plot(int gi, int mi, int grid) // plots how things vary with z1
{
	// mu = -4 or -3. gama = 1
	mu = 0.4; //muval[mi];
	gama = -0.4; //gammaval[gi];
	bool flag = false;
	ofstream z1plot; z1plot.open("z1plot3_" + to_string(gi) + "_" + to_string(mi) + ".txt");
		
	for (int i = 0; i <= 79155; i++)
	{
		double z1 = ((double)i*25.0/(double)grid) - 20.0;
		if (i == 79155) z1 = -0.2116;
		//if (!(abs(M(z1)) > 0.0) || !(abs(M2(z1)))) continue;
		if (flag == true) z1plot << ",";
		flag = true;
		z1plot << z1 << "," << phi(z1) << "," << M(z1) << "," << q(z1) << "," << sigma(z1) << "," << M2(z1) << "," << q2(z1) << "," << sigma2(z1);
	}
	z1plot.close();
}		
		
		
		
/*
		ofstream zdpoints; zdpoints.open("zdpoints3_" + to_string(gi) + "_" + to_string(mi) + ".txt");
		// find any divergent points
		double zlow;
		double zhigh;
		double MM = sqrt(k/mu); // above M doesnt work
		if (mu >= 0.25) // look for lower bound
		{
			for (int i = 0; i <= grid; i++)
			{
				double z1 = 100.0*(double)i/(double)grid - 50.0;
				if (sigma(z1) > 0.0 && sigma(z1) < 1.0/0.0) {zlow = z1; break;}
			}
			if (abs(mud3(zlow) - mu) < pow(10.0, -3.0)) zlow = div3_z1(zlow);
			if (abs(mud3(zlow) - mu) < pow(10.0, -10.0)) zdpoints << zlow << "," << phi(zlow) << "," << second(zlow)*pow(MM/first(zlow), 2.0);
		}
		
				
						double MM = sqrt(k/mu); // above M doesnt work
		double sigmaa = help(z1)*sqrt(2.0/(3.0*pow(second(z1), 2.0)))*first(z1)/MM;
		double qq = second(z1)*pow(MM/first(z1), 2.0);
				
				 zdpoints << sigma(zlow); // check divergent condition
				else points << 1.0/0.0;
				if (mi != 24) points << ",";
*/
		


double maxsig_z1(double zstart) // returns the z1 where sigma reaches a stationary point, depends on mu and gamma
{
	double z1 = zstart;
	double small = pow(10.0, -6.0);
	double y = (sigma(z1 + small) - sigma(z1 - small))/(2.0*small);
	double grad;
	if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
	double ynew;
	double znew;
	while (abs(y) > 0.0)
	{
		grad = (sigma(z1 + small) + sigma(z1 - small) - (2.0*sigma(z1)))/(small*small);
		if (grad == 0.0) break;
		if (info == true) {cout << "small " << small << " grad " << grad << endl;}
		znew = z1 - (y/grad);
		ynew = (sigma(znew + small) - sigma(znew - small))/(2.0*small);
		if (abs(ynew) >= abs(y) && ((ynew*y) > 0.0 || abs(y) < pow(10.0, -5.0))) break;
		y = ynew;
		z1 = znew;
		if (info == true) {cout << "z1 " << z1 << ", y " << y << endl;}
		if (small > pow(10.0, -10.0)) small /= 10.0;
	}
	return z1;	
}


void fiveplot3(int grid, int mi) // plots order parameters for varying sigma
{
	mu = muval[mi];
	
	double zcrit = crit_z1();
	
	for (int gi = 0; gi < 6; gi++)
	{
		cout << mi << "," << gi << endl;
		
		string filename = "measfpb3_" + to_string(mi) + "_" + to_string(gi)+ ".txt";
		ofstream file; file.open(filename);

		gama = gammaval[gi];
	
		double sig_prev = -1.0; // next sigma has to be bigger than previous

		for (int i = 0; i <= grid; i++)
		{
			double z1 = ((double)i*100.0/(double)grid) - 50.0;
			if (sigma(z1) < 0.0) continue;
			//if (sigma(z1) < sig_prev) break;
			sig_prev = sigma(z1);
			if (i != 0) file << ",";
			file << sigma(z1) << "," << z1 << "," << phi(z1) << "," << M(z1) << "," << q(z1) << "," << help(z1);
		}
		file.close();
	}
}

void bunin3(int grid, int gi) // plots phase diagram and finds transition points for fiveplots
{
	gama = gammaval[gi];
	double critz1 = crit_z1();
	double zstart = 0.0;
	double zstop;
	//info = true;
	
	//ofstream zdpoints; zdpoints.open("zdpoints3_" + to_string(gi) + ".txt"); // extract particular values of z1

	
	ofstream points; points.open("points3_" + to_string(gi) + ".txt");
	// points has sigc, sigmax, sigdlower
		
	ofstream fil; fil.open("mu3_" + to_string(gi) + ".txt");
	for (int i = 0; i <= grid; i++)
	{
		mu = (10.0*(double)i/(double)grid) - 7.0;
		if (i != 0) fil << ",";
		zstart = maxsig_z1(zstart); // zstart = zm
		//if (info == true) cin >> carryon;
		if (zstart < critz1) fil << mu << "," << 1.0/0.0 << "," << sigma(zstart); // no crit and maxsig line
		else fil << mu << "," << sigma(critz1) << "," << sigma(zstart); // critical line and maxsig line
		if (abs(mu - 0.4) < 1.0e-15) cout << zstart << ", " << critz1 << endl;
		for (int mi = 0; mi < 25; mi ++) // check if suits any mu
		{
			if (abs(mu - muval[mi]) < pow(10.0, -15.0))
			{
				//if (mi == 0) cout << zstart << ", " << critz1 << endl;
				if (zstart < critz1) points << 1.0/0.0 << "," << sigma(zstart) << ","; // dont observe crit
				else points << sigma(critz1) << "," << sigma(zstart) << ","; // observe crit
				double zlow; // find the lowest value of z1 which has solutions
				for (int i = 0; i <= grid; i++)
				{
					double z1 = 100.0*(double)i/(double)grid - 50.0;
					if (sigma(z1) > 0.0 && sigma(z1) < 1.0/0.0) {zlow = z1; break;}
				}
				if (abs(mud3(zlow) - mu) < 0.001) points << sigma(zlow); // check divergent condition
				else points << 1.0/0.0;
				if (mi != 24) points << ",";
			}
		}
				
		if (sigma(zstart) > 0.0) zstop = zstart; // find zm for highest possible mu
	}
	fil.close();
	points.close();
	ofstream file; file.open("bunin3_" + to_string(gi) + ".txt"); // plot zd phase
	for (int i = 0; i <= grid; i++)
	{
		double z1 = (50.0 + zstop)*(double)i/(double)grid - 50.0; // dont need zd past where we have zm
		//double z1 = 60.0*(double)i/(double)grid - 50.0;
		mu = mud3(z1);
		if (abs(mu - 0.4) < 1.0e-5) cout << mu << ", " << z1 << endl;
		double MM = sqrt(k/mu); // above M doesnt work
		double sigmaa = help(z1)*sqrt(2.0/(3.0*pow(second(z1), 2.0)))*first(z1)/MM;
		double qq = second(z1)*pow(MM/first(z1), 2.0);
		if (!(sigmaa >= 0.0)) break;
		if (i != 0) file << ",";
		file << mu << "," << sigmaa;
	}
	file.close();
}		

int main()
{
	k = 1.0;
	
	muval = vector<double>(25);
	for (int i = 0; i <= 24; i++) muval[i] = (0.25*(double)i) - 4.0;
	//muval = vector<double>(31);
	//for (int i = 0; i <= 30; i++) muval[i] = (0.1*(double)i) - 1.0;
	gammaval = vector<double>{-0.5, -0.4, -0.3, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);
	
	//cout << crit_z1() << ", " << phi(crit_z1()) << endl;
	
	for (int i = 0; i < 6; i++) bunin3(100000, i);
	//for (int i = 0; i < 25; i++) fiveplot3(100000, i);
	
	//muval.push_back(0.4);
	//z1_plot(0, 0, 1000000);
	

	return 0;
}

/*
int main()
{
	p = 3.0;
	k = 1.0;
	//gama = -0.4;
	int grid = 100000;
	//double z1 = -1000.0;
	//mu = 0.2;
	
	info = true;
	
	//z1_plot(0, 1, grid);
	
	
	muval = vector<double>(25);
	for (int i = 0; i <= 24; i++) muval[i] = (0.25*(double)i) - 4.0;
	gammaval = vector<double>{-0.5, -0.4, -0.3, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);
	
	//bunin3(1000000, 1);
	
	gama = -0.4;
	mu = 0.4;
	div3_z1(-0.2116);
	//div3_z1(-1.54754);
	
	
	/*
	ofstream file; file.open("M3sig0.txt");
	
	for (int i = 0; i <= 1000; i++)
	{
		mu = 10.0*(double)i/1000.0 - 5.0;
		if (mu == 0.0) file << mu << "," << 1.0/0.0 << "," << 1.0/0.0 << "," << 1.0/0.0 << "," << 1.0/0.0 << endl;
		else
		{
			if (pow(z1*help(z1)/first(z1), 2.0) - (4.0*mu*k) >= 0.0) file << mu << "," << M(z1) << "," << 0.0 << "," << M2(z1) << "," << 0.0 << endl;
			else file << mu << "," << (z1*help(z1)/first(z1))/(-2.0*mu) << "," << sqrt((4.0*mu*k) - pow(z1*help(z1)/first(z1), 2.0))/(-2.0*mu) << "," << (z1*help(z1)/first(z1))/(-2.0*mu) << "," << sqrt((4.0*mu*k) - pow(z1*help(z1)/first(z1), 2.0))/(2.0*mu) << endl;
		}
	}*/
	
	
	
	/*
	for (int i = 0; i <= grid; i++)
	{
		//double z1 = (50.0 + zstop)*(double)i/(double)grid - 50.0; // dont need zd past where we have zm
		double z1 = 60.0*(double)i/(double)grid - 50.0;
		mu = mud3(z1);
		double MM = sqrt(k/mu); // above M doesnt work
		double sigmaa = help(z1)*sqrt(2.0/(3.0*pow(second(z1), 2.0)))*first(z1)/MM;
		double qq = second(z1)*pow(MM/first(z1), 2.0);
		if (!(sigmaa >= 0.0)) break;
		cout << "z1 " << z1 << ", mud " << mu << ", sigma " << sigmaa << endl;
	}
	
	info = true;
	p = 3.0;
	mu = 0.0;
	k = 1.0;
	gama = -0.5;
	double z1 = 0.01;
	cout << "phi " << phi(z1) << endl;
	cout << "S " << second(z1) << endl;
	cout << "X " << X(z1) << endl;
	cout << "H " << help(z1) << endl;
	cout << "F " << first(z1) << endl;
	cout << "M " << M(z1) << endl;
	cout << "q " << q(z1) << endl;
	cout << "sig " << sigma(z1) << endl;
	cout << "mud " << mud3(z1) << endl;
	
	return 0;
}

*/




