#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iomanip>
#include <stdexcept>

#include <boost/numeric/odeint.hpp>
// See https://stackoverflow.com/questions/22782224/using-boostodeint-with-eigenmatrix-as-state-vector#comment109042140_27781777
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using namespace Eigen;
using namespace std;
using namespace boost::numeric::odeint;
using namespace std::placeholders;

// Parameters
bool output_solutions;
bool output_params;
string ode_method;
double ode_abs_err;
double ode_rel_err;

random_device generator;
mt19937 twist(generator());
//IOFormat p(20, DontAlignCols, ",", ",", "", "", "", "");
IOFormat c(6, DontAlignCols, ",", ",", "", "", "", "");
normal_distribution<double> distributionn(0.0, 1.0);
uniform_real_distribution<double> distributionr(0.0, 1.0);

const int N = 200;
double Nd; // the double version of N
double evolution_time; // The coordinate time to evolve for
double dt; // The timestep size
// these have index starting at 2
ArrayXd mu(3);
ArrayXd gama(3);
ArrayXd sigma(3);
//int p;
double statcheck;
Tensor<double, 2> A;
Tensor<double, 3> B;
Tensor<double, 4> C;
vector<string> code;
vector<double> mu2val;
vector<double> gamma2val;
vector<double> sigma2val;
vector<double> sigma2maxval;
vector<double> mu3val;
vector<double> gamma3val;
vector<double> sigma3val;
vector<double> sigma3maxval;
int v1;
int v2;
int v3;
int r;

string carryon;
bool coutall;
ofstream coutfile;

template<typename T>
using  ArrayType = Eigen::Array<T, Eigen::Dynamic, 1>;

template<typename Scalar, typename sizeType>
auto Tarr(const Eigen::Tensor<Scalar, 1> &tensor,const sizeType rows) // this function maps tensor to array
{
    return Eigen::Map<const ArrayType<Scalar>> (tensor.data(), rows);
}

typedef Eigen::IndexPair<int> IP;
Eigen::array<IP, 1> d10 = { IP(1, 0) }; // second index of T

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar, typename sizeType>
auto Tmat(const Eigen::Tensor<Scalar, 2> &tensor,const sizeType rows, const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows, cols);
}

template<typename Scalar, typename... Dims>
auto Maten(const MatrixType<Scalar> &matrix, Dims... dims)
{
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}

typedef Array<double, N, 1> state_type; // State type for Boost odeint




MatrixXd cholesky(int p) // now corrected for dividing by zero
{
	MatrixXd low = MatrixXd::Zero(p, p);

	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double sum = 0.0;
			if (j == i)
			{
				for (int k = 0; k < j; k++) sum += pow(low(j, k), 2.0);
				low(j, j) = sqrt(1 - sum);
			}
			else
			{
				for (int k = 0; k < j; k++) sum += (low(i, k)*low(j,k));
				if (gama(p-2) - sum == 0.0) low(i, j) = 0.0;
				else low(i, j) = (gama(p-2) - sum)/low(j, j);
				 
			}
		}
	}

	return low;
}


Tensor<double, 2> order2()
{
	Tensor<double, 2> A(N, N);
	A.setZero();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			MatrixXd cho = cholesky(2);
			VectorXd z(2);
			for (int zi = 0; zi < 2; zi++) z(zi) = distributionn(twist);
			MatrixXd list = cho*z;
			A(i,j) = (mu(0)/(Nd - 1.0)) + ((sigma(0)*list(0))/sqrt(Nd - 1.0));
			A(j,i) = (mu(0)/(Nd - 1.0)) + ((sigma(0)*list(1))/sqrt(Nd - 1.0));
		}
	}
	return A;
}


Tensor<double, 3> order3() // now corrected for gamma = 1
{
	Tensor<double, 3> A(N, N, N);
	A.setZero();

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			for (int k = 0; k < j; k++)
			{
				MatrixXd cho(3, 3);
				if (gama(1) == 1.0) {cho.setZero(); cho(0, 0) = 1.0; cho(1, 0) = 1.0; cho(2, 0) = 1.0;}
				else cho = cholesky(3);
				VectorXd z(3);
				for (int zi = 0; zi < 3; zi++) z(zi) = distributionn(twist);
				MatrixXd list = cho*z;
				A(i,j,k) = (2.0*mu(1)/((Nd - 1.0)*(Nd - 2.0))) + (sqrt(3.0)*sigma(1)*list(0)/sqrt((Nd - 1.0)*(Nd - 2.0)));
				A(j,k,i) = (2.0*mu(1)/((Nd - 1.0)*(Nd - 2.0))) + (sqrt(3.0)*sigma(1)*list(1)/sqrt((Nd - 1.0)*(Nd - 2.0)));
				A(k,i,j) = (2.0*mu(1)/((Nd - 1.0)*(Nd - 2.0))) + (sqrt(3.0)*sigma(1)*list(2)/sqrt((Nd - 1.0)*(Nd - 2.0)));
			}
		}
	}
	return A;
}

Tensor<double, 4> order4()
{
	Tensor<double, 4> A(N, N, N, N);
	A.setZero();

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			for (int k = 0; k < j; k++)
			{
				for (int l = 0; l < k; l++)
				{
					MatrixXd cho = cholesky(4);
					VectorXd z(4);
					for (int zi = 0; zi < 4; zi++) z(zi) = distributionn(twist);
					MatrixXd list = cho*z;
					A(i,j,k,l) = (6.0*mu(2)/((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0))) + (sqrt(12.0)*sigma(2)*list(0)/sqrt((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0)));
					A(j,k,l,i) = (6.0*mu(2)/((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0))) + (sqrt(12.0)*sigma(2)*list(1)/sqrt((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0)));
					A(k,l,i,j) = (6.0*mu(2)/((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0))) + (sqrt(12.0)*sigma(2)*list(2)/sqrt((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0)));
					A(l,i,j,k) = (6.0*mu(2)/((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0))) + (sqrt(12.0)*sigma(2)*list(3)/sqrt((Nd - 1.0)*(Nd - 2.0)*(Nd - 3.0)));
				}
			}
		}
	}
	return A;
}
	



struct output_state_and_time
{
	ofstream *trajtfile;
	vector<ArrayXd> *traj;
	vector<double> *times;
	state_type *stat;
	double *stime;
	
	output_state_and_time(ofstream *trajtfile_, vector<ArrayXd> *traj_, vector<double> *times_, state_type *stat_, double *stime_)
	{
		trajtfile = trajtfile_;
		traj = traj_;
		times = times_;
		stat = stat_;
		stime = stime_;
	}
		  
	void operator()(const state_type &x, double t)
	{
		if (coutall)
		{
			//cout << t << ", " << x[0] << endl;
			//cout << t << ", min " << x.minCoeff() << ", max " << x.maxCoeff() << endl;
			coutfile << t << ", min " << x.minCoeff() << ", max " << x.maxCoeff() << endl;
			//cout << "stime " << *stime << endl;
		}
		if (x.maxCoeff() > 1.0e+6 || x.minCoeff() < -0.1) throw string("diverged");
		//if (x.maxCoeff() > 1.0e+6 || !(x.mean() >= 0.0)) throw string("diverged");
		if ((abs(x - *stat) < 0.01).minCoeff() == 0) // if change reset values
		{
			//if (coutall) cout << "changed from time " << *stime << endl;
			*stat = x;
			*stime = t;
		}
		//else if (coutall) cout << "state did not change " << endl;
		
		
		(*trajtfile) << t << "," << x.format(c) << endl;
		(*times).push_back(t);
		(*traj).push_back(x);
		
		if (t - *stime > statcheck) throw string("fixed");
	}
};


class simulation
{
	public:
	state_type xstate;
	Tensor<double, 1> k;
	vector<ArrayXd> *traj; //
	vector<double> *times;
	state_type *stat;
	double stime;
	ofstream *trajtfile;
	
	simulation(vector<ArrayXd> *traj_, vector<double> *times_, state_type *stat_, ofstream *trajtfile_)
	{
		traj = traj_;
		times = times_;
		stat = stat_;
		trajtfile = trajtfile_;
		k = Tensor<double, 1>(N);
		for (int i = 0; i < N; i++) {xstate(i) = distributionr(twist); k(i) = 1.0;}
		*stat = xstate;
		stime = 0.0;


	}
	~simulation() {}

	static void rhs( const state_type &x_array , state_type &dxdt , const double  t, Tensor<double, 2> &A, Tensor<double, 3> &B, Tensor<double, 1> k)
	{
		auto x(TensorMap<Tensor<const double, 1>>(x_array.data(), N));
		Tensor<double, 1> dxdt_tensor = x * (k - x + A.contract(x, d10) + B.contract(x, d10).contract(x, d10)); // change for different orders
		dxdt = Map<Array<double,N,1>>(dxdt_tensor.data());
		if (coutall) coutfile << "max dxdt " << dxdt.maxCoeff() << endl;
		//if (dxdt.maxCoeff() > 100.0 && t < 10.0) throw string("diverged");
	}

	void run(double start, double end)
	{
		if (ode_method == "rk4fixed_boost")
		{
			// Constant timestep RK4
			//state_type xInit = Tarr(x,N);
			runge_kutta4< state_type,double,state_type,double,boost::numeric::odeint::vector_space_algebra > stepper;
			integrate_const(stepper, bind(rhs, _1, _2, _3, A, B, k), xstate, start, end, dt, output_state_and_time(trajtfile, traj, times, stat, &stime)); // change
		}
		else if (ode_method == "adaptive_boost")
		{
			// Default adaptive method
			//state_type xInit = Tarr(x,N);
			integrate(bind(rhs, _1, _2, _3, A, B, k), xstate, start, end, dt, output_state_and_time(trajtfile, traj, times, stat, &stime)); // change
		}
		else if (ode_method == "adaptive_boost_error_controlled")
		{
			// Adaptive RK with error controls
			typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
			//state_type xInit = Tarr(x,N);
			size_t steps = integrate_adaptive(
				make_controlled< error_stepper_type >(ode_abs_err, ode_rel_err),
				bind(rhs, _1, _2, _3, A, B, k), xstate, start, end, dt, // change
				output_state_and_time(trajtfile, traj, times, stat, &stime));
		}
		else
		{
			cerr << "Unrecognised ode_method '" << ode_method << "'" << endl;
			abort();
		}	
	}
};

void datarun()
{
	ofstream trajxtfile;
	ofstream trajytfile;
	ofstream paramsfile;
	ofstream datafile;
	
	coutfile.open("cout23_" + to_string(v1) + "_" + to_string(v2) + "_" + to_string(v3) + "_" + to_string(r) + ".txt");

	if (output_solutions)
	{
		trajxtfile.open("trajxt23_" + to_string(v1) + "_" + to_string(v2) + "_" + to_string(v3) + "_" + to_string(r) + ".txt");
		trajytfile.open("trajyt23_" + to_string(v1) + "_" + to_string(v2) + "_" + to_string(v3) + "_" + to_string(r) + ".txt");		
	}
	if (output_params)
	{
		paramsfile.open("params23_" + to_string(v1) + "_" + to_string(v2) + "_" + to_string(v3) + "_" + to_string(r) + ".yaml");
		paramsfile << "dt: " << dt << endl;
		paramsfile << "N: " << N << endl;
		paramsfile << "mu: " << mu.format(c) << endl;
		paramsfile << "gamma: " << gama.format(c) << endl;
		paramsfile << "sigma: " << sigma.format(c) << endl;
		paramsfile << "ode_method: " << ode_method << endl;
		paramsfile << "ode_abs_err: " << ode_abs_err << endl;
		paramsfile << "ode_rel_err: " << ode_rel_err << endl;
	}
	
			
	bool diverged = false;
	bool unique = false;
	bool xfixed = false;
	bool yfixed = false;
		
	vector<ArrayXd> trajx;
	vector<ArrayXd> trajy;
	vector<double> timesx;
	vector<double> timesy;
	
	state_type xstat;
	state_type ystat;
	
	double start = 0.0;
	
	simulation simx(&trajx, &timesx, &xstat, &trajxtfile);
	try {simx.run(start, evolution_time);}
	catch (string stopped)
	{
		if (stopped == "diverged") {diverged = true; coutfile << "divergedx" << endl;}
		else if (stopped == "fixed") {xfixed = true; coutfile << "fixedx" << endl;}
		else coutfile << "bad catch" << endl;
	}
		
	if (xfixed)
	{
		coutfile << "running y sim" << endl;
		simulation simy(&trajy, &timesy, &ystat, &trajytfile);
		try {simy.run(start, evolution_time);}
		catch (string stopped)
		{
			if (stopped == "diverged") {diverged = true; coutfile << "divergedy" << endl;}
			else if (stopped == "fixed") {yfixed = true; coutfile << "fixedy" << endl;}
			else coutfile << "bad catch" << endl;
		}
	}
	
	datafile.open("data23_" + to_string(v1) + "_" + to_string(v2) + "_" + to_string(v3) + "_" + to_string(r) + ".txt");
	if (diverged == false)
	{
		ArrayXd measures = ArrayXd::Zero(3); // <phi>, <M>, <q>
		ArrayXd sum = ArrayXd::Zero(N);
		double stop = timesx.back();
		int tcount = 0;
		for (int t = 0; t < timesx.size(); t++)
		{
			if (stop - timesx[t] < 500.0)
			{
				tcount ++;
				
				double M = trajx[t].mean();
				double q = (trajx[t]*trajx[t]).mean();
				for (int i = 0; i < N; i++) if (trajx[t](i) > 0.0001) measures(0) ++;

				sum += trajx[t];
				measures(1) += M;
				measures(2) += q;
			}
		}
		

		sum /= (double)tcount; // average x(i) over T
		measures(0) /= Nd*(double)tcount; // average phi T
		measures(1) /= (double)tcount; // average M over T
		measures(2) /= (double)tcount; // average q over T

		datafile << measures.format(c);
	
		if (xfixed && yfixed) // check unique
		{
			ArrayXd ysum = ArrayXd::Zero(N);
			double ystop = timesy.back();
			int ytcount = 0;
			for (int t = 0; t < timesy.size(); t++)
			{
				if (ystop - timesy[t] < 500.0)
				{
					ytcount ++;
					ysum += trajy[t];
				}
			}
			ysum /= (double)ytcount; // average y(i) over T
			if ((abs(sum - ysum) < 0.01).minCoeff() == 1) {unique = true; coutfile << "unique" << endl;}

			if (unique) datafile << ",1"; 
			else datafile << ",2";
		}
		else datafile << ",3";
	}
	else datafile << "0,0,0,0";
	
	trajxtfile.close();
	trajytfile.close();
	paramsfile.close();
	datafile.close();
}


int main(int argc, char** argv) // p = 23
{
	gama = ArrayXd::Zero(3);
	mu = ArrayXd::Zero(3);
	sigma = ArrayXd::Zero(3);
	//p = 3;
	// change matrix in rhs()
	// change matrix in ode methods
	// N = 200;
	Nd = (double)N;
	evolution_time = 10000.0; //1000;
	dt = 0.1; //0.01;
	statcheck = 500.0;
	output_solutions = true;
	output_params = false;
	ode_method = "adaptive_boost_error_controlled";
	ode_abs_err = 1.0e-4;
	ode_rel_err = 0;
	
	code = {"A", "B", "C", "D", "E", "G", "N", "O"};
	mu2val = vector<double>{-2.0, -2.0, -3.5, -4, -3.37, -2.0, -3.38, -3.38};
	mu3val = vector<double>{4.0, -4.0, 6.0, 5.0, 6.0, -4.0, 6.0, 6.0};
	gamma2val = vector<double>{-1.0, -1.0, -0.8, -0.8, -0.8, 0.0, -0.88, -0.79};
	gamma3val = vector<double>{-0.5, -0.5, -0.4, -0.4, -0.4, 0, -0.39, -0.44};
	sigma2maxval = vector<double>{5.0, 5.0, 4.0, 8.0, 3.0, 6.0, 6.0, 3.0};
	sigma3maxval = vector<double>{5.0, 10.0, 5.0, 10.0, 4.0, 3.0, 4.0, 8.0};
	
	v1 = atoi(argv[1]); // code 0-7
	v2 = atoi(argv[2]); // s2i 0-20
	v3 = atoi(argv[3]); // s3i 0-20
	r = atoi(argv[4]); // ri 0-19
	
	coutall = true;
	
	twist.seed((v1+1)*(v2+1)*(v3+1)*(r+1));
	
	double sig2max = sigma2maxval[v1];
	sigma2val = vector<double>(21);
	for (int i = 0; i <= 20; i ++) sigma2val[i] = sig2max*(double)i/20.0;
	
	double sig3max = sigma3maxval[v1];
	sigma3val = vector<double>(21);
	for (int i = 0; i <= 20; i ++) sigma3val[i] = sig3max*(double)i/20.0;
	
	mu(0) = mu2val[v1];
	gama(0) = gamma2val[v1];
	sigma(0) = sigma2val[v2];
	
	mu(1) = mu3val[v1];
	gama(1) = gamma3val[v1];
	sigma(1) = sigma3val[v3];
	
	A = order2();
	B = order3();
	
	datarun();

	return 0;
}


/*
int main(int argc, char** argv) // p = 2
{
	gama = ArrayXd::Zero(3);
	mu = ArrayXd::Zero(3);
	sigma = ArrayXd::Zero(3);
	p = 2;
	// change matrix in rhs()
	// change matrix in ode methods
	// N = 200;
	Nd = (double)N;
	evolution_time = 10000.0; //1000;
	dt = 1.0e-6; //0.01;
	statcheck = 500.0;
	output_solutions = true;
	output_params = true;
	ode_method = "adaptive_boost_error_controlled";
	ode_abs_err = 1.0e-4;
	ode_rel_err = 0;
	
	muval = vector<double>(25);
	for (int i = 0; i <= 24; i++) muval[i] = (0.25*(double)i) - 4.0;
	gammaval = vector<double>{-1.0, -0.8, -0.5, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);
	
	coutall = false;
	
	
	v1 = atoi(argv[1]);
	v2 = atoi(argv[2]);
	v3 = atoi(argv[3]);
	r = atoi(argv[4]);
	
	twist.seed((v1+1)*(v2+1)*(v3+1)*(r+1));
	
	mu(p-2) = muval[v1];
	gama(p-2) = gammaval[v2];
	sigma(p-2) = sigmaval[v3];
	
	A = order2();
	
	datarun();

	return 0;

}*/
/*
int main(int argc, char** argv) // p = 3
{
	gama = ArrayXd::Zero(3);
	mu = ArrayXd::Zero(3);
	sigma = ArrayXd::Zero(3);
	p = 3;
	// change matrix in rhs()
	// change matrix in ode methods
	// N = 200;
	Nd = (double)N;
	evolution_time = 10000.0; //1000;
	dt = 0.1; //0.01;
	statcheck = 500.0;
	output_solutions = false;
	output_params = false;
	ode_method = "adaptive_boost_error_controlled";
	ode_abs_err = 1.0e-4;
	ode_rel_err = 0;
	
	muval = vector<double>(25);
	for (int i = 0; i <= 24; i++) muval[i] = (0.25*(double)i) - 4.0;
	gammaval = vector<double>{-0.5, -0.4, -0.3, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);
	
	v1 = atoi(argv[1]);
	v2 = atoi(argv[2]);
	v3 = atoi(argv[3]);
	r = atoi(argv[4]);
	
	coutall = true;
	
	twist.seed((v1+1)*(v2+1)*(v3+1)*(r+1));
	
	mu(p-2) = muval[v1];
	gama(p-2) = gammaval[v2];
	sigma(p-2) = sigmaval[v3];
	
	B = order3();
	
	datarun();

	return 0;
}
*/				
/*
int main(int argc, char** argv) // p = 4
{
	gama = ArrayXd::Zero(3);
	mu = ArrayXd::Zero(3);
	sigma = ArrayXd::Zero(3);
	p = 4;
	// change matrix in rhs()
	// change matrix in ode methods
	// N = 200;
	Nd = (double)N;
	evolution_time = 2000.0; //1000;
	dt = 0.001; //0.01;
	statcheck = 500.0;
	output_solutions = true;
	output_params = true;
	ode_method = "adaptive_boost_error_controlled";
	ode_abs_err = 1.0e-4;
	ode_rel_err = 0;

	muval = vector<double>(25);
	for (int i = 0; i <= 24; i++) muval[i] = (0.25*(double)i) - 4.0;
	gammaval = vector<double>{-1.0/3.0, -0.3, -0.2, 0.0, 0.5, 1.0};
	sigmaval = vector<double>(21);
	for (int i = 0; i <= 20; i++) sigmaval[i] = pow(10.0, (0.1*(double)i) - 1.0);

	v1 = atoi(argv[1]);
	v2 = atoi(argv[2]);
	v3 = atoi(argv[3]);
	r = atoi(argv[4]);
	
	twist.seed((v1+1)*(v2+1)*(v3+1)*(r+1));
	
	mu(p-2) = muval[v1];
	gama(p-2) = gammaval[v2];
	sigma(p-2) = sigmaval[v3];
	
	C = order3();
	
	datarun();

	return 0;
}*/





