// lotka-volterra dynamics with eigen array class, RK4 instead of euler
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <omp.h>
#include <Eigen/Dense>
#include <iomanip>
using namespace Eigen;
using namespace std;

random_device generator;
mt19937 twist(generator());
IOFormat c(6, DontAlignCols, ",", ",", "", "", "", "");

int runs = 20;
int p = 4;

int main()
{
	
	for (int mi = 0; mi < 25; mi ++)
	{
		//int mi = 11;
		for (int gi = 0; gi < 6; gi ++)
		{
			//int gi = 3;
			for (int si = 0; si < 21; si ++)
			{
				//cout << to_string(mi) << "_" << to_string(gi) << "_" << to_string(si) << endl;
				
				int count1 = 0;
				int count2 = 0;
				int count3 = 0;
				int countbad = 0; // count the number of files not there
	
				ArrayXd measures = ArrayXd::Zero(3);
		
				for (int r = 0; r < runs; r++)
				{
					//cout << r << endl;
					string data;
					ifstream file("data" + to_string(p) + "_" + to_string(mi) + "_" + to_string(gi) + "_" + to_string(si) + "_" + to_string(r) + ".txt");
					if (file.is_open())
					{
						getline(file, data);
						file.close();
					}
					else {data = "0,0,0,0"; countbad ++;}// cout << "no file " << mi << ", " << gi << ", " << si << ", " r << endl;}
					
					if (data.size() == 0) {data = "0,0,0,0"; countbad ++;}
					
					stringstream ss(data);
					ArrayXd meas = ArrayXd::Zero(3);
					int type;
					
					for (int i = 0; i < 4; i ++)
					{
						string s;
						getline(ss, s, ',');
						if (i < 3) meas(i) = stod(s);
						else type = stoi(s);
					}
						
					if (type != 0)
					{
						measures += meas;
						if (type == 1) count1 ++;
						if (type == 2) count2 ++;
						if (type == 3) count3 ++;
					}

					//else for (int i = 0; i < 3; i ++) meas(i) = 1.0/0.0;
					//cout << meas.format(c);
					//if (r != runs-1) cout << ",";
				}
				int countruns = count1 + count2 + count3;
	
				double R = (double)count1/(double)runs; // normalised values
				double G = (double)count2/(double)runs;
				double B = (double)count3/(double)runs;
	
				if (R + G + B == 3.0 || countbad > 0) {R = 0.0; G = 0.0; B = 0.0;}// cout << "countbad " << countbad << " for " << mi << ", " << gi << ", " << si << endl;} // lintrans
				else
				{
					double div = 1.0 - (R + G + B);

					R += div;
					G += div;
					B += div;
				}
				measures /= (double)countruns;
				
				//cout << measures.format(c);
	
				cout << R << "," << G << "," << B;
				
				if (si != 20) cout << ",";
			}
			if (gi != 5) cout << ",";
		}
		if (mi != 24) cout << ",";
	}
	return 0;
}














	
