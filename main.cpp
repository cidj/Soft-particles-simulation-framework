#define _USE_MATH_DEFINES
#include <math.h>
#include "rand48.h"
#include<iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

void dump_particle(ostream & os,
	double x, double y, double phi,
	double vx, double vy, double vw,
                   double radius, double mass, int type,
       			   double Young_mod, double poission_ratio, double etacel,
				   double dsa, double sigma,
				   double force_x, double force_y, double moment_w, int kind, int index)
{
	os << x << "\t" << y << "\t" << phi << "\t" 
		<< vx << "\t" << vy << "\t"<< vw << "\t" 
		<< radius << "\t" << mass << "\t" << type << "\t"
	 << Young_mod << "\t" << poission_ratio << "\t" << etacel << "\t" 
	 << dsa << "\t" << sigma << "\t"
	 << force_x << "\t" << force_y << "\t" << moment_w << "\t" 
	 << kind << "\t" << index << "\t\n";
}

//unit: etacel*length*velocity=force

int main()
{
    ofstream fout("stage0.init");
	fout << "#stage: 0\n";
	fout << "#G: 0 -1e-5 0\n";
	//fout << "#G: 0 0 0\n";
	fout << "#Time: 0\n";

	fout << "#nstep: 100000\n"; 
	fout << "#timestep: 0.01\n";
	fout << "#nprint: 1000\n";
	fout << "#nenergy: 10000\n";

	fout << "#lx: 170e-6\n";
	fout << "#ly: 170e-6\n";
	fout << "#x_0: 0.0\n";
	fout << "#y_0: 0.0\n";

	fout << "#T: 273.15\n";
	fout << "#etamed: " << 1.8e-4 / 9e-6 << "\n";
	fout << "#etasub: " << 0.1 / 9e-6 << "\n";

	const double Rmax = 10e-6, Rmin = 6e-6, Vmax = 9.0e-8, Vmin = 1.0e-8;
	double angle = 0.0, vw = 0.0;
	double rou = 1.1;
	int type = 0, type1=0;
	double Young_mod = 1800, poission_ratio = 0.48;
	double etacel = 0.1 / 9e-6;
	double dsa = 1.0e15, sigma = 2.0e-4;
	int kind = 0;

	const double Rmax2 = 15e-6, Rmin2 = 9e-6;
	double Young_mod2 = 800, etacel2 = 0.88*etacel, rou2=1.05;
	double dsa2 = 0.75*dsa, sigma2 = 0.5*sigma;
	int kind2 = 1;

	double force_x = 0.0, force_y = 0.0, moment_w = 0.0;
    
	int m = 12, n = 12;

    for(int i=0; i<m; i++){
        for(int k=0; k<n; k++){
			double centerx = 170e-6/12*i;
			double centery = 170e-6/12*k;
            double z=drand48();
            double r=Rmin*Rmax/(Rmax-z*(Rmax-Rmin));
			double m = 4 / 3 * M_PI*r*r*r*rou;
            z = drand48();
			double vx = Vmin*Vmax / (Vmax - z*(Vmax - Vmin));
			z = drand48();
			double vy = Vmin*Vmax / (Vmax - z*(Vmax - Vmin));

			z = drand48();
			double R2 = Rmin2*Rmax2 / (Rmax2 - z*(Rmax2 - Rmin2));

			int index = i*n+k;

			double m2 = 4 / 3 * M_PI*R2*R2*R2*rou2;

			if ((i == 3 || i == 8) && (k == 3 || k == 8))
			{
				dump_particle(fout, centerx, centery, angle, vx, vy, vw, R2, m2, type,
					Young_mod2, poission_ratio, etacel2, dsa2, sigma2, force_x, force_y, moment_w, kind2,index);
			}
			else if ((i == 0||k == 0))
			{
				dump_particle(fout, centerx, centery, angle, vx, vy, vw, r, m, type1,
					Young_mod, poission_ratio, etacel, dsa, sigma, force_x, force_y, moment_w, kind,index);
			}
			else
			{
				dump_particle(fout, centerx, centery, angle, vx, vy, vw, r, m, type,
					Young_mod, poission_ratio, etacel, dsa, sigma, force_x, force_y, moment_w, kind,index);
			}
        }
    }
	cout << "The initiation completed successfully!" << endl;
	//system("pause");
    return 0;
}
