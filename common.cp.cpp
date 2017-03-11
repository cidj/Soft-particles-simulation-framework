#include "common.h"
#include "CompositeSphere.h"
#include <stdlib.h>           
using namespace std;

int stage;
Vector G, centre;
double Time, timestep;
int nstep, nprint, nenergy;
int stepth;
double lx, ly, x_0, y_0;
vector<Material>mat;
vector<CompositeSphere> cell;
unsigned int no_of_cells;
double Erad, esb, etamed, etasub;
ofstream fphase("phase.dat"), flastinit("stagenew.init"), fenergy("energy.dat");
void init_system(char * fname);
double total_kinetic_energy();
int main(int argc, char ** argv)
{
	if (argc != 2){
		cerr << "usage: " << argv[0] << " cell_initfile\n";
		exit(0);
	}
	fenergy.precision(10);
	init_system(argv[1]);
	init_algorithm();
	phase_plot(fphase);
	for (int i = 0; i < nstep; i++)
	{
		stepth = i;
		step();
		if ((i + 1) % nprint == 0){
			cout << "t: " << Time << " s\n";
			phase_plot(fphase);
		}
		if ((i + 1) % nenergy == 0){
			fenergy << Time << "\t" << total_kinetic_energy() << endl;
		}
	}
	stage1_init(flastinit);
	cout << "The MD simulation completed successfully!" << endl;
	system("pause");
	return 0;
}
void integrate()
{
	for (unsigned int i = 0; i < cell.size(); i++){
		if (cell[i].ptype() == 0) {			
			cell[i].nextrtd0(timestep);
		}
		else {
			cell[i].boundary_conditions(i, timestep, Time);
		}
		cell[i].periodic_bc(x_0, y_0, lx, ly);
		cell[i].set_force_and_area();
	}
	make_forces();
	for (unsigned int i = 0; i < cell.size(); i++){
		if (cell[i].ptype() == 0)
		{
			cell[i].nextrtd1(timestep);
		}
	}
	Time += timestep;
}
void init_system(char * fname)
{
	ifstream fcell(fname);
	while (fcell.peek() == '#'){
		string type;
		fcell >> type;
		if (type == "#stage:"){
			fcell >> stage;
			fcell.ignore(100, '\n');
			cout << "stage: " << stage << endl;
		}
		else if (type == "#G:"){
			fcell >> G.x() >> G.y() >> G.phi();
			fcell.ignore(100, '\n');
			cout << "G: " << G << endl;
		}
		else if (type == "#Time:"){
			fcell >> Time;
			fcell.ignore(100, '\n');
			cout << "Time: " << Time << endl;
		}
		else if (type == "#nstep:"){
			fcell >> nstep;
			fcell.ignore(100, '\n');
			cout << "nstep: " << nstep << endl;
		}
		else if (type == "#timestep:"){
			fcell >> timestep;
			fcell.ignore(100, '\n');
			cout << "timestep: " << timestep << endl;
		}
		else if (type == "#nprint:"){
			fcell >> nprint;
			fcell.ignore(100, '\n');
			cout << "nprint: " << nprint << endl;
		}
		else if (type == "#nenergy:"){
			fcell >> nenergy;
			fcell.ignore(100, '\n');
			cout << "nenergy: " << nenergy << endl;
		}
		else if (type == "#lx:"){
			fcell >> lx;
			fcell.ignore(100, '\n');
			cout << "lx: " << lx << endl;
		}
		else if (type == "#ly:"){
			fcell >> ly;
			fcell.ignore(100, '\n');
			cout << "ly: " << ly << endl;
		}
		else if (type == "#x_0:"){
			fcell >> x_0;
			fcell.ignore(100, '\n');
			cout << "x_0: " << x_0 << endl;
		}
		else if (type == "#y_0:"){
			fcell >> y_0;
			fcell.ignore(100, '\n');
			cout << "y_0: " << y_0 << endl;
		}
		else if (type == "#Erad:"){
			fcell >> Erad;
			fcell.ignore(100, '\n');
			cout << "Erad: " << Erad << endl;
		}
		else if (type == "#esb:"){
			fcell >> esb;
			fcell.ignore(100, '\n');
			cout << "esb: " << esb << endl;
		}
		else if (type == "#etamed:"){
			fcell >> etamed;
			fcell.ignore(100, '\n');
			cout << "etamed: " << etamed << endl;
		}
		else if (type == "#etasub:"){
			fcell >> etasub;
			fcell.ignore(100, '\n');
			cout << "etasub: " << etasub << endl;
		}
		else if (type == "#material:"){
			Material mm;
			fcell >> mm;
			if (fcell){
				mat.push_back(mm);
			}
			fcell.ignore(100, '\n');
			cout << "material " << mat.size() - 1 << ": \n"
				<<"rou\tY\tnu\teta1\teta2\tdsa\tsigma\trpl"
				<<"\n" << mm;
		}
		else {
			cerr << "Init: unknown global property: " << type << endl;
			abort();
		}
	}
	while (fcell){
		CompositeSphere pp;
		fcell >> pp;
		if (fcell){
			cell.push_back(pp);
		}
	}
	//centre = Vector(x_0 + 0.5*lx, y_0 + 0.5*ly, 0);
	centre = Vector(0.5e-5, 0.5e-5, 0);
	no_of_cells = cell.size();
	cout << no_of_cells << " cells read\n" << flush;
}
double total_kinetic_energy()
{
	double sum = 0;
	for (unsigned int i = 0; i < cell.size(); i++){
			sum += cell[i].kinetic_energy();
	}
	return sum;
}
void phase_plot(ostream & os)
{
	os << "#Time: " << "\t" << Time << endl;
	//  os << "#no_of_cells: " << "\t" << no_of_cells << endl;
	for (unsigned int i = 0; i < cell.size(); i++)
	{
		os << cell[i];
	}
	os << flush;
}
void stage1_init(ostream & os)
{
	os << "#stage: " << stage + 1 << "\n";
	os << "#G: 0 0 0\n";
	os << "#Time: " << 0 << "\n";

	os << "#nstep: 8640000\n";
	os << "#timestep: " << timestep << "\n";
	os << "#nprint: " << 30000 << "\n";
	os << "#nenergy: " << nenergy << "\n";

	os << "#lx: " << lx << "\n";
	os << "#ly: " << ly << "\n";
	os << "#x_0: " << x_0 << "\n";
	os << "#y_0: " << y_0 << "\n";

	//os << "#lx: " << 2 * lx << "\n";
	//os << "#ly: " << 2 * ly << "\n";
	//os << "#x_0: " << x_0 - 0.5*lx << "\n";
	//os << "#y_0: " << y_0 - 0.5*ly << "\n";

	os << "#Erad: " << 0 << "\n";
	os << "#esb: " << esb << "\n";
	os << "#etamed: " << etamed << "\n";
	os << "#etasub: " << etasub << "\n";
	
	for (unsigned int i = 0; i < mat.size(); i++)
	{
		os << "#material: " << mat[i];
	}

	for (unsigned int i = 0; i < cell.size(); i++)
	{
		//if (cell[i].ptype() != 0)
		//{
		//	cell[i].set_ptype0();// Set all cells mobilable if there are some cells fixed in the previous simulation stage.
		//}
		os << cell[i];
	}
	os << flush;
}