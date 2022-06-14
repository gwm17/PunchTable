#include "MassLookup.h"
#include "PunchTable.h"
#include "GenerateTable.h"
#include "catima/gwm_integrators.h"
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

int main(int argc, char** argv)
{

	std::string options = "";
	if(argc == 2)
	{
		options = argv[1];
	}

	PunchTable::TableParameters params;

	params.projectileA = 1;
	params.projectileZ = 1;

	params.minKE = 8.23;
	params.maxKE = 25.0;
	params.stepKE = 0.01;

	params.minTheta = 0.0;
	params.maxTheta = 75.0;
	params.stepTheta = 0.1;

	double thicknessSABRE = 500.0 * 1.0e-4 * 2.3216 * 1.0e6; //500 um thick times density of Si-> ug/cm^2

	params.targetZ = {{14}};
	params.targetA = {{28}};
	params.targetS = {{1}};
	params.targetThickness = {thicknessSABRE};

	params.filename = "tables/test.ptab";

	if(options == "--make-table" || options == "")
	{
		PunchTable::GenerateTable(params);
	}

	if(options == "--test" || options == "")
	{
		std::cout<<"-------------Testing---------"<<std::endl;
		PunchTable::PunchTable table(params.filename);
		std::cout<<"Is the table valid? "<<table.IsValid()<<std::endl;
		double ke_test = 14.5; //MeV
		double theta_test = 35.5*M_PI/180.0;
		double ke_dep = PunchTable::GetEnergyDeposited(params, theta_test, ke_test);
		double recov_ke = table.GetInitialKineticEnergy(theta_test, ke_dep);
		std::cout<<"For a "<<PunchTable::MassLookup::GetInstance().FindSymbol(params.projectileZ, params.projectileA)<<" with kinetic energy "<<ke_test<<" MeV "
				 <<" the percent error on recovery is "<<(ke_test-recov_ke)/ke_test*100.0<<" where the recovered energy is "<<recov_ke<<" MeV"
				 <<" and the deposited energy is "<<ke_dep<<" MeV"<<std::endl;
		std::cout<<"-----------------------------"<<std::endl;
	}
}