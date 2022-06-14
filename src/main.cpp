#include "MassLookup.h"
#include "PunchTable.h"
#include "ElossTable.h"
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

	bool test = true;
	if(argc != 2)
	{
		std::cerr<<"PunchTable requires an input file!"<<std::endl;
		return 1;
	}

	PunchTable::TableParameters params = PunchTable::GetTableParameters(argv[1]);

	if(params.tableType == "PunchTable")
	{
		PunchTable::GeneratePunchTable(params);
	}
	else if(params.tableType == "ElossTable")
	{
		PunchTable::GenerateElossTable(params);
	}
	else
	{
		std::cerr<<"Unrecognized table type "<<params.tableType<<std::endl;
		return 1;
	}

	if(test && params.tableType == "PunchTable")
	{
		std::cout<<std::endl;
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
	if(test && params.tableType == "ElossTable")
	{
		std::cout<<std::endl;
		std::cout<<"-------------Testing---------"<<std::endl;
		PunchTable::ElossTable table(params.filename);
		std::cout<<"Is the table valid? "<<table.IsValid()<<std::endl;
		double ke_test = 14.5; //MeV
		double theta_test = 35.5*M_PI/180.0;
		double eloss = table.GetEnergyLoss(theta_test, ke_test);
		std::cout<<"For a "<<PunchTable::MassLookup::GetInstance().FindSymbol(params.projectileZ, params.projectileA)<<" with kinetic energy "<<ke_test<<" MeV "
				 <<"calculate an energy loss of "<<eloss<<" MeV. Please compare to favorite tool (LISE, SRIM, etc.)"<<std::endl;
		std::cout<<"-----------------------------"<<std::endl;
	}

	return 0;
}