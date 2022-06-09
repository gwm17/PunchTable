#include "MassLookup.h"
#include "PunchTable.h"
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

	int Zp = 1;
	int Ap = 1;

	float keStep = 0.01;
	float thickness = 500.0 * 1e-4 * 2.3216 * 1e6; //1000 um thick times density of Si-> ug/cm^2
	float deg2rad = M_PI/180.0;
	float thetaStep = 0.1*deg2rad;
	float precision = 1.0e-6;

	std::vector<int> ZT(1), AT(1), ST(1);
	ZT[0] = 14;
	AT[0] = 28;
	ST[0] = 1;

	EnergyLoss energyLoss;
	energyLoss.SetTargetComponents(ZT, AT, ST);

	float keMin = 8.23f;
	float keMax = 25.0f;

	float thetaMin = 0.0f*deg2rad;
	float thetaMax = 85.0f*deg2rad;

	int nebins = int((keMax-keMin)/keStep);
	int ntbins = int ((thetaMax-thetaMin)/thetaStep);

	if(options == "--make-table" || options == "")
	{

		std::cout<<"nebins: "<<nebins<<" ntbins: "<<ntbins<<std::endl;
	
		float keCur;
		float thetaCur;
		float energy;
	
		std::vector<std::vector<float>> energy_dep(ntbins);
		std::vector<std::vector<float>> energy_in(ntbins);
	
		for(int i=0; i<ntbins; i++)
		{
			std::cout<<"\rWorking on theta bin "<<i+1<<" of "<<ntbins<<std::flush;
			thetaCur = thetaMin + i*thetaStep;
			for(int j=0; j<nebins; j++)
			{
				keCur = keMax - j*keStep;
				energy = energyLoss.GetEnergyLoss(Zp, Ap, keCur, thickness / std::fabs(std::cos(thetaCur)));
				if(std::fabs(energy-keCur) > precision)
				{
					energy_dep[i].push_back(energy);
					energy_in[i].push_back(keCur);
				}
			}
		}
		std::cout<<std::endl;
	
		std::ofstream output("tables/test.ptab");
		output<<std::setprecision(5);
		output<<"Material: ";
		for(size_t i=0; i<ZT.size(); i++)
			output<<MassLookup::GetInstance()->FindSymbol(ZT[i], AT[i]);
		output<<std::endl;
		output<<"Thickness: "<<thickness<<std::endl;
		output<<"IncidentAngleRange(deg): "<<thetaMin/deg2rad<<" to "<<thetaMax/deg2rad<<" stepSize: "<<thetaStep/deg2rad<<std::endl;
		output<<std::setw(16)<<"Energy Deposited"<<"|"<<std::setw(16)<<"Intial Energy"<<std::endl;
		output<<"---------------------------------"<<std::endl;
		for(int i=0; i<ntbins; i++)
		{
			thetaCur = thetaMin + i*thetaStep;
			output<<"begin_theta "<<thetaCur/deg2rad<<std::endl;
			for(unsigned int j=0; j<energy_dep[i].size(); j++) 
			{
				output<<std::setw(16)<<energy_dep[i][j]<<" "<<std::setw(16)<<energy_in[i][j]<<std::endl;
			}
			output<<"end_theta"<<std::endl;
		}
		output.close();
	}

	if(options == "--test" || options == "")
	{
		std::cout<<"-------------Testing---------"<<std::endl;
		PunchTable table("tables/test.ptab");
		std::cout<<"Is the table valid? "<<table.IsValid()<<std::endl;
		double ke_test = 14.5; //MeV
		double theta_test = 35.0*deg2rad;
		double ke_dep = energyLoss.GetEnergyLoss(Zp, Ap, ke_test, thickness/std::fabs(std::cos(theta_test)));
		double recov_ke = table.GetInitialKineticEnergy(theta_test, ke_dep);
		std::cout<<"For a "<<MassLookup::GetInstance()->FindSymbol(Zp, Ap)<<" with kinetic energy "<<ke_test<<" MeV "
				 <<" the percent error on recovery is "<<(ke_test-recov_ke)/ke_test*100.0<<" where the recovered energy is "<<recov_ke<<" MeV"
				 <<" and the deposited energy is "<<ke_dep<<" MeV"<<std::endl;
		std::cout<<"-----------------------------"<<std::endl;
	}
}