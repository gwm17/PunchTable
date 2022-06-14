#include "GenerateTable.h"
#include "MassLookup.h"
#include <fstream>
#include <iostream>
#include <cmath>

#include "catima/gwm_integrators.h"

namespace PunchTable {

    void GenerateTable(const TableParameters& params)
    {
        static double s_deg2rad = M_PI/180.0; //deg -> radians
        static double s_energyPrecision = 0.001; // keV precision
        static double s_ug2g = 1.0e-6; //ug -> g
        MassLookup& masses = MassLookup::GetInstance();

        catima::Projectile projectile(masses.FindMassU(params.projectileZ, params.projectileA), params.projectileZ, 0.0, 0.0);
        std::vector<catima::Material> materials;

        if(params.targetZ.size() != params.targetA.size() ||
           params.targetZ.size() != params.targetS.size() ||
           params.targetZ.size() != params.targetThickness.size())
        {
            std::cerr<<"ERR --  Invalid target parameters passed to GenerateTable. Mismatching target arrays."<<std::endl;
            return;
        }

        for(size_t i=0; i<params.targetZ.size(); i++)
        {
            auto& tZ = params.targetZ[i];
            auto& tA = params.targetA[i];
            auto& tS = params.targetS[i];
            materials.emplace_back();
            for(size_t j=0; j<tZ.size(); j++)
            {
                materials[i].add_element(masses.FindMassU(tZ[j], tA[j]), tZ[j], tS[j]);
            }
        }

        size_t depLayer = params.targetZ.size()-1;
        size_t thetaBins = (params.maxTheta - params.minTheta)/params.stepTheta;
        size_t energyBins = (params.maxKE - params.minKE)/params.stepKE;

        std::ofstream output(params.filename);
        if(!output.is_open())
        {
            std::cerr<<"ERR -- Unable to open output file "<<params.filename<<" at GenerateTable"<<std::endl;
            return;
        }
		output<<std::setprecision(5);

		output<<"Materials: "<<std::endl;;
		for(size_t i=0; i<depLayer; i++)
        {
            output<<"\tGoing through: ";
            for(size_t j=0; j<params.targetZ[i].size(); j++)
			    output<<masses.FindSymbol(params.targetZ[i][j], params.targetA[i][j])<<params.targetS[i][j];
            output<<" ("<<params.targetThickness[i]<<" ug/cm2)"<<std::endl;
        }
        output<<"\tDeposting into: ";
        for(size_t i=0; i<params.targetZ[depLayer].size(); i++)
            output<<masses.FindSymbol(params.targetZ[depLayer][i], params.targetA[depLayer][i])<<params.targetS[depLayer][i];
        output<<" ("<<params.targetThickness[depLayer]<<" ug/cm2)"<<std::endl;
		output<<"---------------------------------"<<std::endl;

		output<<"IncidentAngleRange(deg): "<<params.minTheta<<" to "<<params.maxTheta<<" stepSize: "<<params.stepTheta<<std::endl;
		output<<"---------------------------------"<<std::endl;

		output<<std::setw(16)<<"Energy Deposited"<<"|"<<std::setw(16)<<"Intial Energy"<<std::endl;
		output<<"---------------------------------"<<std::endl;

        std::vector<double> energyInData;
        std::vector<double> energyDepositedData;

        double theta, ke;
        double edep;
        double totalDep;
        bool stopped;

        size_t totalBins = thetaBins*energyBins;
        double flush_percent=0.01;
        size_t count=0, flush_val = flush_percent*totalBins, flush_count=0;
        for(size_t i=0; i < thetaBins; i++)
        {
            theta = (params.minTheta + i*params.stepTheta)*s_deg2rad;
            energyInData.clear();
            energyDepositedData.clear();
            for(size_t j=0; j<energyBins; j++)
            {
                count++;
                if(count == flush_val)
                {
                    count=0;
                    flush_count++;
                    std::cout<<"\rPercent of data generated: "<<flush_count*flush_percent*100.0<<"%"<<std::flush;
                }
                ke = params.maxKE - j*params.stepKE;
                projectile.T = ke/projectile.A;
                totalDep = 0.0;
                stopped = false;
                for(size_t k=0; k<depLayer; k++)
                {
                    materials[k].thickness(params.targetThickness[k]*s_ug2g/std::fabs(std::cos(theta)));
                    totalDep += catima::integrate_energyloss(projectile, materials[k]);
                    if((ke - totalDep) < s_energyPrecision)
                    {
                        stopped = true;
                        break;
                    }
                }

                if(stopped)
                    continue;

                materials[depLayer].thickness(params.targetThickness[depLayer]*s_ug2g/std::fabs(std::cos(theta)));
                edep = catima::integrate_energyloss(projectile, materials[depLayer]);
                if(std::abs(ke-edep) > s_energyPrecision)
                {
                    energyInData.push_back(ke);
                    energyDepositedData.push_back(edep);
                }
            }

            output<<"begin_theta "<<theta<<std::endl;
			for(size_t j=0; j<energyDepositedData.size(); j++) 
			{
				output<<std::setw(16)<<energyDepositedData[j]<<" "<<std::setw(16)<<energyInData[j]<<std::endl;
			}
			output<<"end_theta"<<std::endl;   
        }

        
		output.close();
    }

    double GetEnergyDeposited(const TableParameters& params, double thetaIncident, double initialKE)
    {
        static double s_energyPrecision = 0.001; // keV precision
        static double s_ug2g = 1.0e-6; //ug -> g
        MassLookup& masses = MassLookup::GetInstance();


        catima::Projectile projectile(masses.FindMassU(params.projectileZ, params.projectileA), params.projectileZ, 0.0, 0.0);
        std::vector<catima::Material> materials;

         if(params.targetZ.size() != params.targetA.size() ||
           params.targetZ.size() != params.targetS.size() ||
           params.targetZ.size() != params.targetThickness.size())
        {
            std::cerr<<"ERR --  Invalid target parameters passed to GenerateTable. Mismatching target arrays."<<std::endl;
            return 0.0;
        }

        for(size_t i=0; i<params.targetZ.size(); i++)
        {
            auto& tZ = params.targetZ[i];
            auto& tA = params.targetA[i];
            auto& tS = params.targetS[i];
            materials.emplace_back();
            for(size_t j=0; j<tZ.size(); j++)
            {
                materials[i].add_element(masses.FindMassU(tZ[j], tA[j]), tZ[j], tS[j]);
            }
        }

        size_t depLayer = params.targetZ.size() - 1;
        projectile.T = initialKE/projectile.A;
        double totalDep = 0.0;
        for(size_t k=0; k<depLayer; k++)
        {
            materials[k].thickness(params.targetThickness[k]*s_ug2g/std::fabs(std::cos(thetaIncident)));
            totalDep += catima::integrate_energyloss(projectile, materials[k]);
            if((initialKE - totalDep) < s_energyPrecision)
            {
                return 0.0;
            }
        }

        materials[depLayer].thickness(params.targetThickness[depLayer]*s_ug2g/std::fabs(std::cos(thetaIncident)));
        return catima::integrate_energyloss(projectile, materials[depLayer]);
    }
}