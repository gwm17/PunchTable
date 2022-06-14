#ifndef GENERATE_TABLE_H
#define GENERATE_TABLE_H

#include <string>
#include <vector>

namespace PunchTable {

    struct TableParameters
    {
        double minKE; //MeV
        double maxKE; //MeV
        double stepKE; //MeV
        double minTheta; //deg
        double maxTheta; //deg
        double stepTheta; //deg
        int projectileZ;
        int projectileA;
        std::vector<std::vector<int>> targetZ; //Last element is the deposition layer
        std::vector<std::vector<int>> targetA;
        std::vector<std::vector<int>> targetS; //Stoichiometry
        std::vector<double> targetThickness; //ug/cm2
        std::string filename;
    };

    void GenerateTable(const TableParameters& params);

    //For testing (thetaIncident = radians, initialKE = MeV)
    double GetEnergyDeposited(const TableParameters& params, double thetaIncident, double initialKE);
}

#endif