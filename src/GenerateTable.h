#ifndef GENERATE_TABLE_H
#define GENERATE_TABLE_H

#include <string>
#include <vector>

namespace PunchTable {

    struct TableParameters
    {
        std::string tableType; //Either PunchTable or ElossTable
        double minKE; //MeV
        double maxKE; //MeV
        double stepKE; //MeV
        double minTheta; //deg
        double maxTheta; //deg
        double stepTheta; //deg
        int projectileZ;
        int projectileA;
        std::vector<std::vector<int>> targetZ; //Last element is the deposition layer for PunchTable
        std::vector<std::vector<int>> targetA;
        std::vector<std::vector<int>> targetS; //Stoichiometry
        std::vector<double> targetThickness; //ug/cm2
        std::string filename;
    };

    TableParameters GetTableParameters(const std::string& filename);

    void GeneratePunchTable(const TableParameters& params);

    void GenerateElossTable(const TableParameters& params);

    //For testing (thetaIncident = radians, initialKE = MeV)
    double GetEnergyDeposited(const TableParameters& params, double thetaIncident, double initialKE);
}

#endif