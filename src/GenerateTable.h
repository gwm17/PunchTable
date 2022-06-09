#ifndef GENERATE_TABLE_H
#define GENERATE_TABLE_H

#include <string>
#include <vector>

namespace PunchTable {

    struct TableParameters
    {
        double minKE;
        double maxKE;
        double stepKE;
        double minTheta;
        double maxTheta;
        double stepTheta;
        int projectileZ;
        int projectileA;
        std::vector<std::vector<int>> targetZ;
        std::vector<std::vector<int>> targetA;
        std::vector<std::vector<int>> targetS;
        double targetThickness;
        std::string filename;
    };

    void GenerateTable(const TableParameters& params);
}

#endif