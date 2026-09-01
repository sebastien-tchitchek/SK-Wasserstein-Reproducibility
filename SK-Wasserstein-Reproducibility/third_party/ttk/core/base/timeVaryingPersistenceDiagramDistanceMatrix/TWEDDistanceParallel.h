#pragma once

#include <vector>
#include "PersistenceDiagramUtils.h"
#include "PersistenceDiagramDistanceMatrix.h"

void TWEDDistance2(double weight, double step, std::vector<std::vector<double>> const& costMatrix, std::vector<ttk::DiagramType> const& GeodesicSeti, std::vector<ttk::DiagramType> const& GeodesicSetj, std::vector<double>& GeodesicSetiTime, std::vector<double>& GeodesicSetjTime, double &distanceBetweenCurves, int &threadNumber, bool &Hilbert, int &choiceHilbertDistance, std::vector<double> &intraDistancePOTI, std::vector<double> &intraDistancePOTJ, std::vector<double> &intraDistanceSKOTI, std::vector<double> &intraDistanceSKOTJ, std::vector<double> &intraDistanceWSKI, std::vector<double> &intraDistanceWSKJ);
