#pragma once
#include <vector>
#include <Eigen/Core>
#include "FiniteElementCore.h"

struct StepResultModerno {
    Eigen::VectorXd u;
    double lambda;
};

class ArcLengthSolver {
private:
    int numSteps;
    int maxIter;
    double tol;
    double dL0;
    double Nd;

public:
    ArcLengthSolver(int steps = 30, int iterations = 50, double tolerance = 1e-6, double arcLength = 0.025, double targetIter = 5.0);

    std::vector<StepResultModerno> solve(const EstruturaModerna& est, const DOFManager& dofManager, const Eigen::VectorXd& fExterno);
};
