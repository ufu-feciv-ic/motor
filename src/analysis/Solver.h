#pragma once
#include <vector>
#include <Eigen/Core>
#include "FiniteElementCore.h"

/**
 * @brief Resultado de um passo da simulação.
 */
struct StepResultModerno {
    Eigen::VectorXd u;
    double lambda;
};

/**
 * @brief Classe base abstrata para todos os solvers.
 */
class SolverModerno {
public:
    virtual ~SolverModerno() = default;
    virtual std::vector<StepResultModerno> solve(
        const EstruturaModerna& est, 
        const DOFManager& dofManager, 
        const Eigen::VectorXd& fExterno
    ) = 0;
};

/**
 * @brief Solver Linear: resolve F = K * u em um único passo (Lambda = 1.0).
 */
class LinearSolver : public SolverModerno {
public:
    std::vector<StepResultModerno> solve(
        const EstruturaModerna& est, 
        const DOFManager& dofManager, 
        const Eigen::VectorXd& fExterno
    ) override;
};

/**
 * @brief Solver Arc-Length: resolve caminhos não lineares complexos.
 */
class ArcLengthSolver : public SolverModerno {
private:
    int numSteps;
    int maxIter;
    double tol;
    double dL0;
    double Nd;

public:
    ArcLengthSolver(int steps = 30, int iterations = 50, double tolerance = 1e-6, double arcLength = 0.025, double targetIter = 5.0);

    std::vector<StepResultModerno> solve(
        const EstruturaModerna& est, 
        const DOFManager& dofManager, 
        const Eigen::VectorXd& fExterno
    ) override;
};
