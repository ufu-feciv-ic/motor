#include "Solver.h"
#include "SparseAssembler.h"
#include <Eigen/SparseLU>
#include <iostream>

// --- LINEAR SOLVER ---

std::vector<StepResultModerno> LinearSolver::solve(const EstruturaModerna& est, const DOFManager& dofManager, const Eigen::VectorXd& fExterno) {
    int totalDofs = dofManager.getNumTotalDofs();
    if (totalDofs <= 0) return {{Eigen::VectorXd::Zero(0), 0.0}};

    Eigen::VectorXd uZero = Eigen::VectorXd::Zero(totalDofs);
    
    std::vector<SparseAssembler::Triplet> triplets;
    for (const auto& el : est.elementos) {
        auto indices = el->getIndicesGlobais(dofManager);
        auto Kloc = el->getMatrizRigidez(uZero, dofManager);
        SparseAssembler::adicionarContribuicao(triplets, indices, Kloc);
    }
    
    if (triplets.empty()) {
        return {{uZero, 0.0}};
    }

    auto K = SparseAssembler::construirMatriz(totalDofs, triplets);
    Eigen::VectorXd F = fExterno;
    est.aplicarCondicoesContorno(K, F, dofManager);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);
    if (solver.info() != Eigen::Success) {
        std::cerr << "ERRO: Matriz de rigidez linear singular ou mal condicionada!" << std::endl;
        return {{uZero, 0.0}};
    }

    Eigen::VectorXd uFinal = solver.solve(F);
    if (solver.info() != Eigen::Success) {
        std::cerr << "ERRO: Falha ao resolver sistema linear!" << std::endl;
        return {{uZero, 0.0}};
    }

    return {{uZero, 0.0}, {uFinal, 1.0}};
}

// --- ARC LENGTH SOLVER ---

ArcLengthSolver::ArcLengthSolver(int steps, int iterations, double tolerance, double arcLength, double targetIter)
    : numSteps(steps), maxIter(iterations), tol(tolerance), dL0(arcLength), Nd(targetIter) {}

std::vector<StepResultModerno> ArcLengthSolver::solve(const EstruturaModerna& est, const DOFManager& dofManager, const Eigen::VectorXd& fExterno) {
    int totalDofs = dofManager.getNumTotalDofs();
    if (totalDofs <= 0) return {{Eigen::VectorXd::Zero(0), 0.0}};

    Eigen::VectorXd u = Eigen::VectorXd::Zero(totalDofs);
    Eigen::VectorXd DU = Eigen::VectorXd::Zero(totalDofs);
    double lambda = 0.0;
    double dL = dL0;

    std::vector<StepResultModerno> history;
    history.push_back({u, 0.0});

    for (int step = 1; step <= numSteps; ++step) {
        std::vector<SparseAssembler::Triplet> triplets;
        for (const auto& el : est.elementos) {
            SparseAssembler::adicionarContribuicao(triplets, el->getIndicesGlobais(dofManager), el->getMatrizRigidez(u, dofManager));
        }
        
        if (triplets.empty()) break;
        auto Kt = SparseAssembler::construirMatriz(totalDofs, triplets);
        
        Eigen::VectorXd fTemp = fExterno;
        est.aplicarCondicoesContorno(Kt, fTemp, dofManager);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(Kt);
        if (solver.info() != Eigen::Success) {
            std::cerr << "ERRO: Fatoracao falhou no passo " << step << std::endl;
            break;
        }

        Eigen::VectorXd dUr = solver.solve(fExterno);
        double nDur = dUr.norm();
        if (nDur < 1e-12) nDur = 1e-12;
        double dLambda = dL / nDur;

        if (step > 1 && DU.dot(dUr) < 0) dLambda = -dLambda;

        Eigen::VectorXd dU = dLambda * dUr;
        Eigen::VectorXd DU0 = dU;
        DU = dU;
        double lambdaIter = lambda + dLambda;

        int iter = 0;
        bool convergiu = false;
        while (iter < maxIter) {
            Eigen::VectorXd fInt = Eigen::VectorXd::Zero(totalDofs);
            for (const auto& el : est.elementos) {
                auto indices = el->getIndicesGlobais(dofManager);
                Eigen::VectorXd fi = el->getForcasInternas(u + DU, dofManager);
                for(int i=0; i<(int)indices.size(); ++i) {
                    if (indices[i] >= 0 && indices[i] < totalDofs) fInt(indices[i]) += fi(i);
                }
            }

            Eigen::VectorXd g = lambdaIter * fExterno - fInt;
            for (auto const& [nodeId, restDofs] : est.restricoes) {
                auto indices = dofManager.obterIndicesGlobais(nodeId);
                for (int d : restDofs) {
                    if (d >= 0 && d < (int)indices.size()) {
                        int idx = indices[d];
                        if (idx >= 0 && idx < totalDofs) g(idx) = 0.0;
                    }
                }
            }

            if (g.norm() < tol) { convergiu = true; break; }

            triplets.clear();
            for (const auto& el : est.elementos) {
                SparseAssembler::adicionarContribuicao(triplets, el->getIndicesGlobais(dofManager), el->getMatrizRigidez(u + DU, dofManager));
            }
            Kt = SparseAssembler::construirMatriz(totalDofs, triplets);
            Eigen::VectorXd fDummy = fExterno;
            est.aplicarCondicoesContorno(Kt, fDummy, dofManager);
            
            solver.compute(Kt);
            if (solver.info() != Eigen::Success) break;

            Eigen::VectorXd dUg = solver.solve(g);
            dUr = solver.solve(fExterno);

            double denom = DU0.dot(dUr);
            if (std::abs(denom) < 1e-18) break;

            double dlambda = -DU0.dot(dUg) / denom;
            DU += dUg + dlambda * dUr;
            dLambda += dlambda;
            lambdaIter = lambda + dLambda;
            iter++;
        }

        if (!convergiu) {
            dL *= 0.5;
            if (dL < 1e-8) break;
            continue;
        }

        u += DU;
        lambda = lambdaIter;
        if (iter > 0) dL = dL * std::sqrt(Nd / (double)iter);
        if (dL > dL0 * 2.0) dL = dL0 * 2.0;

        history.push_back({u, lambda});
    }

    return history;
}
