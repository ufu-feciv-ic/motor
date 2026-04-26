#include "Solver.h"
#include "SparseAssembler.h"
#include <Eigen/SparseLU>
#include <iostream>

ArcLengthSolver::ArcLengthSolver(int steps, int iterations, double tolerance, double arcLength, double targetIter)
    : numSteps(steps), maxIter(iterations), tol(tolerance), dL0(arcLength), Nd(targetIter) {}

std::vector<StepResultModerno> ArcLengthSolver::solve(const EstruturaModerna& est, const DOFManager& dofManager, const Eigen::VectorXd& fExterno) {
    int totalDofs = dofManager.getNumTotalDofs();
    Eigen::VectorXd u = Eigen::VectorXd::Zero(totalDofs);
    Eigen::VectorXd DU = Eigen::VectorXd::Zero(totalDofs);
    double lambda = 0.0;
    double dL = dL0;

    std::vector<StepResultModerno> history;
    history.push_back({u, 0.0});

    for (int step = 1; step <= numSteps; ++step) {
        // --- PREDITOR ---
        std::vector<SparseAssembler::Triplet> triplets;
        for (const auto& el : est.elementos) {
            SparseAssembler::adicionarContribuicao(triplets, el->getIndicesGlobais(dofManager), el->getMatrizRigidez(u, dofManager));
        }
        auto Kt = SparseAssembler::construirMatriz(totalDofs, triplets);
        Eigen::VectorXd fTemp = fExterno;
        est.aplicarCondicoesContorno(Kt, fTemp, dofManager);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(Kt);
        Eigen::VectorXd dUr = solver.solve(fExterno);

        double nDur = dUr.norm();
        if (nDur < 1e-12) nDur = 1e-12;
        double dLambda = dL / nDur;

        if (step > 1 && DU.dot(dUr) < 0) dLambda = -dLambda;

        Eigen::VectorXd dU = dLambda * dUr;
        Eigen::VectorXd DU0 = dU;
        DU = dU;
        double lambdaIter = lambda + dLambda;

        // --- CORRETOR ---
        int iter = 0;
        while (iter < maxIter) {
            Eigen::VectorXd fInt = Eigen::VectorXd::Zero(totalDofs);
            for (const auto& el : est.elementos) {
                auto indices = el->getIndicesGlobais(dofManager);
                Eigen::VectorXd fi = el->getForcasInternas(u + DU, dofManager);
                for(int i=0; i<(int)indices.size(); ++i) fInt(indices[i]) += fi(i);
            }

            Eigen::VectorXd g = lambdaIter * fExterno - fInt;
            // Zerar resíduo nas restrições (manualmente para o vetor)
            for (auto const& [nodeId, dofs] : est.restricoes) {
                auto indices = dofManager.obterIndicesGlobais(nodeId);
                for (int d : dofs) g(indices[d]) = 0.0;
            }

            if (g.norm() < tol) break;

            triplets.clear();
            for (const auto& el : est.elementos) {
                SparseAssembler::adicionarContribuicao(triplets, el->getIndicesGlobais(dofManager), el->getMatrizRigidez(u + DU, dofManager));
            }
            Kt = SparseAssembler::construirMatriz(totalDofs, triplets);
            fTemp = fExterno;
            est.aplicarCondicoesContorno(Kt, fTemp, dofManager);
            
            solver.compute(Kt);
            Eigen::VectorXd dUg = solver.solve(g);
            dUr = solver.solve(fExterno);

            double dlambda = -DU0.dot(dUg) / DU0.dot(dUr);
            DU += dUg + dlambda * dUr;
            dLambda += dlambda;
            lambdaIter = lambda + dLambda;
            iter++;
        }

        u += DU;
        lambda = lambdaIter;
        if (iter > 0) dL = dL0 * std::sqrt(Nd / iter);
        history.push_back({u, lambda});
        std::cout << "Passo " << step << " concluido. Lambda: " << lambda << std::endl;
    }

    return history;
}
