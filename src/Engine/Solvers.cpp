#include "Engine/Solvers.h"
#include <Eigen/Eigenvalues>
#include <algorithm>

namespace Engine {

// --- AnaliseBuckling ---

std::vector<AnaliseBuckling::ModoBuckling> AnaliseBuckling::executar(Estrutura& est) {
    std::cout << "Iniciando Analise de Flambagem Linear...\n";

    Eigen::VectorXd uZero = Eigen::VectorXd::Zero(est.NumGDLs);
    Eigen::MatrixXd K = Construtor::montarMatrizRigidezGlobal(est, uZero);
    Eigen::VectorXd F = Construtor::montarVetorForcasReferencia(est);
    
    Eigen::MatrixXd K_solver = K;
    Eigen::VectorXd F_solver = F;
    est.aplicarCondicoesContorno(K_solver, F_solver);
    Eigen::VectorXd uRef = K_solver.colPivHouseholderQr().solve(F_solver);

    Eigen::VectorXd Nref = Eigen::VectorXd::Zero(est.Elementos.size());
    for (size_t i = 0; i < est.Elementos.size(); ++i) {
        EsforcosLocais esf = est.Elementos[i]->getEsforcosLocais(uRef);
        Nref(i) = esf.N;
    }

    Eigen::MatrixXd KG = Construtor::montarMatrizRigidezGeometricaGlobal(est, Nref);

    std::vector<int> livres;
    std::vector<bool> eFixo(est.NumGDLs, false);
    for (int dof : est.NosFixos) eFixo[dof] = true;
    for (int i = 0; i < est.NumGDLs; ++i) if (!eFixo[i]) livres.push_back(i);

    int nf = livres.size();
    Eigen::MatrixXd Kff(nf, nf), KGff(nf, nf);
    for (int i = 0; i < nf; ++i) {
        for (int j = 0; j < nf; ++j) {
            Kff(i, j) = K(livres[i], livres[j]);
            KGff(i, j) = KG(livres[i], livres[j]);
        }
    }

    Eigen::MatrixXd C = Kff.ldlt().solve(-KGff);
    Eigen::EigenSolver<Eigen::MatrixXd> es(C);

    std::vector<ModoBuckling> modos;
    for (int i = 0; i < es.eigenvalues().size(); ++i) {
        double mu = es.eigenvalues()(i).real();
        if (std::abs(es.eigenvalues()(i).imag()) < 1e-10 && mu > 1e-10) {
            double lambda = 1.0 / mu;
            Eigen::VectorXd modoF = es.eigenvectors().col(i).real();
            Eigen::VectorXd modoG = Eigen::VectorXd::Zero(est.NumGDLs);
            for (int j = 0; j < nf; ++j) modoG(livres[j]) = modoF(j);
            modos.push_back({modoG, lambda});
        }
    }

    std::sort(modos.begin(), modos.end(), [](const ModoBuckling& a, const ModoBuckling& b) {
        return a.lambda < b.lambda;
    });

    return modos;
}

// --- AnaliseLinear ---

std::vector<Resultado> AnaliseLinear::executar(Estrutura& est) {
    std::cout << "Iniciando solucao do sistema linear...\n";

    Eigen::VectorXd uZero = Eigen::VectorXd::Zero(est.NumGDLs);
    Eigen::MatrixXd KGlobalOriginal = Construtor::montarMatrizRigidezGlobal(est, uZero);
    Eigen::MatrixXd KGlobalEst = KGlobalOriginal;
    Eigen::VectorXd FRef = Construtor::montarVetorForcasReferencia(est);
    Eigen::VectorXd FGlobal = FRef;

    est.aplicarCondicoesContorno(KGlobalEst, FGlobal);
    Eigen::VectorXd uGlobal = KGlobalEst.colPivHouseholderQr().solve(FGlobal);
    Eigen::VectorXd R = KGlobalOriginal * uGlobal - FRef;

    Resultado resultado;
    resultado.u = uGlobal;
    resultado.reacoes = R;
    resultado.FatorCarga = 1.0;

    return {resultado};
}

// --- AnaliseNaoLinearNR ---

AnaliseNaoLinearNR::AnaliseNaoLinearNR(int passos, int iteracoes, double tolerancia)
    : numPassos(passos), maxIter(iteracoes), tol(tolerancia) {}

std::vector<Resultado> AnaliseNaoLinearNR::executar(Estrutura& est) {
    std::cout << "--- Iniciando Solver Newton-Raphson ---\n";
    std::vector<Resultado> historico;

    Eigen::VectorXd uAtual = Eigen::VectorXd::Zero(est.NumGDLs);
    historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

    Eigen::VectorXd FRef = Construtor::montarVetorForcasReferencia(est);

    for (int passo = 1; passo <= numPassos; ++passo) {
        double lambda = (double)passo / numPassos;
        Eigen::VectorXd PassoCarga = lambda * FRef;

        int iter = 0;
        double erro = 1.0;

        while (iter < maxIter) {
            Eigen::MatrixXd KTangente = Construtor::montarMatrizRigidezGlobal(est, uAtual);
            Eigen::VectorXd Fint = Construtor::montarForcasInternasGlobais(est, uAtual);
            Eigen::VectorXd g = PassoCarga - Fint;

            for (int dof : est.NosFixos) g(dof) = 0.0;
            erro = g.norm();
            if (erro < tol) break;

            est.aplicarCondicoesContorno(KTangente, g);
            Eigen::VectorXd deltaU = KTangente.ldlt().solve(g);
            uAtual += deltaU;
            iter++;
        }

        Eigen::VectorXd FintFinal = Construtor::montarForcasInternasGlobais(est, uAtual);
        Eigen::VectorXd R = FintFinal - PassoCarga;

        historico.push_back(Resultado{uAtual, PassoCarga, R, lambda});
    }
    return historico;
}

// --- AnaliseNaoLinearCompArco ---

AnaliseNaoLinearCompArco::AnaliseNaoLinearCompArco(
    int passos, int iteracoes, double tolerancia, double comprimentoArcoInicial, double iterDesejadas)
    : numPassos(passos), maxIter(iteracoes), tol(tolerancia), deltal0(comprimentoArcoInicial), Nd(iterDesejadas) {}

std::vector<Resultado> AnaliseNaoLinearCompArco::executar(Estrutura& est) {
    std::cout << "--- Iniciando Solver Arc-Length ---\n";
    std::vector<Resultado> historico;
    Eigen::VectorXd uAtual = Eigen::VectorXd::Zero(est.NumGDLs);
    Eigen::VectorXd DELTAU = Eigen::VectorXd::Zero(est.NumGDLs);
    double lambda = 0.0;
    double deltal = deltal0;

    historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

    Eigen::VectorXd FRef = Construtor::montarVetorForcasReferencia(est);

    for (int passo = 1; passo <= numPassos; ++passo) {
        Eigen::MatrixXd Kt = Construtor::montarMatrizRigidezGlobal(est, uAtual);
        Eigen::VectorXd FRef_solver = FRef;
        est.aplicarCondicoesContorno(Kt, FRef_solver);
        Eigen::VectorXd deltaUr = Kt.colPivHouseholderQr().solve(FRef_solver);

        double normDur = deltaUr.norm();
        if (normDur < 1e-12) normDur = 1e-12;
        double DeltaLambda = deltal / normDur;
        if (passo > 1 && DELTAU.dot(deltaUr) < 0.0) DeltaLambda = -DeltaLambda;

        Eigen::VectorXd deltaU = DeltaLambda * deltaUr;
        Eigen::VectorXd DELTAU0 = deltaU;
        DELTAU = deltaU;
        double lambdaTentativo = lambda + DeltaLambda;

        int iter = 0;
        while (iter < maxIter) {
            Eigen::VectorXd Fint = Construtor::montarForcasInternasGlobais(est, uAtual + DELTAU);
            Eigen::VectorXd g = lambdaTentativo * FRef - Fint;
            for (int dof : est.NosFixos) g(dof) = 0.0;
            if (g.norm() <= tol) break;

            Kt = Construtor::montarMatrizRigidezGlobal(est, uAtual + DELTAU);
            Eigen::VectorXd FRef_iter = FRef;
            est.aplicarCondicoesContorno(Kt, FRef_iter);
            Eigen::VectorXd deltaUg = Kt.colPivHouseholderQr().solve(g);
            deltaUr = Kt.colPivHouseholderQr().solve(FRef_iter);

            double topo = DELTAU0.dot(deltaUg);
            double base = DELTAU0.dot(deltaUr);
            double dlambda = (std::abs(base) < 1e-12) ? 0.0 : -topo / base;

            DELTAU += deltaUg + dlambda * deltaUr;
            DeltaLambda += dlambda;
            lambdaTentativo = lambda + DeltaLambda;
            iter++;
        }

        uAtual += DELTAU;
        lambda = lambdaTentativo;
        if (iter > 0) deltal = deltal0 * std::sqrt(Nd / static_cast<double>(iter));

        Eigen::VectorXd FintFinal = Construtor::montarForcasInternasGlobais(est, uAtual);
        Eigen::VectorXd R = FintFinal - (lambda * FRef);
        historico.push_back({uAtual, lambda * FRef, R, lambda});
    }
    return historico;
}

} // namespace Engine
