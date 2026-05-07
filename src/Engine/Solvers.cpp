#include "Engine/Solvers.h"
#include <Eigen/Eigenvalues>
#include <Eigen/SparseLU>
#include <algorithm>

namespace Engine {

// --- AnaliseBuckling ---

std::vector<AnaliseBuckling::ModoBuckling> AnaliseBuckling::executar(Estrutura& est) {
    std::cout << "Iniciando Analise de Flambagem Linear (Esparsa)...\n";

    Eigen::VectorXd uZero = Eigen::VectorXd::Zero(est.NumGDLs);
    Eigen::SparseMatrix<double> K = Construtor::montarMatrizRigidezGlobal(est, uZero);
    Eigen::VectorXd F = Construtor::montarVetorForcasReferencia(est);
    
    Eigen::SparseMatrix<double> K_solver = K;
    Eigen::VectorXd F_solver = F;
    est.aplicarCondicoesContorno(K_solver, F_solver);
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(K_solver);
    Eigen::VectorXd uRef = lu.solve(F_solver);

    Eigen::VectorXd Nref = Eigen::VectorXd::Zero(est.Elementos.size());
    for (size_t i = 0; i < est.Elementos.size(); ++i) {
        EsforcosLocais esf = est.Elementos[i]->getEsforcosLocais(uRef);
        Nref(i) = esf.N;
    }

    Eigen::SparseMatrix<double> KG = Construtor::montarMatrizRigidezGeometricaGlobal(est, Nref);

    std::vector<int> livres;
    std::vector<bool> eFixo(est.NumGDLs, false);
    for (int dof : est.NosFixos) eFixo[dof] = true;
    for (int i = 0; i < est.NumGDLs; ++i) if (!eFixo[i]) livres.push_back(i);

    int nf = livres.size();
    // Para autovalores, como o problema costuma ser pequeno nos testes, 
    // convertemos as submatrizes livres para denso.
    Eigen::MatrixXd Kff = Eigen::MatrixXd::Zero(nf, nf);
    Eigen::MatrixXd KGff = Eigen::MatrixXd::Zero(nf, nf);
    
    for (int i = 0; i < nf; ++i) {
        for (int j = 0; j < nf; ++j) {
            Kff(i, j) = K.coeff(livres[i], livres[j]);
            KGff(i, j) = KG.coeff(livres[i], livres[j]);
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
    std::cout << "Iniciando solucao do sistema linear (Esparsa)...\n";

    Eigen::VectorXd uZero = Eigen::VectorXd::Zero(est.NumGDLs);
    Eigen::SparseMatrix<double> KGlobal = Construtor::montarMatrizRigidezGlobal(est, uZero);
    Eigen::VectorXd FRef = Construtor::montarVetorForcasReferencia(est);

    Eigen::SparseMatrix<double> K_solver = KGlobal;
    Eigen::VectorXd F_solver = FRef;
    est.aplicarCondicoesContorno(K_solver, F_solver);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(K_solver);
    Eigen::VectorXd uGlobal = lu.solve(F_solver);
    
    // Reacoes: R = Fint - Fext. Para linear, Fint = K * u.
    Eigen::VectorXd R = KGlobal * uGlobal - FRef;

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
    std::cout << "--- Iniciando Solver Newton-Raphson (Esparsa) ---\n";
    std::vector<Resultado> historico;

    Eigen::VectorXd uAtual = Eigen::VectorXd::Zero(est.NumGDLs);
    historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

    Eigen::VectorXd FRef = Construtor::montarVetorForcasReferencia(est);

    for (int passo = 1; passo <= numPassos; ++passo) {
        double lambda = (double)passo / numPassos;
        Eigen::VectorXd PassoCarga = lambda * FRef;

        int iter = 0;
        while (iter < maxIter) {
            Eigen::SparseMatrix<double> KTangente = Construtor::montarMatrizRigidezGlobal(est, uAtual);
            Eigen::VectorXd Fint = Construtor::montarForcasInternasGlobais(est, uAtual);
            Eigen::VectorXd g = PassoCarga - Fint;

            for (int dof : est.NosFixos) g(dof) = 0.0;
            if (g.norm() < tol) break;

            est.aplicarCondicoesContorno(KTangente, g);
            
            Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
            lu.compute(KTangente);
            Eigen::VectorXd deltaU = lu.solve(g);
            
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
    std::cout << "--- Iniciando Solver Arc-Length (Esparsa) ---\n";
    std::vector<Resultado> historico;
    Eigen::VectorXd uAtual = Eigen::VectorXd::Zero(est.NumGDLs);
    Eigen::VectorXd DELTAU = Eigen::VectorXd::Zero(est.NumGDLs);
    double lambda = 0.0;
    double deltal = deltal0;

    historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

    Eigen::VectorXd FRef = Construtor::montarVetorForcasReferencia(est);

    for (int passo = 1; passo <= numPassos; ++passo) {
        Eigen::SparseMatrix<double> Kt = Construtor::montarMatrizRigidezGlobal(est, uAtual);
        Eigen::VectorXd FRef_solver = FRef;
        est.aplicarCondicoesContorno(Kt, FRef_solver);
        
        Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
        lu.compute(Kt);
        Eigen::VectorXd deltaUr = lu.solve(FRef_solver);

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
            
            lu.compute(Kt);
            Eigen::VectorXd deltaUg = lu.solve(g);
            deltaUr = lu.solve(FRef_iter);

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
