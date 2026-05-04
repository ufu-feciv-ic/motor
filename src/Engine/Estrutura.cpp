#include "Engine/Estrutura.h"
#include <algorithm>
#include <cmath>

namespace Engine {

void Estrutura::adicionarNo(std::shared_ptr<No> no) {
    Nos.push_back(no);
    NumGDLs += no->gdlGlobais.size();
}

void Estrutura::adicionarElemento(std::shared_ptr<ElementoFinito> elemento) {
    Elementos.push_back(elemento);
}

void Estrutura::aplicarCondicoesContorno(Eigen::MatrixXd& K, Eigen::VectorXd& F) const {
    for (int GDL : NosFixos) {
        K.row(GDL).setZero();
        K.col(GDL).setZero();
        K(GDL, GDL) = 1.0;
        F(GDL) = 0.0;
    }
}

void Estrutura::perturbarGeometria(const Eigen::VectorXd& modo, double amplitude) {
    double maxDesl = 0.0;
    for (const auto& no : Nos) {
        double dx = modo(no->gdlGlobais[0]);
        double dy = modo(no->gdlGlobais[1]);
        maxDesl = std::max(maxDesl, std::sqrt(dx * dx + dy * dy));
    }
    if (maxDesl < 1e-12) return;

    for (auto& no : Nos) {
        no->x += amplitude * (modo(no->gdlGlobais[0]) / maxDesl);
        no->y += amplitude * (modo(no->gdlGlobais[1]) / maxDesl);
    }

    for (auto& elem : Elementos) {
        elem->atualizarGeometria();
    }
}

Eigen::MatrixXd Construtor::montarMatrizRigidezGlobal(const Estrutura& est, const Eigen::VectorXd& uGlobal) {
    Eigen::MatrixXd KGlobalEst = Eigen::MatrixXd::Zero(est.NumGDLs, est.NumGDLs);

    for (const auto& elemento : est.Elementos) {
        Eigen::MatrixXd KGlobalElem = elemento->getMatrizRigidezGlobal(uGlobal);
        std::vector<int> GDLs = elemento->getGDLsGlobais();

        for (size_t i = 0; i < GDLs.size(); ++i)
            for (size_t j = 0; j < GDLs.size(); ++j)
                KGlobalEst(GDLs[i], GDLs[j]) += KGlobalElem(i, j);
    }

    return KGlobalEst;
}

Eigen::VectorXd Construtor::montarForcasInternasGlobais(const Estrutura& est, const Eigen::VectorXd& uGlobal) {
    Eigen::VectorXd FGlobal = Eigen::VectorXd::Zero(est.NumGDLs);

    for (const auto& elemento : est.Elementos) {
        Eigen::VectorXd FGlobalElem = elemento->getForcasInternasGlobais(uGlobal);
        std::vector<int> GDLs = elemento->getGDLsGlobais();

        for (size_t i = 0; i < GDLs.size(); ++i)
            FGlobal(GDLs[i]) += FGlobalElem(i);
    }

    return FGlobal;
}

Eigen::MatrixXd Construtor::montarMatrizRigidezGeometricaGlobal(const Estrutura& est, const Eigen::VectorXd& esforcosNormais) {
    Eigen::MatrixXd KGGlobalEst = Eigen::MatrixXd::Zero(est.NumGDLs, est.NumGDLs);

    for (size_t i = 0; i < est.Elementos.size(); ++i) {
        const auto& elemento = est.Elementos[i];
        Eigen::MatrixXd KGElem = elemento->getMatrizRigidezGeometrica(esforcosNormais(i));
        std::vector<int> GDLs = elemento->getGDLsGlobais();

        for (size_t row = 0; row < GDLs.size(); ++row)
            for (size_t col = 0; col < GDLs.size(); ++col)
                KGGlobalEst(GDLs[row], GDLs[col]) += KGElem(row, col);
    }

    return KGGlobalEst;
}

Eigen::VectorXd Construtor::montarVetorForcasReferencia(const Estrutura& est) {
    Eigen::VectorXd FRef = est.ForcasExternas;

    for (const auto& elemento : est.Elementos) {
        Eigen::VectorXd fEquiv = elemento->getForcasEquivalentes();
        std::vector<int> GDLs = elemento->getGDLsGlobais();

        for (size_t i = 0; i < GDLs.size(); ++i) {
            FRef(GDLs[i]) += fEquiv(i);
        }
    }

    return FRef;
}

} // namespace Engine
