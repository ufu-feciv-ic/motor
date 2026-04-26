#pragma once
#include "FiniteElementCore.h"
#include "StructuralModel.h"

class Viga2DCorrotacionalModerna : public ElementoFinitoModerno {
private:
    std::shared_ptr<NoModerno> n1, n2;
    double E, A, I;
    double L0;

public:
    Viga2DCorrotacionalModerna(std::shared_ptr<NoModerno> node1, std::shared_ptr<NoModerno> node2, const MaterialData& mat);

    std::vector<int> getIndicesGlobais(const DOFManager& dofManager) const override;
    std::pair<int, int> getNodeIds() const override { return {n1->id, n2->id}; }
    Eigen::MatrixXd getMatrizRigidez(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const override;
    Eigen::VectorXd getForcasInternas(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const override;
    EsforcosLocaisModernos getEsforcosLocais(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const override;
    
    static double normalizarAngulo(double angulo);
};

class Viga2DLinearModerna : public ElementoFinitoModerno {
private:
    std::shared_ptr<NoModerno> n1, n2;
    double E, A, I;
    double L0, cos0, sin0;

public:
    Viga2DLinearModerna(std::shared_ptr<NoModerno> node1, std::shared_ptr<NoModerno> node2, const MaterialData& mat);

    std::vector<int> getIndicesGlobais(const DOFManager& dofManager) const override;
    Eigen::MatrixXd getMatrizRigidez(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const override;
    Eigen::VectorXd getForcasInternas(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const override;
    EsforcosLocaisModernos getEsforcosLocais(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const override;
};
