#pragma once
#include <vector>
#include <memory>
#include <map>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "DOFManager.h"

class NoModerno {
public:
    int id;
    double x, y;
    NoModerno(int id, double x, double y) : id(id), x(x), y(y) {}
};

struct EsforcosLocaisModernos {
    double N, V, M1, M2;
};

class ElementoFinitoModerno {
public:
    virtual ~ElementoFinitoModerno() = default;
    virtual std::vector<int> getIndicesGlobais(const DOFManager& dofManager) const = 0;
    virtual Eigen::MatrixXd getMatrizRigidez(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const = 0;
    virtual Eigen::VectorXd getForcasInternas(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const = 0;
    virtual EsforcosLocaisModernos getEsforcosLocais(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const = 0;
};

class EstruturaModerna {
public:
    std::vector<std::shared_ptr<NoModerno>> nos;
    std::vector<std::shared_ptr<ElementoFinitoModerno>> elementos;
    std::map<int, std::vector<int>> restricoes;

    void aplicarCondicoesContorno(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const DOFManager& dofManager) const;
};
