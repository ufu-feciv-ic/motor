#ifndef ENGINE_ELEMENTOFINITO_H
#define ENGINE_ELEMENTOFINITO_H

#include <memory>
#include <vector>
#include "Engine/Common.h"
#include "Engine/No.h"

namespace Engine {

class ElementoFinito {
public:
    virtual ~ElementoFinito() = default;
    virtual std::vector<int> getGDLsGlobais() const = 0;
    virtual Eigen::MatrixXd getMatrizRigidezGlobal(const Eigen::VectorXd& uGlobal) const = 0;
    virtual Eigen::VectorXd getForcasInternasGlobais(const Eigen::VectorXd& uGlobal) const = 0;
    virtual EsforcosLocais getEsforcosLocais(const Eigen::VectorXd& uGlobal) const = 0;
    virtual Eigen::MatrixXd getMatrizRigidezGeometrica(double N) const = 0;
    virtual void atualizarGeometria() = 0;
};

class Viga2DLinear : public ElementoFinito {
private:
    std::shared_ptr<No> n1, n2;
    PropriedadesMaterial material;
    double Linicial, cosInicial, senInicial;

public:
    Viga2DLinear(std::shared_ptr<No> n1, std::shared_ptr<No> n2, PropriedadesMaterial material);
    void atualizarGeometria() override;
    std::vector<int> getGDLsGlobais() const override;
    Eigen::MatrixXd getMatrizRigidezGlobal(const Eigen::VectorXd& uGlobal) const override;
    Eigen::VectorXd getForcasInternasGlobais(const Eigen::VectorXd& uGlobal) const override;
    EsforcosLocais getEsforcosLocais(const Eigen::VectorXd& uGlobal) const override;
    Eigen::MatrixXd getMatrizRigidezGeometrica(double N) const override;

private:
    Eigen::MatrixXd calcularkLocal() const;
    Eigen::MatrixXd calcularTransformacao() const;
};

class Viga2DCorrotacional : public ElementoFinito {
private:
    std::shared_ptr<No> n1, n2;
    PropriedadesMaterial mat;
    double L0;

public:
    Viga2DCorrotacional(std::shared_ptr<No> no1, std::shared_ptr<No> no2, PropriedadesMaterial m);
    void atualizarGeometria() override;
    std::vector<int> getGDLsGlobais() const override;
    Eigen::MatrixXd getMatrizRigidezGlobal(const Eigen::VectorXd& uGlobal) const override;
    Eigen::VectorXd getForcasInternasGlobais(const Eigen::VectorXd& uGlobal) const override;
    EsforcosLocais getEsforcosLocais(const Eigen::VectorXd& uGlobal) const override;
    Eigen::MatrixXd getMatrizRigidezGeometrica(double N) const override;
};

} // namespace Engine

#endif // ENGINE_ELEMENTOFINITO_H
