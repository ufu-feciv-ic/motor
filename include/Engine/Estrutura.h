#ifndef ENGINE_ESTRUTURA_H
#define ENGINE_ESTRUTURA_H

#include <vector>
#include <memory>
#include <Eigen/Sparse>
#include "Engine/No.h"
#include "Engine/ElementoFinito.h"

namespace Engine {

class Estrutura {
public:
    int NumGDLs = 0;
    std::vector<std::shared_ptr<No>> Nos;
    std::vector<std::shared_ptr<ElementoFinito>> Elementos;
    Eigen::VectorXd ForcasExternas;
    std::vector<int> NosFixos;

    void adicionarNo(std::shared_ptr<No> no);
    void adicionarElemento(std::shared_ptr<ElementoFinito> elemento);
    void aplicarCondicoesContorno(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const;
    void perturbarGeometria(const Eigen::VectorXd& modo, double amplitude);
};

class Construtor {
public:
    static Eigen::SparseMatrix<double> montarMatrizRigidezGlobal(const Estrutura& est, const Eigen::VectorXd& uGlobal);
    static Eigen::VectorXd montarForcasInternasGlobais(const Estrutura& est, const Eigen::VectorXd& uGlobal);
    static Eigen::SparseMatrix<double> montarMatrizRigidezGeometricaGlobal(const Estrutura& est, const Eigen::VectorXd& esforcosNormais);
    static Eigen::VectorXd montarVetorForcasReferencia(const Estrutura& est);
};

} // namespace Engine

#endif // ENGINE_ESTRUTURA_H
