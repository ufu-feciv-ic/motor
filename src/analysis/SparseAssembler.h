#pragma once
#include <vector>
#include <Eigen/Sparse>

/**
 * @brief Classe utilitária para montagem de matrizes e vetores esparsos.
 */
class SparseAssembler {
public:
    typedef Eigen::Triplet<double> Triplet;

    static void adicionarContribuicao(
        std::vector<Triplet>& listaTriplets,
        const std::vector<int>& indicesGlobais,
        const Eigen::MatrixXd& matrizLocal);

    static Eigen::SparseMatrix<double> construirMatriz(int totalDofs, const std::vector<Triplet>& listaTriplets);
};
