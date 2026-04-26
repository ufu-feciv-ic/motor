#include "SparseAssembler.h"

void SparseAssembler::adicionarContribuicao(
    std::vector<Triplet>& listaTriplets,
    const std::vector<int>& indicesGlobais,
    const Eigen::MatrixXd& matrizLocal) 
{
    int size = static_cast<int>(indicesGlobais.size());
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            listaTriplets.emplace_back(indicesGlobais[i], indicesGlobais[j], matrizLocal(i, j));
        }
    }
}

Eigen::SparseMatrix<double> SparseAssembler::construirMatriz(int totalDofs, const std::vector<Triplet>& listaTriplets) {
    Eigen::SparseMatrix<double> K(totalDofs, totalDofs);
    K.setFromTriplets(listaTriplets.begin(), listaTriplets.end());
    return K;
}
