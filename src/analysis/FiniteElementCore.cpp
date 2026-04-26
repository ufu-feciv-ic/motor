#include "FiniteElementCore.h"
#include <iostream>
#include <vector>

void EstruturaModerna::aplicarCondicoesContorno(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const DOFManager& dofManager) const {
    int totalDofs = (int)F.size();
    if (totalDofs <= 0) return;

    // 1. Identificar todos os DOFs fixos de forma eficiente
    std::vector<bool> isFixed(totalDofs, false);
    for (auto const& [nodeId, dofs] : restricoes) {
        auto indices = dofManager.obterIndicesGlobais(nodeId);
        for (int d : dofs) {
            if (d >= 0 && d < (int)indices.size()) {
                int idx = indices[d];
                if (idx >= 0 && idx < totalDofs) isFixed[idx] = true;
            }
        }
    }

    // 2. Zerar linhas e colunas dos DOFs fixos na matriz esparsa
    // Itera apenas sobre os elementos não-zero existentes
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            if (isFixed[it.row()] || isFixed[it.col()]) {
                it.valueRef() = 0.0;
            }
        }
    }

    // 3. Setar diagonal como 1.0 e zerar o termo de carga correspondente
    for (int i = 0; i < totalDofs; ++i) {
        if (isFixed[i]) {
            K.coeffRef(i, i) = 1.0;
            F(i) = 0.0;
        }
    }

    K.makeCompressed();
}
