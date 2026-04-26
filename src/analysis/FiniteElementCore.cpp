#include "FiniteElementCore.h"

void EstruturaModerna::aplicarCondicoesContorno(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const DOFManager& dofManager) const {
    for (auto const& [nodeId, dofs] : restricoes) {
        auto indices = dofManager.obterIndicesGlobais(nodeId);
        if (indices.empty()) continue;

        for (int localDof : dofs) {
            if (localDof >= (int)indices.size()) continue;
            int globalIdx = indices[localDof];
            
            K.prune([globalIdx](int row, int col, double value) {
                return row != globalIdx; 
            });
            K.coeffRef(globalIdx, globalIdx) = 1.0;
            F(globalIdx) = 0.0;
        }
    }
}
