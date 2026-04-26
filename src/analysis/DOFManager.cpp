#include "DOFManager.h"

DOFManager::DOFManager(int dofsPorNo) 
    : maiorIdEncontrado(-1), dofsPorNoDefault(dofsPorNo) {}

void DOFManager::registrarNo(int idNo, int /*numDofs*/) {
    if (idNo > maiorIdEncontrado) {
        maiorIdEncontrado = idNo;
    }
}

std::vector<int> DOFManager::obterIndicesGlobais(int idNo, int numDofs) const {
    int nDofs = (numDofs == -1) ? dofsPorNoDefault : numDofs;
    std::vector<int> indices;
    indices.reserve(nDofs);
    
    int base = idNo * nDofs;
    for (int i = 0; i < nDofs; ++i) {
        indices.push_back(base + i);
    }
    return indices;
}

int DOFManager::getNumTotalDofs() const {
    if (maiorIdEncontrado == -1) return 0;
    return (maiorIdEncontrado + 1) * dofsPorNoDefault;
}

void DOFManager::limpar() {
    maiorIdEncontrado = -1;
}
