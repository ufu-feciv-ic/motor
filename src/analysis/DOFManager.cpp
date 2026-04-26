#include "DOFManager.h"
#include <iostream>

DOFManager::DOFManager(int dofsPorNo) 
    : proximoIndiceLivre(0), dofsPorNoDefault(dofsPorNo) {}

void DOFManager::registrarNo(int idExterno, int numDofs) {
    if (idExternoParaIndiceBase.find(idExterno) != idExternoParaIndiceBase.end()) {
        return; // Já registrado
    }

    int nDofs = (numDofs == -1) ? dofsPorNoDefault : numDofs;
    idExternoParaIndiceBase[idExterno] = proximoIndiceLivre;
    proximoIndiceLivre += nDofs;
}

std::vector<int> DOFManager::obterIndicesGlobais(int idExterno, int numDofs) const {
    auto it = idExternoParaIndiceBase.find(idExterno);
    if (it == idExternoParaIndiceBase.end()) {
        std::cerr << "ERRO: No ID " << idExterno << " nao foi registrado no DOFManager!" << std::endl;
        return {};
    }

    int nDofs = (numDofs == -1) ? dofsPorNoDefault : numDofs;
    std::vector<int> indices;
    indices.reserve(nDofs);
    
    int base = it->second;
    for (int i = 0; i < nDofs; ++i) {
        indices.push_back(base + i);
    }
    return indices;
}

int DOFManager::getNumTotalDofs() const {
    return proximoIndiceLivre;
}

void DOFManager::limpar() {
    idExternoParaIndiceBase.clear();
    proximoIndiceLivre = 0;
}
