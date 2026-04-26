#include "StructureBuilder.h"
#include "Elements.h"
#include <iostream>
#include <unordered_map>

void StructureBuilder::construir(const StructuralModel& model, DOFManager& outDofManager, EstruturaModerna& outEst) {
    outDofManager.limpar();
    outEst.nos.clear();
    outEst.elementos.clear();
    outEst.restricoes.clear();

    if (model.nodes.empty()) return;

    // 1. Criar e registrar nós
    std::unordered_map<int, std::shared_ptr<NoModerno>> nodeMap;
    for (const auto& n : model.nodes) {
        auto no = std::make_shared<NoModerno>(n.id, n.x, n.y);
        outEst.nos.push_back(no);
        nodeMap[n.id] = no;
        outDofManager.registrarNo(n.id);
    }

    // 2. Mapeamento de Materiais
    std::unordered_map<int, MaterialData> matMap;
    for (const auto& m : model.materials) matMap[m.id] = m;

    // 3. Criar Elementos
    for (const auto& e : model.elements) {
        if (nodeMap.find(e.node1_id) == nodeMap.end() || nodeMap.find(e.node2_id) == nodeMap.end()) continue;
        
        auto n1 = nodeMap[e.node1_id];
        auto n2 = nodeMap[e.node2_id];
        
        if (matMap.find(e.material_id) == matMap.end()) continue;
        auto mat = matMap[e.material_id];

        if (e.tipo == "Corrotacional") {
            outEst.elementos.push_back(std::make_shared<Viga2DCorrotacionalModerna>(n1, n2, mat));
        } else {
            outEst.elementos.push_back(std::make_shared<Viga2DLinearModerna>(n1, n2, mat));
        }
    }

    // 4. Condições de Contorno
    for (const auto& bc : model.boundaries) {
        outEst.restricoes[bc.node_id] = bc.fixed_dofs;
    }
}

Eigen::VectorXd StructureBuilder::montarVetorCargas(int totalDofs, const StructuralModel& model, const DOFManager& dofManager) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(totalDofs);
    if (totalDofs <= 0) return F;

    for (const auto& load : model.loads) {
        auto indices = dofManager.obterIndicesGlobais(load.node_id);
        if (!indices.empty() && load.dof < (int)indices.size()) {
            int idxGlobal = indices[load.dof];
            if (idxGlobal >= 0 && idxGlobal < totalDofs) {
                F(idxGlobal) = load.value;
            }
        }
    }
    return F;
}
