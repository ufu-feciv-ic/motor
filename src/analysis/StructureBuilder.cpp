#include "StructureBuilder.h"
#include "Elements.h"
#include <iostream>
#include <unordered_map>

void StructureBuilder::construir(const StructuralModel& model, DOFManager& outDofManager, EstruturaModerna& outEst) {
    outDofManager.limpar();
    outEst.nos.clear();
    outEst.elementos.clear();
    outEst.restricoes.clear();

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
        auto n1 = nodeMap[e.node1_id];
        auto n2 = nodeMap[e.node2_id];
        auto mat = matMap[e.material_id];

        if (e.tipo == "Corrotacional") {
            outEst.elementos.push_back(std::make_shared<Viga2DCorrotacionalModerna>(n1, n2, mat));
        } else if (e.tipo == "Linear") {
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
    for (const auto& load : model.loads) {
        auto indices = dofManager.obterIndicesGlobais(load.node_id);
        if (!indices.empty() && load.dof < (int)indices.size()) {
            F(indices[load.dof]) = load.value;
        }
    }
    return F;
}
