#pragma once
#include <vector>
#include <string>

/**
 * @brief Estruturas de dados puras para representar o modelo estrutural.
 * Estas structs são ideais para serialização/deserialização (JSON).
 */

struct MaterialData {
    int id;
    double E; // Módulo de elasticidade
    double A; // Área da seção
    double I; // Momento de inércia
    std::string nome;
};

struct NodeData {
    int id;
    double x;
    double y;
};

struct ElementData {
    int id;
    int node1_id;
    int node2_id;
    int material_id;
    std::string tipo; // "Linear", "Corrotacional", etc.
};

struct BoundaryCondition {
    int node_id;
    std::vector<int> fixed_dofs; // Ex: {0, 1} para fixar X e Y
};

struct LoadData {
    int node_id;
    int dof; // 0=X, 1=Y, 2=Z_Rot
    double value;
};

/**
 * @brief O "Projeto" completo, pronto para ser construído ou salvo.
 */
struct StructuralModel {
    std::vector<NodeData> nodes;
    std::vector<ElementData> elements;
    std::vector<MaterialData> materials;
    std::vector<BoundaryCondition> boundaries;
    std::vector<LoadData> loads;

    void clear() {
        nodes.clear();
        elements.clear();
        materials.clear();
        boundaries.clear();
        loads.clear();
    }
};
