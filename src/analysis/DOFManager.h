#pragma once
#include <vector>

/**
 * @brief Gerenciador de Graus de Liberdade (DOFs).
 * Agora simplificado: assume que o ID do nó define sua posição na matriz global.
 * O mapeamento é direto: indiceBase = idNo * dofsPorNo.
 */
class DOFManager {
private:
    int maiorIdEncontrado;
    int dofsPorNoDefault;

public:
    DOFManager(int dofsPorNo = 3);

    /**
     * @brief Registra um nó para determinar o tamanho total da matriz.
     */
    void registrarNo(int idNo, int numDofs = -1);

    /**
     * @brief Retorna os índices globais de um nó (cálculo direto: id * nDofs).
     */
    std::vector<int> obterIndicesGlobais(int idNo, int numDofs = -1) const;

    int getNumTotalDofs() const;

    void limpar();
};
