#pragma once
#include <vector>
#include <map>

/**
 * @brief Gerenciador de Graus de Liberdade (DOFs).
 * Responsável por mapear IDs externos (JSON/UI) para índices internos da matriz global.
 */
class DOFManager {
private:
    // Mapeia ID Externo (ex: do JSON) -> Índice Base na Matriz Global
    std::map<int, int> idExternoParaIndiceBase;
    int proximoIndiceLivre;
    int dofsPorNoDefault;

public:
    DOFManager(int dofsPorNo = 3);

    /**
     * @brief Registra um nó no sistema e reserva seus graus de liberdade.
     */
    void registrarNo(int idExterno, int numDofs = -1);

    /**
     * @brief Retorna os índices globais de um nó para montagem da matriz.
     */
    std::vector<int> obterIndicesGlobais(int idExterno, int numDofs = -1) const;

    int getNumTotalDofs() const;

    void limpar();
};
