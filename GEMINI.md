# Project Overview: IC - Análise Não Linear (Motor)

Este projeto é um "motor" de análise estrutural por Elementos Finitos (MEF) focado em problemas não lineares, especificamente para pórticos 2D. O software foi refatorado para uma arquitetura profissional desacoplada, permitindo alta performance e facilidade de expansão.

## Tecnologias Principais
- **Linguagem:** C++17
- **Álgebra Linear:** [Eigen](https://eigen.tuxfamily.org/) (utilizando Solvers Esparsos)
- **Gráficos & Janelamento:** [Raylib](https://www.raylib.com/)
- **Interface de Usuário:** [Dear ImGui](https://github.com/ocornut/imgui)
- **Plotagem:** [ImPlot](https://github.com/epezent/implot)

## Arquitetura Desacoplada (`src/analysis/`)
A lógica de engenharia está isolada da interface gráfica:

- **`StructuralModel.h`**: Container de dados puros (Nós, Elementos, Materiais, Cargas). Preparado para integração com arquivos JSON.
- **`DOFManager.h/.cpp`**: Gerenciador de Graus de Liberdade. Mapeia IDs arbitrários para índices internos, permitindo geometrias complexas e não sequenciais.
- **`SparseAssembler.h/.cpp`**: Montador de alto desempenho. Utiliza a técnica de *Triplets* para construir matrizes esparsas do Eigen de forma eficiente.
- **`FiniteElementCore.h/.cpp`**: Define as interfaces base e o container da `EstruturaModerna`.
- **`Elements.h/.cpp`**: Implementação física dos elementos:
  - `Viga2DLinearModerna`: Teoria de pequenas deformações.
  - `Viga2DCorrotacionalModerna`: Formulação corrotacional para grandes deslocamentos e rotações.
- **`Solver.h/.cpp`**: Algoritmos de solução:
  - `ArcLengthSolver`: Implementação robusta do método de Comprimento de Arco para traçar curvas de equilíbrio pós-flambagem (Snap-through), utilizando `Eigen::SparseLU`.
- **`StructureBuilder.h/.cpp`**: Orquestrador que transforma dados brutos do modelo em objetos de cálculo prontos para a simulação.

## Funcionalidades Atuais
- **Seleção de Modelos**: Suporte nativo para Pórtico de Williams, Viga em Balanço e Pórtico Simples (P-Delta).
- **Alternância de Teorias**: Possibilidade de comparar resultados Lineares vs. Não Lineares em tempo real.
- **Visualização Avançada**:
  - Renderização da estrutura original e deformada.
  - Diagramas de Esforços Internos (Normal, Cortante e Momento) com preenchimento colorido e tratamento de inversão de sinal.
  - Gráfico de resposta dinâmico (Força vs. Deslocamento) utilizando ImPlot.

## Construção e Execução

### Pré-requisitos
- **Raylib**: Deve estar instalado em `C:/raylib/raylib` (ou ajuste o `RAYLIB_PATH` no Makefile).
- **Compilador**: MinGW-w64 ou w64devkit com suporte a C++17.

### Comandos de Build (Windows)
O Makefile é recursivo e detecta automaticamente arquivos em `src/`.

```powershell
# Compilar projeto completo
mingw32-make RAYLIB_PATH=C:/raylib/raylib PROJECT_NAME=main

# Limpar arquivos objeto
mingw32-make clean
```

## Próximos Passos
- [ ] Implementar `JSONParser` para carregar estruturas de arquivos externos.
- [ ] Adicionar suporte a cargas distribuídas.
- [ ] Implementar novos tipos de elementos (Treliças, Elementos de Mola).
- [ ] Refinar a interface para permitir edição dinâmica de nós e elementos na tela.
