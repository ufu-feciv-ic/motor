# Project Overview: IC - Análise Não Linear (Motor)

Este projeto é um "motor" de análise estrutural por Elementos Finitos (MEF) focado em problemas não lineares, especificamente para pórticos 2D (vigas e colunas). Ele utiliza uma formulação corrotacional para lidar com grandes deslocamentos e rotações.

## Main Technologies
- **Language:** C++14
- **Graphics & Windowing:** [Raylib](https://www.raylib.com/)
- **User Interface:** [Dear ImGui](https://github.com/ocornut/imgui)
- **Plotting:** [ImPlot](https://github.com/epezent/implot)
- **Linear Algebra:** [Eigen](https://eigen.tuxfamily.org/)
- **Integration:** [rlImGui](https://github.com/raylib-extras/rlImGui) (Raylib + ImGui bridge)

## Architecture
- **Core Logic (`main.cpp`):** Contém as definições de classes para nós (`No`), materiais (`PropriedadesMaterial`) e elementos finitos (`ElementoFinito`, `Viga2DLinear`, `Viga2DCorrotacional`).
- **Formulação Corrotacional:** Implementada na classe `Viga2DCorrotacional`, permitindo a análise de estruturas com grandes deformações geométricas (P-Delta, etc.).
- **Visualização:** Utiliza Raylib para renderizar a estrutura deformada em tempo real e ImPlot para exibir gráficos de resposta (ex: Força vs. Deslocamento).
- **Dependencies:** 
  - Eigen está incluído em `include/Eigen`.
  - ImGui e ImPlot estão integrados diretamente no diretório raiz.
  - Raylib é esperado no sistema (caminho padrão `C:/raylib/raylib` no Windows).

## Building and Running

### Prerequisites
- **Windows:** [w64devkit](https://github.com/skeeto/w64devkit) ou Mingw-w64 com `make`.
- **Raylib:** Deve estar instalado em `C:/raylib/raylib` ou o caminho deve ser ajustado no `Makefile` ou tarefas do VS Code.

### Build Commands (Windows)
O projeto utiliza um `Makefile` baseado no template do Raylib. No VS Code, existem tarefas pré-configuradas em `.vscode/tasks.json`.

**Manual Build:**
```powershell
# Usando mingw32-make (incluído no w64devkit)
mingw32-make RAYLIB_PATH=C:/raylib/raylib PROJECT_NAME=main OBJS=*.cpp
```

**VS Code:**
- `Ctrl+Shift+B` -> Selecionar `build debug` ou `build release`.

### Running
Após o build, execute o arquivo `main.exe`:
```powershell
./main.exe
```

## Development Conventions

### Coding Style
- O código utiliza uma mistura de Português e Inglês para nomes de variáveis e classes (ex: `No`, `ElementoFinito`, `getGDLsGlobais`).
- Utiliza-se `Eigen` para todas as operações matriciais.
- Preferência por Smart Pointers (`std::shared_ptr`) para gerenciamento de memória de objetos compartilhados como Nós.

### Project Structure
- `include/`: Headers de bibliotecas externas (Eigen, ImGui, ImPlot, Raylib wrappers).
- `src/`: Destinado a arquivos de lógica (atualmente contém o subdiretório `analysis`, que pode ser usado para refatoração do `main.cpp`).
- `resources/`: Assets como fontes (`.ttf`) e imagens.
- Raiz: Contém os arquivos fonte do ImGui, ImPlot e o `main.cpp` principal.

## TODO / Roadmap
- [ ] Refatorar o `main.cpp` movendo as classes de análise para `src/analysis`.
- [ ] Implementar novos tipos de elementos (ex: Treliças, Cascas).
- [ ] Adicionar suporte a outros algoritmos de solução não linear (ex: Newton-Raphson, Comprimento de Arco - Arc-Length).
- [ ] Melhorar a interface de entrada de dados para permitir a criação de modelos dinamicamente.
