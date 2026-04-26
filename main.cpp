#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>

#include "eigenpch.h"

#include <raylib.h>
#include "imgui.h"
#include "implot.h"
#include "rlImGui.h" 

#include "src/analysis/StructuralModel.h"
#include "src/analysis/StructureBuilder.h"
#include "src/analysis/Solver.h"
#include "src/analysis/Elements.h"

// Funções utilitárias
Vector2 WorldToScreen(double x, double y, float screenW, float screenH, float scale = 35.0f, float offX = 100.0f, float offY = 0.7f) {
    return { static_cast<float>(offX + x * scale), static_cast<float>(screenH * offY - y * scale) };
}

struct SimulacaoResultado {
    std::vector<StepResultModerno> history;
    std::vector<double> plot_u;
    std::vector<double> plot_f;
    DOFManager dofManager;
    EstruturaModerna est;
    int dofInteresse;
    double refDim; // Dimensão de referência para normalização
};

SimulacaoResultado executarSimulacao(bool naoLinear, int tipoEstrutura) {
    StructuralModel projeto;
    SimulacaoResultado res;
    std::string tipoEl = naoLinear ? "Corrotacional" : "Linear";
    
    // Material padrão (Williams)
    projeto.materials.push_back({1, 1.0, 1.885e6, 9.274e3, "Padrao"});

    if (tipoEstrutura == 0) { // PORTICO DE WILLIAMS
        std::vector<std::pair<double, double>> coords = {
            {0.0, 0.0}, {2.5872, 0.0736}, {5.1744, 0.1472}, {7.7616, 0.2208},
            {10.3488, 0.2944}, {12.936, 0.368}, {15.5232, 0.2944}, {18.1104, 0.2208},
            {20.6976, 0.1472}, {23.2848, 0.0736}, {25.872, 0.0}
        };
        for(int i=0; i<11; ++i) projeto.nodes.push_back({i, coords[i].first, coords[i].second});
        for(int i=0; i<10; ++i) projeto.elements.push_back({i, i, i+1, 1, tipoEl});
        projeto.boundaries.push_back({0, {0, 1, 2}});
        projeto.boundaries.push_back({10, {0, 1, 2}});
        projeto.loads.push_back({5, 1, -1.0});
        res.refDim = 0.368; // Altura do ápice
        res.dofInteresse = 5; // Nó central
    } 
    else if (tipoEstrutura == 1) { // VIGA EM BALANÇO (GRANDES ROTAÇÕES)
        for(int i=0; i<6; ++i) projeto.nodes.push_back({i, (double)i * 2.0, 0.0});
        for(int i=0; i<5; ++i) projeto.elements.push_back({i, i, i+1, 1, tipoEl});
        projeto.boundaries.push_back({0, {0, 1, 2}});
        projeto.loads.push_back({5, 1, -500.0}); // Carga vertical alta
        res.refDim = 1.0; 
        res.dofInteresse = 5;
    }
    else if (tipoEstrutura == 2) { // PORTICO SIMPLES (P-DELTA)
        projeto.nodes.push_back({0, 0.0, 0.0});
        projeto.nodes.push_back({1, 0.0, 5.0});
        projeto.nodes.push_back({2, 5.0, 5.0});
        projeto.nodes.push_back({3, 5.0, 0.0});
        projeto.elements.push_back({0, 0, 1, 1, tipoEl});
        projeto.elements.push_back({1, 1, 2, 1, tipoEl});
        projeto.elements.push_back({2, 2, 3, 1, tipoEl});
        projeto.boundaries.push_back({0, {0, 1, 2}});
        projeto.boundaries.push_back({3, {0, 1, 2}});
        projeto.loads.push_back({1, 0, 100.0}); // Carga horizontal no topo
        res.refDim = 1.0;
        res.dofInteresse = 1;
    }

    StructureBuilder::construir(projeto, res.dofManager, res.est);
    Eigen::VectorXd fExt = StructureBuilder::montarVetorCargas(res.dofManager.getNumTotalDofs(), projeto, res.dofManager);
    
    ArcLengthSolver solver(40, 50, 1e-6, (tipoEstrutura==0 ? 0.025 : 0.5), 5.0);
    res.history = solver.solve(res.est, res.dofManager, fExt);

    int idxGlobal = res.dofManager.obterIndicesGlobais(res.dofInteresse)[(tipoEstrutura==2 ? 0 : 1)]; 
    for (const auto& step : res.history) {
        res.plot_u.push_back(std::abs(step.u(idxGlobal)) / res.refDim);
        res.plot_f.push_back(step.lambda);
    }
    return res;
}

int main()
{
    bool usarNaoLinear = true;
    int tipoEstrutura = 0; // 0: Williams, 1: Balanço, 2: Pórtico
    SimulacaoResultado sim = executarSimulacao(usarNaoLinear, tipoEstrutura);

    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(1280, 720, "Motor de Analise Estrutural - IC UFU");
    SetTargetFPS(60);

    rlImGuiSetup(true);
    ImPlot::CreateContext();

    int current_step = 0;
    int tipoDiagrama = 0;
    float escalaDiagrama = 0.002f;

    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(Color{30, 30, 35, 255});

        const auto& state = sim.history[current_step];

        // RENDERIZAÇÃO
        for (size_t i = 0; i < sim.est.elementos.size(); ++i) {
            auto el = sim.est.elementos[i];
            auto indices1 = el->getIndicesGlobais(sim.dofManager); // Pega os 6 índices
            
            // Precisamos encontrar os nós n1 e n2 do elemento para as coordenadas originais
            // Vamos usar uma abordagem genérica buscando pelos IDs via DOFManager (simplificado aqui)
            // No futuro, o elemento deve guardar seus IDs de nós de forma acessível.
            
            // Hack temporário para visualização: assume que nós estão em ordem na lista da estrutura
            // que condiz com a ordem de criação no executarSimulacao
            int idxN1 = (tipoEstrutura == 2 && i == 2) ? 2 : (int)i;
            int idxN2 = (tipoEstrutura == 2 && i == 2) ? 3 : (int)i + 1;
            if(tipoEstrutura == 2 && i == 1) { idxN1 = 1; idxN2 = 2; }

            auto n1 = sim.est.nos[idxN1];
            auto n2 = sim.est.nos[idxN2];
            auto g1 = sim.dofManager.obterIndicesGlobais(n1->id);
            auto g2 = sim.dofManager.obterIndicesGlobais(n2->id);

            double x1_def = n1->x + state.u(g1[0]);
            double y1_def = n1->y + state.u(g1[1]);
            double x2_def = n2->x + state.u(g2[0]);
            double y2_def = n2->y + state.u(g2[1]);

            float vScale = (tipoEstrutura == 0 ? 35.0f : 50.0f);
            Vector2 p1 = WorldToScreen(x1_def, y1_def, (float)GetScreenWidth(), (float)GetScreenHeight(), vScale);
            Vector2 p2 = WorldToScreen(x2_def, y2_def, (float)GetScreenWidth(), (float)GetScreenHeight(), vScale);

            DrawLineEx(p1, p2, 3.5f, SKYBLUE);
            DrawCircleV(p1, 4.0f, WHITE);
            DrawCircleV(p2, 4.0f, WHITE);

            if (tipoDiagrama != 0) {
                double dx = x2_def - x1_def; double dy = y2_def - y1_def;
                double comp = std::hypot(dx, dy);
                double nx = -dy / comp; double ny = dx / comp;  

                EsforcosLocaisModernos esf = el->getEsforcosLocais(state.u, sim.dofManager);
                double val1 = 0.0, val2 = 0.0;
                Color corBorda = BLACK; Color corPreench = BLANK;

                if (tipoDiagrama == 1) { val1 = val2 = esf.N * escalaDiagrama; corBorda = GREEN; corPreench = ColorAlpha(GREEN, 0.4f); }
                else if (tipoDiagrama == 2) { val1 = val2 = esf.V * escalaDiagrama; corBorda = ORANGE; corPreench = ColorAlpha(ORANGE, 0.4f); }
                else if (tipoDiagrama == 3) { val1 = esf.M1 * escalaDiagrama; val2 = -esf.M2 * escalaDiagrama; corBorda = RED; corPreench = ColorAlpha(RED, 0.4f); }

                Vector2 d1 = WorldToScreen(x1_def + nx * val1, y1_def + ny * val1, (float)GetScreenWidth(), (float)GetScreenHeight(), vScale);
                Vector2 d2 = WorldToScreen(x2_def + nx * val2, y2_def + ny * val2, (float)GetScreenWidth(), (float)GetScreenHeight(), vScale);
                auto DrawTri = [](Vector2 v1, Vector2 v2, Vector2 v3, Color col) { DrawTriangle(v1, v2, v3, col); DrawTriangle(v1, v3, v2, col); };
                if ((val1 > 0 && val2 < 0) || (val1 < 0 && val2 > 0)) {
                    float t = (float)(std::abs(val1) / (std::abs(val1) + std::abs(val2)));
                    Vector2 pz = { p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y) };
                    DrawTri(p1, d1, pz, corPreench); DrawLineEx(p1, d1, 2.0f, corBorda); DrawLineEx(d1, pz, 2.0f, corBorda);
                    DrawTri(pz, d2, p2, corPreench); DrawLineEx(pz, d2, 2.0f, corBorda); DrawLineEx(d2, p2, 2.0f, corBorda);
                } else {
                    DrawTri(p1, d1, d2, corPreench); DrawTri(p1, d2, p2, corPreench);
                    DrawLineEx(p1, d1, 2.0f, corBorda); DrawLineEx(d1, d2, 2.0f, corBorda); DrawLineEx(d2, p2, 2.0f, corBorda);
                }
            }
        }

        rlImGuiBegin();
        ImGui::Begin("Painel de Controle");
        
        ImGui::Text("Estrutura:");
        const char* itens[] = { "Portico de Williams", "Viga em Balanco", "Portico Simples" };
        if (ImGui::Combo("Modelo", &tipoEstrutura, itens, 3)) {
            sim = executarSimulacao(usarNaoLinear, tipoEstrutura);
            current_step = 0;
        }

        if (ImGui::Checkbox("Analise Nao Linear (Corrotacional)", &usarNaoLinear)) {
            sim = executarSimulacao(usarNaoLinear, tipoEstrutura);
            current_step = 0;
        }

        ImGui::Separator();
        ImGui::SliderInt("Passo", &current_step, 0, (int)sim.history.size() - 1);
        ImGui::SliderFloat("Escala Diagrama", &escalaDiagrama, 0.0001f, 0.05f, "%.5f");
        
        ImGui::Text("Diagramas:");
        ImGui::RadioButton("Nenhum", &tipoDiagrama, 0); ImGui::SameLine();
        ImGui::RadioButton("Normal", &tipoDiagrama, 1); ImGui::SameLine();
        ImGui::RadioButton("Cortante", &tipoDiagrama, 2); ImGui::SameLine();
        ImGui::RadioButton("Momento", &tipoDiagrama, 3);

        ImGui::Separator();
        ImGui::Text("Fator Carga (Lambda): %.4f", state.lambda);
        ImGui::End();

        ImGui::Begin("Grafico de Resposta");
        if (ImPlot::BeginPlot("Curva de Equilibrio", ImVec2(-1, -1))) {
            ImPlot::SetupAxes("Deslocamento (Normalizado)", "Fator de Carga (Lambda)");
            ImPlot::PlotLine("Caminho", sim.plot_u.data(), sim.plot_f.data(), (int)sim.plot_u.size());
            double cu = sim.plot_u[current_step], cf = sim.plot_f[current_step];
            ImPlot::PlotScatter("Ponto Atual", &cu, &cf, 1);
            ImPlot::EndPlot();
        }
        ImGui::End();
        rlImGuiEnd();
        EndDrawing();
    }

    ImPlot::DestroyContext();
    rlImGuiShutdown();
    CloseWindow();
    return 0;
}
