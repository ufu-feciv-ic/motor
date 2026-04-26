#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <unordered_map>

#include "eigenpch.h"

#include <raylib.h>
#include <raymath.h>
#include "imgui.h"
#include "implot.h"
#include "rlImGui.h" 

#include "src/analysis/StructuralModel.h"
#include "src/analysis/StructureBuilder.h"
#include "src/analysis/Solver.h"
#include "src/analysis/Elements.h"

struct SimulacaoResultado {
    std::vector<StepResultModerno> history;
    std::vector<double> plot_u;
    std::vector<double> plot_f;
    DOFManager dofManager;
    EstruturaModerna est;
    int dofInteresse;
    double refDim; 
    StructuralModel projetoBase;
};

SimulacaoResultado executarSimulacao(bool naoLinear, int tipoEstrutura) {
    StructuralModel projeto;
    SimulacaoResultado res;
    std::string tipoEl = naoLinear ? "Corrotacional" : "Linear";
    projeto.materials.push_back({1, 1.0, 1.885e6, 9.274e3, "Padrao"});

    if (tipoEstrutura == 0) { // WILLIAMS
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
        res.refDim = 0.368; res.dofInteresse = 5;
    } 
    else if (tipoEstrutura == 1) { // BALANÇO
        for(int i=0; i<6; ++i) projeto.nodes.push_back({i, (double)i * 2.0, 0.0});
        for(int i=0; i<5; ++i) projeto.elements.push_back({i, i, i+1, 1, tipoEl});
        projeto.boundaries.push_back({0, {0, 1, 2}});
        projeto.loads.push_back({5, 1, -500.0}); 
        res.refDim = 1.0; res.dofInteresse = 5;
    }
    else if (tipoEstrutura == 2) { // PORTICO
        projeto.nodes.push_back({0, 0.0, 0.0}); projeto.nodes.push_back({1, 0.0, 5.0});
        projeto.nodes.push_back({2, 5.0, 5.0}); projeto.nodes.push_back({3, 5.0, 0.0});
        projeto.elements.push_back({0, 0, 1, 1, tipoEl}); projeto.elements.push_back({1, 1, 2, 1, tipoEl});
        projeto.elements.push_back({2, 2, 3, 1, tipoEl});
        projeto.boundaries.push_back({0, {0, 1, 2}}); projeto.boundaries.push_back({3, {0, 1, 2}});
        projeto.loads.push_back({1, 0, 100.0}); 
        res.refDim = 1.0; res.dofInteresse = 1;
    }

    res.projetoBase = projeto;
    StructureBuilder::construir(projeto, res.dofManager, res.est);
    Eigen::VectorXd fExt = StructureBuilder::montarVetorCargas(res.dofManager.getNumTotalDofs(), projeto, res.dofManager);
    std::unique_ptr<SolverModerno> solver;
    if (naoLinear) solver = std::make_unique<ArcLengthSolver>(40, 50, 1e-6, (tipoEstrutura==0 ? 0.025 : 0.5), 5.0);
    else solver = std::make_unique<LinearSolver>();
    res.history = solver->solve(res.est, res.dofManager, fExt);

    if (!res.history.empty()) {
        int idxGlobal = res.dofManager.obterIndicesGlobais(res.dofInteresse)[(tipoEstrutura==2 ? 0 : 1)]; 
        for (const auto& step : res.history) {
            res.plot_u.push_back(std::abs(step.u(idxGlobal)) / res.refDim);
            res.plot_f.push_back(step.lambda);
        }
    }
    return res;
}

void DesenharSetaForca(Vector2 pos, float fx, float fy) {
    float mag = std::hypot(fx, fy);
    if (mag < 1e-6) return;
    Vector2 dir = { fx / mag, -fy / mag };
    float comprimentoSeta = 1.5f; float tamanhoCabeca = 0.3f;
    Vector2 tip = pos;
    Vector2 shaftEnd = { pos.x - dir.x * comprimentoSeta, pos.y - dir.y * comprimentoSeta };
    DrawLineEx(shaftEnd, tip, 0.08f, RED);
    float angulo = std::atan2(dir.y, dir.x);
    Vector2 p1 = { tip.x - tamanhoCabeca * std::cos(angulo - 0.5f), tip.y - tamanhoCabeca * std::sin(angulo - 0.5f) };
    Vector2 p2 = { tip.x - tamanhoCabeca * std::cos(angulo + 0.5f), tip.y - tamanhoCabeca * std::sin(angulo + 0.5f) };
    DrawTriangle(tip, p1, p2, RED);
    DrawTriangle(tip, p2, p1, RED);
}

void DesenharDiagramaV2(Vector2 p1, Vector2 p2, double val1, double val2, Color corBorda, Color corPreench) {
    double dx = p2.x - p1.x; double dy = p2.y - p1.y;
    double comp = std::hypot(dx, dy);
    if (comp < 1e-6) return;
    float v1 = (float)val1; float v2 = (float)val2;
    float nx = (float)(-dy / comp); float ny = (float)(dx / comp);
    Vector2 d1 = { p1.x + nx * v1, p1.y + ny * v1 };
    Vector2 d2 = { p2.x + nx * v2, p2.y + ny * v2 };
    auto DrawTri = [](Vector2 v1, Vector2 v2, Vector2 v3, Color col) {
        DrawTriangle(v1, v2, v3, col); DrawTriangle(v1, v3, v2, col);
    };
    if ((v1 > 0 && v2 < 0) || (v1 < 0 && v2 > 0)) {
        float t = std::abs(v1) / (std::abs(v1) + std::abs(v2));
        Vector2 pz = { p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y) };
        DrawTri(p1, d1, pz, corPreench); DrawLineEx(p1, d1, 0.05f, corBorda); DrawLineEx(d1, pz, 0.05f, corBorda);
        DrawTri(pz, d2, p2, corPreench); DrawLineEx(pz, d2, 0.05f, corBorda); DrawLineEx(d2, p2, 0.05f, corBorda);
    } else {
        DrawTri(p1, d1, d2, corPreench); DrawTri(p1, d2, p2, corPreench);
        DrawLineEx(p1, d1, 0.05f, corBorda); DrawLineEx(d1, d2, 0.05f, corBorda); DrawLineEx(d2, p2, 0.05f, corBorda);
    }
}

int main()
{
    bool usarNaoLinear = true;
    int tipoEstrutura = 0; 
    SimulacaoResultado sim = executarSimulacao(usarNaoLinear, tipoEstrutura);

    if (sim.history.empty()) {
        std::cerr << "ERRO FATAL: Falha na simulacao inicial!" << std::endl;
        return 1;
    }

    bool exibirOriginal = true; bool exibirDeformada = true;
    bool diagNaOriginal = false; bool diagNaDeformada = true;
    bool exibirForcas = true;

    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(1280, 720, "Motor de Analise Estrutural - IC UFU");
    
    Camera2D camera = { 0 };
    camera.target = { 13.0f, 0.0f };
    camera.offset = { 1280/2.0f, 720/2.0f + 100.0f };
    camera.rotation = 0.0f;
    camera.zoom = 25.0f;

    SetTargetFPS(60);
    rlImGuiSetup(true);
    ImPlot::CreateContext();

    int current_step = 0;
    int tipoDiagrama = 0;
    float escalaDiagrama = 1.0f;

    while (!WindowShouldClose())
    {
        if (!ImGui::GetIO().WantCaptureMouse) {
            if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
                Vector2 delta = GetMouseDelta();
                delta = Vector2Scale(delta, -1.0f/camera.zoom);
                camera.target = Vector2Add(camera.target, delta);
            }
            float wheel = GetMouseWheelMove();
            if (wheel != 0) {
                Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), camera);
                camera.offset = GetMousePosition();
                camera.target = mouseWorldPos;
                camera.zoom += wheel * 2.0f;
                if (camera.zoom < 1.0f) camera.zoom = 1.0f;
            }
        }

        BeginDrawing();
        ClearBackground(Color{20, 20, 25, 255});
        
        if (current_step >= (int)sim.history.size()) current_step = (int)sim.history.size() - 1;
        if (current_step < 0) current_step = 0;
        const auto& state = sim.history[current_step];

        // Mapa de nós para busca segura por ID
        std::unordered_map<int, std::shared_ptr<NoModerno>> nodeMap;
        for(const auto& n : sim.est.nos) nodeMap[n->id] = n;

        BeginMode2D(camera);
            DrawLineEx({-100, 0}, {100, 0}, 0.02f, DARKGRAY);
            DrawLineEx({0, -100}, {0, 100}, 0.02f, DARKGRAY);

            for (const auto& el : sim.est.elementos) {
                auto [id1, id2] = el->getNodeIds();
                if (nodeMap.find(id1) == nodeMap.end() || nodeMap.find(id2) == nodeMap.end()) continue;

                auto n1 = nodeMap[id1]; auto n2 = nodeMap[id2];
                auto g1 = sim.dofManager.obterIndicesGlobais(id1);
                auto g2 = sim.dofManager.obterIndicesGlobais(id2);

                if (g1.size() < 2 || g2.size() < 2) continue;

                Vector2 p1_0 = { (float)n1->x, (float)-n1->y };
                Vector2 p2_0 = { (float)n2->x, (float)-n2->y };
                Vector2 p1_def = { (float)(n1->x + state.u(g1[0])), (float)-(n1->y + state.u(g1[1])) };
                Vector2 p2_def = { (float)(n2->x + state.u(g2[0])), (float)-(n2->y + state.u(g2[1])) };

                if (exibirOriginal) DrawLineEx(p1_0, p2_0, 0.05f, Color{100, 100, 100, 150});
                if (exibirDeformada) {
                    DrawLineEx(p1_def, p2_def, 0.12f, SKYBLUE);
                    DrawCircleV(p1_def, 0.08f, WHITE); DrawCircleV(p2_def, 0.08f, WHITE);
                }

                if (tipoDiagrama != 0) {
                    EsforcosLocaisModernos esf = el->getEsforcosLocais(state.u, sim.dofManager);
                    double v1 = 0, v2 = 0; Color cb = BLACK, cp = BLANK;
                    if (tipoDiagrama == 1) { v1 = v2 = esf.N * escalaDiagrama; cb = GREEN; cp = ColorAlpha(GREEN, 0.4f); }
                    else if (tipoDiagrama == 2) { v1 = v2 = esf.V * escalaDiagrama; cb = ORANGE; cp = ColorAlpha(ORANGE, 0.4f); }
                    else if (tipoDiagrama == 3) { v1 = -esf.M1 * escalaDiagrama; v2 = esf.M2 * escalaDiagrama; cb = RED; cp = ColorAlpha(RED, 0.4f); }
                    
                    if (diagNaOriginal && exibirOriginal) DesenharDiagramaV2(p1_0, p2_0, v1, v2, cb, cp);
                    if (diagNaDeformada && exibirDeformada) DesenharDiagramaV2(p1_def, p2_def, v1, v2, cb, cp);
                }
            }

            if (exibirForcas) {
                for (const auto& load : sim.projetoBase.loads) {
                    if (nodeMap.find(load.node_id) == nodeMap.end()) continue;
                    auto n = nodeMap[load.node_id];
                    auto g = sim.dofManager.obterIndicesGlobais(load.node_id);
                    if (g.size() < 2) continue;
                    Vector2 pos = { (float)(n->x + state.u(g[0])), (float)-(n->y + state.u(g[1])) };
                    DesenharSetaForca(pos, (load.dof==0?(float)load.value:0), (load.dof==1?(float)load.value:0));
                }
            }
        EndMode2D();

        rlImGuiBegin();
        ImGui::Begin("Painel de Controle");
        const char* itens[] = { "Portico de Williams", "Viga em Balanco", "Portico Simples" };
        if (ImGui::Combo("Modelo", &tipoEstrutura, itens, 3)) { 
            sim = executarSimulacao(usarNaoLinear, tipoEstrutura); current_step = (int)sim.history.size() - 1; 
        }
        if (ImGui::Checkbox("Analise Nao Linear", &usarNaoLinear)) { 
            sim = executarSimulacao(usarNaoLinear, tipoEstrutura); current_step = (int)sim.history.size() - 1; 
        }
        ImGui::Checkbox("Exibir Forcas", &exibirForcas);
        ImGui::Separator();
        ImGui::Checkbox("Estrutura Original", &exibirOriginal);
        ImGui::Checkbox("Estrutura Deformada", &exibirDeformada);
        ImGui::Checkbox("Diagrama na Original", &diagNaOriginal);
        ImGui::Checkbox("Diagrama na Deformada", &diagNaDeformada);
        ImGui::Separator();
        ImGui::SliderInt("Passo", &current_step, 0, (int)sim.history.size() - 1);
        ImGui::SliderFloat("Escala Diagramas", &escalaDiagrama, 0.01f, 10.0f, "%.2f");
        ImGui::Text("Esforço:");
        ImGui::RadioButton("Nenhum", &tipoDiagrama, 0); ImGui::SameLine();
        ImGui::RadioButton("Normal", &tipoDiagrama, 1); ImGui::SameLine();
        ImGui::RadioButton("Cortante", &tipoDiagrama, 2); ImGui::SameLine();
        ImGui::RadioButton("Momento", &tipoDiagrama, 3);
        ImGui::End();

        ImGui::Begin("Resposta Estrutural");
        if (!sim.plot_u.empty() && ImPlot::BeginPlot("Curva de Equilibrio")) {
            ImPlot::SetupAxes("Desl. Normalizado", "Lambda");
            ImPlot::PlotLine("Caminho", sim.plot_u.data(), sim.plot_f.data(), (int)sim.plot_u.size());
            if (current_step < (int)sim.plot_u.size()) {
                double cu = sim.plot_u[current_step], cf = sim.plot_f[current_step];
                ImPlot::PlotScatter("Ponto Atual", &cu, &cf, 1);
            }
            ImPlot::EndPlot();
        }
        ImGui::End();
        rlImGuiEnd();
        EndDrawing();
    }
    ImPlot::DestroyContext(); rlImGuiShutdown(); CloseWindow();
    return 0;
}
