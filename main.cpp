#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>

#include "Engine/Common.h"
#include "Engine/No.h"
#include "Engine/ElementoFinito.h"
#include "Engine/Estrutura.h"
#include "Engine/Solvers.h"

#include <raylib.h>
#include "imgui.h"
#include "implot.h"
#include "rlImGui.h" 

using namespace Engine;

struct ResultadoPassoUI
{
    Eigen::VectorXd udesl;   
    Eigen::VectorXd reacoes; 
    double lambda;           
    double u_apex;           
    double f_apex;           
};

struct ArestaRender
{
    int n1;
    int n2;
};

// --- FUNÇÕES DE CONFIGURAÇÃO DE MODELOS ---

void SetupWilliams(Estrutura& est, std::vector<ArestaRender>& arestas, double& h_apex, int& dofApexY, int analiseTipo) {
    est = Estrutura();
    arestas.clear();
    PropriedadesMaterial mat = {1.0, 1.885e6, 9.274e3};
    std::vector<std::pair<double, double>> coords = {
        {0.0, 0.0}, {2.5872, 0.0736}, {5.1744, 0.1472}, {7.7616, 0.2208},
        {10.3488, 0.2944}, {12.936, 0.368}, {15.5232, 0.2944}, {18.1104, 0.2208},
        {20.6976, 0.1472}, {23.2848, 0.0736}, {25.872, 0.0}
    };
    int dof_count = 0;
    for (int i = 0; i < 11; ++i) {
        auto no = std::make_shared<No>(i, coords[i].first, coords[i].second, std::vector<int>{dof_count, dof_count + 1, dof_count + 2});
        est.adicionarNo(no);
        dof_count += 3;
    }
    for (int i = 0; i < 10; ++i) {
        std::shared_ptr<ElementoFinito> viga;
        if (analiseTipo == 0) viga = std::make_shared<Viga2DLinear>(est.Nos[i], est.Nos[i+1], mat);
        else viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i + 1], mat);
        est.adicionarElemento(viga);
        arestas.push_back({i, i + 1});
    }
    est.NosFixos = {0, 1, 2, 30, 31, 32};
    est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    est.ForcasExternas(16) = -1.0;
    h_apex = 0.368;
    dofApexY = 16;
}

void SetupPilarBiengastado(Estrutura& est, std::vector<ArestaRender>& arestas, double& h_apex, int& dofApexY, int analiseTipo) {
    est = Estrutura();
    arestas.clear();
    PropriedadesMaterial mat = {210e9, 0.01, 0.0001};
    double L = 10.0;
    int nElem = 4;
    int dof_count = 0;
    for (int i = 0; i <= nElem; ++i) {
        auto no = std::make_shared<No>(i, 0.0, (L / nElem) * i, std::vector<int>{dof_count, dof_count + 1, dof_count + 2});
        est.adicionarNo(no);
        dof_count += 3;
    }
    for (int i = 0; i < nElem; ++i) {
        std::shared_ptr<ElementoFinito> viga;
        if (analiseTipo == 0) viga = std::make_shared<Viga2DLinear>(est.Nos[i], est.Nos[i+1], mat);
        else viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i + 1], mat);
        est.adicionarElemento(viga);
        arestas.push_back({i, i + 1});
    }
    est.NosFixos = {0, 1, 2, 12, 13, 14}; 
    est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    est.ForcasExternas(6) = 1000.0; 
    h_apex = 1.0;
    dofApexY = 6; 
}

void SetupVigaEngastada(Estrutura& est, std::vector<ArestaRender>& arestas, double& h_apex, int& dofApexY, int analiseTipo) {
    est = Estrutura();
    arestas.clear();
    PropriedadesMaterial mat = {210e9, 0.01, 0.0001};
    double L = 10.0;
    int nElem = 4;
    int dof_count = 0;
    for (int i = 0; i <= nElem; ++i) {
        auto no = std::make_shared<No>(i, (L / nElem) * i, 0.0, std::vector<int>{dof_count, dof_count + 1, dof_count + 2});
        est.adicionarNo(no);
        dof_count += 3;
    }
    for (int i = 0; i < nElem; ++i) {
        std::shared_ptr<ElementoFinito> viga;
        if (analiseTipo == 0) viga = std::make_shared<Viga2DLinear>(est.Nos[i], est.Nos[i+1], mat);
        else viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i + 1], mat);
        est.adicionarElemento(viga);
        arestas.push_back({i, i + 1});
    }
    est.NosFixos = {0, 1, 2}; 
    est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    est.ForcasExternas(13) = -1000.0; 
    h_apex = 1.0;
    dofApexY = 13;
}

void SetupVigaBiapoiada(Estrutura& est, std::vector<ArestaRender>& arestas, double& h_apex, int& dofApexY, int analiseTipo) {
    est = Estrutura();
    arestas.clear();
    PropriedadesMaterial mat = {210e9, 0.01, 0.0001};
    double L = 10.0;
    int nElem = 4;
    int dof_count = 0;
    for (int i = 0; i <= nElem; ++i) {
        auto no = std::make_shared<No>(i, (L / nElem) * i, 0.0, std::vector<int>{dof_count, dof_count + 1, dof_count + 2});
        est.adicionarNo(no);
        dof_count += 3;
    }
    for (int i = 0; i < nElem; ++i) {
        std::shared_ptr<ElementoFinito> viga;
        if (analiseTipo == 0) viga = std::make_shared<Viga2DLinear>(est.Nos[i], est.Nos[i+1], mat);
        else viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i + 1], mat);
        est.adicionarElemento(viga);
        arestas.push_back({i, i + 1});
    }
    est.NosFixos = {0, 1, 13}; 
    est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    est.ForcasExternas(7) = -1000.0; 
    h_apex = 1.0;
    dofApexY = 7;
}

void SetupPilarEngastado(Estrutura& est, std::vector<ArestaRender>& arestas, double& h_apex, int& dofApexY, int analiseTipo) {
    est = Estrutura();
    arestas.clear();
    PropriedadesMaterial mat = {210e9, 0.01, 0.0001};
    double L = 10.0;
    int nElem = 10;
    int dof_count = 0;
    for (int i = 0; i <= nElem; ++i) {
        auto no = std::make_shared<No>(i, 0.0, (L / nElem) * i, std::vector<int>{dof_count, dof_count + 1, dof_count + 2});
        est.adicionarNo(no);
        dof_count += 3;
    }
    for (int i = 0; i < nElem; ++i) {
        std::shared_ptr<ElementoFinito> viga;
        if (analiseTipo == 0) viga = std::make_shared<Viga2DLinear>(est.Nos[i], est.Nos[i+1], mat);
        else viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i + 1], mat);
        est.adicionarElemento(viga);
        arestas.push_back({i, i + 1});
    }
    est.NosFixos = {0, 1, 2}; // Engastado na base
    est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    est.ForcasExternas(nElem * 3 + 1) = -1000.0; // Carga axial no topo
    h_apex = L;
    dofApexY = nElem * 3; // Deslocamento X no topo para monitorar flambagem
}

void SetupVigaDistribuida(Estrutura& est, std::vector<ArestaRender>& arestas, double& h_apex, int& dofApexY, int analiseTipo) {
    est = Estrutura();
    arestas.clear();
    PropriedadesMaterial mat = {210e9, 0.01, 0.0001};
    double L = 10.0;
    int nElem = 10;
    int dof_count = 0;
    for (int i = 0; i <= nElem; ++i) {
        auto no = std::make_shared<No>(i, (L / nElem) * i, 0.0, std::vector<int>{dof_count, dof_count + 1, dof_count + 2});
        est.adicionarNo(no);
        dof_count += 3;
    }
    for (int i = 0; i < nElem; ++i) {
        std::shared_ptr<ElementoFinito> viga;
        if (analiseTipo == 0) viga = std::make_shared<Viga2DLinear>(est.Nos[i], est.Nos[i+1], mat);
        else viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i + 1], mat);
        
        // Aplicando carga distribuída transversal de 1000 N/m
        viga->setCargaDistribuida(0.0, -1000.0);
        
        est.adicionarElemento(viga);
        arestas.push_back({i, i + 1});
    }
    // Bi-apoiada (Pinned-Roller): No 0 fixo em X e Y, Ultimo No fixo em Y
    est.NosFixos = {
        est.Nos.front()->gdlGlobais[0], 
        est.Nos.front()->gdlGlobais[1], 
        est.Nos.back()->gdlGlobais[1]
    };
    est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    h_apex = 1.0;
    dofApexY = est.Nos[nElem / 2]->gdlGlobais[1]; // Centro da viga
}

void DrawArrow(Vector2 start, Vector2 end, float thickness, Color color) {
    float dx = end.x - start.x;
    float dy = end.y - start.y;
    float length = sqrtf(dx*dx + dy*dy);
    if (length < 0.1f) return;

    DrawLineEx(start, end, thickness, color);
    
    float angle = atan2f(dy, dx);
    float headLen = 15.0f;
    float headAngle = 0.4f; 
    
    Vector2 p1 = { end.x - headLen * cosf(angle - headAngle), end.y - headLen * sinf(angle - headAngle) };
    Vector2 p2 = { end.x - headLen * cosf(angle + headAngle), end.y - headLen * sinf(angle + headAngle) };
    
    DrawTriangle(end, p1, p2, color);
}

Vector2 WorldToScreen(double x, double y, float screenW, float screenH) {
    float scale = 30.0f; 
    float offsetX = screenW * 0.2f;
    float offsetY = screenH * 0.6f; 
    return { static_cast<float>(offsetX + x * scale), static_cast<float>(offsetY - y * scale) };
}

int main()
{
    // Variáveis de Controle
    int modeloSelecionado = 0; 
    int analiseSelecionada = 1; 
    
    // Imperfeição Inicial
    bool usarImperfeicao = false;
    int modoImperfeicao = 0;
    float amplitudeImperfeicao = 0.05f;
    std::vector<AnaliseBuckling::ModoBuckling> modosBuckling;

    Estrutura est;
    std::vector<ArestaRender> arestas;
    double h_apex = 0.368;
    int dofApexY = 16;
    std::vector<Resultado> history;
    std::vector<ResultadoPassoUI> historyUI;
    std::vector<double> plot_u, plot_f;
    int current_step = 0;
    int max_steps = 0;
    int tipoDiagrama = 0; 
    float escalaDiagrama = 0.001f;
    bool mostrarReacoes = false;
    float escalaReacoes = 0.001f;

    auto ExecutarSimulacao = [&]() {
        if (modeloSelecionado == 0) SetupWilliams(est, arestas, h_apex, dofApexY, analiseSelecionada);
        else if (modeloSelecionado == 1) SetupPilarBiengastado(est, arestas, h_apex, dofApexY, analiseSelecionada);
        else if (modeloSelecionado == 2) SetupVigaEngastada(est, arestas, h_apex, dofApexY, analiseSelecionada);
        else if (modeloSelecionado == 3) SetupVigaBiapoiada(est, arestas, h_apex, dofApexY, analiseSelecionada);
        else if (modeloSelecionado == 4) SetupPilarEngastado(est, arestas, h_apex, dofApexY, analiseSelecionada);
        else if (modeloSelecionado == 5) SetupVigaDistribuida(est, arestas, h_apex, dofApexY, analiseSelecionada);

        modosBuckling = AnaliseBuckling::executar(est);

        if (usarImperfeicao && !modosBuckling.empty()) {
            if (modoImperfeicao >= (int)modosBuckling.size()) modoImperfeicao = 0;
            est.perturbarGeometria(modosBuckling[modoImperfeicao].modo, (double)amplitudeImperfeicao);
        }

        if (analiseSelecionada == 0) {
            AnaliseLinear solver;
            history = solver.executar(est);
            std::vector<Resultado> h_lin;
            h_lin.push_back({Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), 0.0});
            h_lin.push_back(history[0]);
            history = h_lin;
        } else {
            if (modeloSelecionado == 0) {
                AnaliseNaoLinearCompArco solver(33, 50, 1e-6, 0.025, 5.0);
                history = solver.executar(est);
            } else {
                AnaliseNaoLinearCompArco solver(30, 50, 1e-6, 0.2, 5.0);
                history = solver.executar(est);
            }
        }

        historyUI.clear();
        for (const auto& step : history) {
            ResultadoPassoUI s;
            s.udesl = step.u;
            s.reacoes = step.reacoes;
            s.lambda = step.FatorCarga;
            s.u_apex = std::abs(step.u(dofApexY)) / h_apex;
            s.f_apex = std::abs(step.FatorCarga);
            historyUI.push_back(s);
        }
        plot_u.clear(); plot_f.clear();
        for (const auto& step : historyUI) {
            plot_u.push_back(step.u_apex);
            plot_f.push_back(step.f_apex);
        }
        current_step = (int)historyUI.size() - 1;
        max_steps = (int)historyUI.size() - 1;
    };

    ExecutarSimulacao();

    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(1280, 720, "Análise Estrutural 2D - Raylib + ImGui");
    SetTargetFPS(60);

    rlImGuiSetup(true);
    ImPlot::CreateContext();

    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(Color{30, 30, 35, 255});

        if (!historyUI.empty()) {
            const ResultadoPassoUI& state = historyUI[current_step];

            for (const auto& e : arestas) {
                const auto& no1 = est.Nos[e.n1];
                const auto& no2 = est.Nos[e.n2];
                Vector2 p1 = WorldToScreen(no1->x + state.udesl(no1->gdlGlobais[0]), no1->y + state.udesl(no1->gdlGlobais[1]), (float)GetScreenWidth(), (float)GetScreenHeight());
                Vector2 p2 = WorldToScreen(no2->x + state.udesl(no2->gdlGlobais[0]), no2->y + state.udesl(no2->gdlGlobais[1]), (float)GetScreenWidth(), (float)GetScreenHeight());
                DrawLineEx(p1, p2, 3.0f, SKYBLUE);
                DrawCircleV(p1, 4.0f, WHITE);
                DrawCircleV(p2, 4.0f, WHITE);
            }

            // Desenhar Cargas Distribuídas
            for (size_t i = 0; i < arestas.size(); ++i) {
                const auto& elem = est.Elementos[i];
                if (std::abs(elem->carga.qx) > 1e-6 || std::abs(elem->carga.qy) > 1e-6) {
                    const auto& e = arestas[i];
                    const auto& n1 = est.Nos[e.n1];
                    const auto& n2 = est.Nos[e.n2];
                    
                    Vector2 p1 = WorldToScreen(n1->x + state.udesl(n1->gdlGlobais[0]), n1->y + state.udesl(n1->gdlGlobais[1]), (float)GetScreenWidth(), (float)GetScreenHeight());
                    Vector2 p2 = WorldToScreen(n2->x + state.udesl(n2->gdlGlobais[0]), n2->y + state.udesl(n2->gdlGlobais[1]), (float)GetScreenWidth(), (float)GetScreenHeight());
                    
                    double dx = (n2->x + state.udesl(n2->gdlGlobais[0])) - (n1->x + state.udesl(n1->gdlGlobais[0]));
                    double dy = (n2->y + state.udesl(n2->gdlGlobais[1])) - (n1->y + state.udesl(n1->gdlGlobais[1]));
                    double L = std::hypot(dx, dy);
                    if (L < 1e-6) L = 1e-6;
                    double nx = -dy / L;
                    double ny = dx / L;
                    double tx = dx / L;
                    double ty = dy / L;

                    int numArrows = 3;
                    for (int j = 0; j <= numArrows; ++j) {
                        float t = (float)j / (float)numArrows;
                        Vector2 p = { p1.x + (p2.x - p1.x) * t, p1.y + (p2.y - p1.y) * t };
                        
                        double loadX = (elem->carga.qx * tx + elem->carga.qy * nx) * state.lambda;
                        double loadY = (elem->carga.qx * ty + elem->carga.qy * ny) * state.lambda;
                        
                        float vx = (float)loadX;
                        float vy = (float)(-loadY);
                        float mag = sqrtf(vx*vx + vy*vy);
                        if (mag < 1e-6) continue;
                        vx /= mag; vy /= mag;

                        float arrowLen = 30.0f;
                        Vector2 pStart = { p.x - vx * arrowLen, p.y - vy * arrowLen };
                        DrawArrow(pStart, p, 2.0f, ORANGE);
                    }
                }
            }

            for (const auto& no : est.Nos) {
                double fx = est.ForcasExternas(no->gdlGlobais[0]) * state.lambda;
                double fy = est.ForcasExternas(no->gdlGlobais[1]) * state.lambda;

                if (std::abs(fx) > 1e-6 || std::abs(fy) > 1e-6) {
                    Vector2 pNode = WorldToScreen(no->x + state.udesl(no->gdlGlobais[0]), 
                                                  no->y + state.udesl(no->gdlGlobais[1]), 
                                                  (float)GetScreenWidth(), (float)GetScreenHeight());
                    
                    float vx = (float)fx;
                    float vy = (float)(-fy); 
                    float mag = sqrtf(vx*vx + vy*vy);
                    vx /= mag; vy /= mag;

                    float arrowLen = 50.0f;
                    Vector2 pStart = { pNode.x - vx * arrowLen, pNode.y - vy * arrowLen };
                    DrawArrow(pStart, pNode, 3.0f, RED);
                    DrawText(TextFormat("%.2f", std::sqrt(fx*fx + fy*fy)), (int)pStart.x - 10, (int)pStart.y - 20, 18, YELLOW);
                }
            }

            if (mostrarReacoes) {
                for (const auto& no : est.Nos) {
                    double rx = state.reacoes(no->gdlGlobais[0]);
                    double ry = state.reacoes(no->gdlGlobais[1]);
                    double mz = state.reacoes(no->gdlGlobais[2]);

                    if (std::abs(rx) > 1e-4 || std::abs(ry) > 1e-4) {
                        Vector2 pNode = WorldToScreen(no->x + state.udesl(no->gdlGlobais[0]), 
                                                      no->y + state.udesl(no->gdlGlobais[1]), 
                                                      (float)GetScreenWidth(), (float)GetScreenHeight());
                        
                        float vx = (float)rx;
                        float vy = (float)(-ry); 
                        float mag = sqrtf(vx*vx + vy*vy);
                        vx /= mag; vy /= mag;

                        float arrowLen = 40.0f;
                        Vector2 pStart = { pNode.x - vx * arrowLen, pNode.y - vy * arrowLen };
                        DrawArrow(pStart, pNode, 2.5f, MAGENTA);
                        DrawText(TextFormat("R:%.1f", mag), (int)pStart.x + 5, (int)pStart.y + 5, 15, MAGENTA);
                    }
                    
                    if (std::abs(mz) > 1e-4) {
                         Vector2 pNode = WorldToScreen(no->x + state.udesl(no->gdlGlobais[0]), 
                                                      no->y + state.udesl(no->gdlGlobais[1]), 
                                                      (float)GetScreenWidth(), (float)GetScreenHeight());
                         DrawCircleLines((int)pNode.x, (int)pNode.y, 15.0f, MAGENTA);
                         DrawText(TextFormat("M:%.1f", mz), (int)pNode.x - 20, (int)pNode.y + 20, 15, MAGENTA);
                    }
                }
            }

            if (tipoDiagrama != 0) {
                for (size_t i = 0; i < arestas.size(); ++i) {
                    const auto& e = arestas[i];
                    const auto& n1 = est.Nos[e.n1]; const auto& n2 = est.Nos[e.n2];
                    double x1_def = n1->x + state.udesl(n1->gdlGlobais[0]); double y1_def = n1->y + state.udesl(n1->gdlGlobais[1]);
                    double x2_def = n2->x + state.udesl(n2->gdlGlobais[0]); double y2_def = n2->y + state.udesl(n2->gdlGlobais[1]);
                    double dx = x2_def - x1_def; double dy = y2_def - y1_def; double comp = std::hypot(dx, dy);
                    double nx = -dy / comp; double ny = dx / comp;
                    EsforcosLocais esf = est.Elementos[i]->getEsforcosLocais(state.udesl);
                    double v1 = 0, v2 = 0; Color c = RED;
                    if (tipoDiagrama == 1) { v1 = v2 = esf.N * escalaDiagrama; c = GREEN; }
                    else if (tipoDiagrama == 2) { v1 = v2 = esf.V * escalaDiagrama; c = ORANGE; }
                    else if (tipoDiagrama == 3) { v1 = esf.M1 * escalaDiagrama; v2 = -esf.M2 * escalaDiagrama; c = RED; }
                    Vector2 p1 = WorldToScreen(x1_def, y1_def, (float)GetScreenWidth(), (float)GetScreenHeight());
                    Vector2 p2 = WorldToScreen(x2_def, y2_def, (float)GetScreenWidth(), (float)GetScreenHeight());
                    Vector2 d1 = WorldToScreen(x1_def + nx * v1, y1_def + ny * v1, (float)GetScreenWidth(), (float)GetScreenHeight());
                    Vector2 d2 = WorldToScreen(x2_def + nx * v2, y2_def + ny * v2, (float)GetScreenWidth(), (float)GetScreenHeight());
                    DrawLineEx(p1, d1, 2.0f, c); DrawLineEx(p2, d2, 2.0f, c); DrawLineEx(d1, d2, 2.0f, c);
                }
            }
        }

        rlImGuiBegin();
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(350, 480), ImGuiCond_Once);
        ImGui::Begin("Configurações");
        
        ImGui::Text("Modelo Estrutural:");
        ImGui::RadioButton("Pórtico de Williams", &modeloSelecionado, 0);
        ImGui::RadioButton("Pilar Bi-engastado", &modeloSelecionado, 1);
        ImGui::RadioButton("Pilar Engastado", &modeloSelecionado, 4);
        ImGui::RadioButton("Viga Engastada", &modeloSelecionado, 2);
        ImGui::RadioButton("Viga Bi-apoiada", &modeloSelecionado, 3);
        ImGui::RadioButton("Viga Carga Distr.", &modeloSelecionado, 5);
        
        ImGui::Separator();
        ImGui::Text("Tipo de Análise:");
        ImGui::RadioButton("Linear", &analiseSelecionada, 0);
        ImGui::RadioButton("Não Linear (Arc-Length)", &analiseSelecionada, 1);
        
        ImGui::Separator();
        ImGui::Checkbox("Aplicar Imperfeição Inicial", &usarImperfeicao);
        if (usarImperfeicao) {
            ImGui::SliderFloat("Amplitude", &amplitudeImperfeicao, 0.0f, 0.5f, "%.3f m");
            if (ImGui::BeginCombo("Modo de Flambagem", TextFormat("Modo %d (L=%.2f)", modoImperfeicao + 1, modosBuckling.empty() ? 0.0 : modosBuckling[modoImperfeicao].lambda))) {
                for (int i = 0; i < (int)modosBuckling.size(); i++) {
                    bool is_selected = (modoImperfeicao == i);
                    if (ImGui::Selectable(TextFormat("Modo %d (L=%.2f)", i + 1, modosBuckling[i].lambda), is_selected))
                        modoImperfeicao = i;
                    if (is_selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        }

        ImGui::Separator();
        if (ImGui::Button("Calcular Nova Simulação", ImVec2(-1, 40))) ExecutarSimulacao();
        
        ImGui::Separator();
        ImGui::SliderInt("Passo de Carga", &current_step, 0, max_steps);
        ImGui::Text("Diagramas:");
        ImGui::RadioButton("Nenhum", &tipoDiagrama, 0); ImGui::SameLine();
        ImGui::RadioButton("N", &tipoDiagrama, 1); ImGui::SameLine();
        ImGui::RadioButton("V", &tipoDiagrama, 2); ImGui::SameLine();
        ImGui::RadioButton("M", &tipoDiagrama, 3);
        if (tipoDiagrama != 0) ImGui::SliderFloat("Escala Diag.", &escalaDiagrama, 0.00001f, 0.01f, "%.5f");
        
        ImGui::Separator();
        ImGui::Checkbox("Mostrar Reações de Apoio", &mostrarReacoes);
        if (mostrarReacoes) ImGui::SliderFloat("Escala Reações", &escalaReacoes, 0.00001f, 0.01f, "%.5f");

        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2((float)GetScreenWidth() - 520, 10), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiCond_Once);
        ImGui::Begin("Gráfico de Resposta");
        if (ImPlot::BeginPlot("Curva de Equilíbrio", ImVec2(-1, -1))) {
            ImPlot::SetupAxes("Deslocamento (v/h)", "Fator de Carga (Lambda)");
            ImPlot::PlotLine("Caminho", plot_u.data(), plot_f.data(), (int)plot_u.size());
            if (!historyUI.empty()) {
                double cu = historyUI[current_step].u_apex; double cf = historyUI[current_step].f_apex;
                ImPlot::PlotScatter("Atual", &cu, &cf, 1);
            }
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
