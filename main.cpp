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

double normalizaAngulo (double angulo)
{
    double a = std::fmod(angulo + M_PI, 2.0 * M_PI);
    if (a < 0.0) a += 2 * M_PI;
    return a - M_PI;
};

struct PropriedadesMaterial
{
    double E;
    double A;
    double I;
};

class No
{
public:
    int id;
    double x, y;
    std::vector<int> gdlGlobais;

    No (int id, double x, double y, std::vector<int> gdls)
    : id(id), x(x), y(y), gdlGlobais(gdls) {}
};

struct EsforcosLocais
{
    double N; // Força normal
    double V; // força cortante (constante para cargar nodais)
    double M1; // Momento fletor no nó 1
    double M2; // Momento fletor no nó 2
};

class ElementoFinito
{
public:
    virtual ~ElementoFinito() = default;
    virtual std::vector<int> getGDLsGlobais() const = 0;
    virtual Eigen::MatrixXd getMatrizRigidezGlobal(const Eigen::VectorXd& uGlobal) const = 0;
    virtual Eigen::VectorXd getForcasInternasGlobais(const Eigen::VectorXd& uGlobal) const = 0;
    virtual EsforcosLocais getEsforcosLocais(const Eigen::VectorXd& uGlobal) const = 0;
};

class Viga2DLinear : public ElementoFinito
{
private: 
    std::shared_ptr<No> n1, n2;
    PropriedadesMaterial material;
    double Linicial, cosInicial, senInicial;

public:
    Viga2DLinear(std::shared_ptr<No> n1, std::shared_ptr<No> n2, PropriedadesMaterial material)
    : n1(n1), n2(n2), material(material) {

        double dx = n2->x - n1->x;
        double dy = n2->y - n1->y;
        Linicial = std::sqrt(dx * dx + dy * dy);
        cosInicial = dx / Linicial;
        senInicial = dy / Linicial;
    }

    std::vector<int> getGDLsGlobais() const override
    {
        std::vector<int> GDLs;
        GDLs.insert(GDLs.end(), n1->gdlGlobais.begin(), n1->gdlGlobais.end());
        GDLs.insert(GDLs.end(), n2->gdlGlobais.begin(), n2->gdlGlobais.end());
        return GDLs;
    }

    Eigen::MatrixXd calcularkLocal() const
    {
        // Parâmetros 
        double EAL = (material.E * material.A) / Linicial;
        double EIL = (material.E * material.I) / Linicial;
        double EIL2 = (material.E * material.I) / (Linicial * Linicial);
        double EIL3 = (material.E * material.I) / (Linicial * Linicial * Linicial);

        // Matriz de rigidez local da barra
        Eigen::MatrixXd kLocal = Eigen::MatrixXd::Zero(6, 6);
        kLocal << 
            EAL, 0, 0, -EAL, 0, 0,
            0, 12*EIL3, 6*EIL2, 0, -12*EIL3, 6*EIL2, 
            0, 6*EIL2, 4*EIL, 0, -6*EIL2, 2*EIL,
            -EAL, 0, 0, EAL, 0, 0,
            0, -12*EIL3, -6*EIL2, 0, 12*EIL3, -6*EIL2,
            0, 6*EIL2, 2*EIL, 0, -6*EIL2, 4*EIL;

        return kLocal;
    }

    Eigen::MatrixXd calcularTransformacao() const
    {
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);
        T << 
            cosInicial, senInicial, 0, 0, 0, 0,
            -senInicial, cosInicial, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, cosInicial, senInicial, 0,
            0, 0, 0, -senInicial, cosInicial, 0,
            0, 0, 0, 0, 0, 1;

        return T;
    }

    Eigen::MatrixXd getMatrizRigidezGlobal (const Eigen::VectorXd& uGlobal) const override
    {
        // Matriz de rigidez local da barra
        Eigen::MatrixXd kLocal = calcularkLocal();

        // Matriz de transformação das coordenadas locais para globais
        Eigen::MatrixXd T = calcularTransformacao();

        // Cálculo da matriz de rigidez em coordenadas globais
        return T.transpose() * kLocal * T;    
    }

    Eigen::VectorXd getForcasInternasGlobais (const Eigen::VectorXd& uGlobal) const override
    {
        std::vector<int> GLDs = getGDLsGlobais();
        Eigen::VectorXd uGlobalElem = Eigen::VectorXd::Zero(6);

        for (int i = 0; i < 6; ++i)
        {
            uGlobalElem(i) = uGlobal(GLDs[i]);
        }

        // Matriz de rigidez local da barra
        Eigen::MatrixXd kLocal = calcularkLocal();

        // Matriz de transformação das coordenadas locais para globais
        Eigen::MatrixXd T = calcularTransformacao();

        // Global para local
        Eigen::VectorXd uLocal = T * uGlobalElem;
        Eigen::VectorXd fLocal = kLocal * uLocal;
        
        // Forças internas para global
        return T.transpose() * fLocal;
    }

    EsforcosLocais getEsforcosLocais(const Eigen::VectorXd& uGlobal) const override 
    {
        std::vector<int> GLDs = getGDLsGlobais();
        Eigen::VectorXd uGlobalElem = Eigen::VectorXd::Zero(6);
        for (int i = 0; i < 6; ++i) uGlobalElem(i) = uGlobal(GLDs[i]);

        Eigen::MatrixXd kLocal = calcularkLocal();
        Eigen::MatrixXd T = calcularTransformacao();

        Eigen::VectorXd uLocal = T * uGlobalElem;
        Eigen::VectorXd fLocal = kLocal * uLocal;

        // Na matriz local: fLocal = [N1, V1, M1, N2, V2, M2]
        // N é constante (assumimos tração positiva pegando de N2), M1 = fLocal(2), M2 = fLocal(5)
        double N = fLocal(3); 
        double V = fLocal(1); // ou -fLocal(4) dependendo da convenção
        double M1 = fLocal(2);
        double M2 = fLocal(5);

        return {N, V, M1, M2};
    }
};

class Viga2DCorrotacional : public ElementoFinito
{
private: 
    std::shared_ptr<No> n1, n2;
    PropriedadesMaterial mat;
    double L0;

public:
    Viga2DCorrotacional(std::shared_ptr<No> no1, std::shared_ptr<No> no2, PropriedadesMaterial m)
        : n1(no1), n2(no2), mat(m) {
        
        // Calcula o comprimento original (indeslocado) no momento da criação
        double dx = n2->x - n1->x;
        double dy = n2->y - n1->y;
        L0 = std::sqrt(dx * dx + dy * dy);
    }

    std::vector<int> getGDLsGlobais() const override
    {
        std::vector<int> GDLs;
        GDLs.insert(GDLs.end(), n1->gdlGlobais.begin(), n1->gdlGlobais.end());
        GDLs.insert(GDLs.end(), n2->gdlGlobais.begin(), n2->gdlGlobais.end());
        return GDLs;
    }

    // Calcula a matriz de rigidez global do elemento
    Eigen::MatrixXd getMatrizRigidezGlobal(const Eigen::VectorXd& uGlobal) const override
    {
        std::vector<int> GDLs = getGDLsGlobais();
        Eigen::VectorXd uGlobalElem = Eigen::VectorXd::Zero(6);

        for (int i = 0; i < 6; ++i) uGlobalElem(i) = uGlobal(GDLs[i]);

        double X1 = n1->x; double Y1 = n1->y;
        double X2 = n2->x; double Y2 = n2->y;

        // Geometria atualizada
        double dx = (X2 + uGlobalElem(3)) - (X1 + uGlobalElem(0));
        double dy = (Y2 + uGlobalElem(4)) - (Y1 + uGlobalElem(1));
        double L = std::sqrt(dx * dx + dy * dy);
        if (L < 1e-12) L = 1e-12;

        double C = dx / L;
        double S = dy / L;
        double beta0 = std::atan2(Y2 - Y1, X2 - X1);
        double beta = std::atan2(dy, dx);

        double teta1 = normalizaAngulo(uGlobalElem(2) + beta0 - beta);
        double teta2 = normalizaAngulo(uGlobalElem(5) + beta0 - beta);

        // Deformação de Green-Lagrange
        double ul = (L * L - L0 * L0) / (L + L0);

        // Esforços internos locais
        double N = mat.E * mat.A * ul / L0;
        double constanteFlexao = 2.0 * mat.E * mat.I / L0;
        double M1 = constanteFlexao * (2.0 * teta1 + teta2);
        double M2 = constanteFlexao * (teta1 + 2.0 * teta2);

        // Matriz B de transformação (Trabalha direto do sistema básico 3DOF pro Global 6DOF)
        Eigen::Matrix<double, 3, 6> B;
        B << -C,   -S,   0.0,  C,    S,   0.0,
             -S/L,  C/L, 1.0,  S/L, -C/L, 0.0,
             -S/L,  C/L, 0.0,  S/L, -C/L, 1.0;

        Eigen::Matrix3d D;
        D << mat.E*mat.A/L0, 0, 0,
             0, 4.0*mat.E*mat.I/L0, 2.0*mat.E*mat.I/L0,
             0, 2.0*mat.E*mat.I/L0, 4.0*mat.E*mat.I/L0;
        
        // Matriz tangente material 
        Eigen::Matrix<double, 6, 6> KM = B.transpose() * D * B;

        Eigen::Vector<double, 6> rVec, zVec;
        rVec << -C, -S, 0, C, S, 0;
        zVec << S, -C, 0, -S, C, 0;

        // Matrizes Tangentes Geométricas (K_1 e K_2) - Dependem da força atual!
        Eigen::Matrix<double, 6, 6> K1 = (N / L) * (zVec * zVec.transpose());
        Eigen::Matrix<double, 6, 6> K2 = ((M1 + M2) / (L * L)) * (rVec * zVec.transpose() + zVec * rVec.transpose());

        return KM + K1 + K2;
    }

    // Calcula o Vetor de Forças Internas Global (F_int)
    Eigen::VectorXd getForcasInternasGlobais(const Eigen::VectorXd& uGlobal) const override {
        std::vector<int> GDLs = getGDLsGlobais();
        Eigen::Vector<double, 6> uGlobalElem;
        for (int i = 0; i < 6; ++i) uGlobalElem(i) = uGlobal(GDLs[i]);

        double X1 = n1->x; double Y1 = n1->y;
        double X2 = n2->x; double Y2 = n2->y;

        double dx = (X2 + uGlobalElem(3)) - (X1 + uGlobalElem(0));
        double dy = (Y2 + uGlobalElem(4)) - (Y1 + uGlobalElem(1));
        double L = std::sqrt(dx * dx + dy * dy);
        if (L < 1e-12) L = 1e-12;

        double C = dx / L; 
        double S = dy / L;
        double beta0 = std::atan2(Y2 - Y1, X2 - X1);
        double beta = std::atan2(dy, dx);

        double teta1 = normalizaAngulo(uGlobalElem(2) + beta0 - beta);
        double teta2 = normalizaAngulo(uGlobalElem(5) + beta0 - beta);
        double ul = (L * L - L0 * L0) / (L + L0);

        // Esforços internos locais
        double N = mat.E * mat.A * ul / L0;
        double constanteFlexao = 2.0 * mat.E * mat.I / L0;
        double M1 = constanteFlexao * (2.0 * teta1 + teta2);
        double M2 = constanteFlexao * (teta1 + 2.0 * teta2);

        Eigen::Vector3d fLocal(N, M1, M2);

        Eigen::Matrix<double, 3, 6> B;
        B << -C,   -S,   0.0,  C,    S,   0.0,
             -S/L,  C/L, 1.0,  S/L, -C/L, 0.0,
             -S/L,  C/L, 0.0,  S/L, -C/L, 1.0;

        // O B^T já transforma a força de volta para os eixos X e Y globais reais
        Eigen::Vector<double, 6> FGlobal = B.transpose() * fLocal;
        return FGlobal;
    }

    EsforcosLocais getEsforcosLocais(const Eigen::VectorXd& uGlobal) const override 
    {
        std::vector<int> GDLs = getGDLsGlobais();
        Eigen::Vector<double, 6> uGlobalElem;
        for (int i = 0; i < 6; ++i) uGlobalElem(i) = uGlobal(GDLs[i]);

        double X1 = n1->x; double Y1 = n1->y;
        double X2 = n2->x; double Y2 = n2->y;

        // Comprimento atualizado
        double dx = (X2 + uGlobalElem(3)) - (X1 + uGlobalElem(0));
        double dy = (Y2 + uGlobalElem(4)) - (Y1 + uGlobalElem(1));
        double L = std::sqrt(dx * dx + dy * dy);
        if (L < 1e-12) L = 1e-12;

        double beta0 = std::atan2(Y2 - Y1, X2 - X1);
        double beta = std::atan2(dy, dx);

        double teta1 = normalizaAngulo(uGlobalElem(2) + beta0 - beta);
        double teta2 = normalizaAngulo(uGlobalElem(5) + beta0 - beta);
        double ul = (L * L - L0 * L0) / (L + L0);

        // Esforços internos locais
        double N = mat.E * mat.A * ul / L0;
        double constanteFlexao = 2.0 * mat.E * mat.I / L0;
        double M1 = constanteFlexao * (2.0 * teta1 + teta2);
        double M2 = constanteFlexao * (teta1 + 2.0 * teta2);
        
        // Cortante obtido por equilíbrio no estado deformado (somatório de momentos = 0)
        double V = (M1 + M2) / L; 

        return {N, V, M1, M2};
    }
};

class Estrutura
{
public:
    int NumGDLs = 0;
    std::vector<std::shared_ptr<No>> Nos;
    std::vector<std::shared_ptr<ElementoFinito>> Elementos;
    Eigen::VectorXd ForcasExternas;
    std::vector<int> NosFixos;

    // Adicioanar um nó e atualizar a quantidade total de graus de liberdade
    void adicionarNo (std::shared_ptr<No> no)
    {
        Nos.push_back(no);
        NumGDLs += no->gdlGlobais.size();
    }

    void adicionarElemento (std::shared_ptr<ElementoFinito> elemento)
    {
        Elementos.push_back(elemento);
    }

    void aplicarCondicoesContorno(Eigen::MatrixXd& K, Eigen::VectorXd& F) const
    {
        for (int GDL : NosFixos)
        {
            K.row(GDL).setZero();
            K.col(GDL).setZero();
            K(GDL, GDL) = 1.0;
            F(GDL) = 0.0;
        }
    }
};

class Construtor
{
public: 
    static Eigen::MatrixXd montarMatrizRigidezGlobal (const Estrutura& est, const Eigen::VectorXd& uGlobal)
    {
        Eigen::MatrixXd KGlobalEst = Eigen::MatrixXd::Zero(est.NumGDLs, est.NumGDLs);

        for (const auto& elemento : est.Elementos)
        {
            // Pega a matriz do elemento já em coordenadas globais
            Eigen::MatrixXd KGlobalElem = elemento->getMatrizRigidezGlobal(uGlobal);

            // Pega os graus de liberdade globais 
            std::vector<int> GDLs = elemento->getGDLsGlobais();

            for (size_t i = 0; i < GDLs.size(); ++i)
                for (size_t j = 0; j < GDLs.size(); ++j)
                    KGlobalEst(GDLs[i], GDLs[j]) += KGlobalElem(i, j);
        }

        return KGlobalEst;
    }

    static Eigen::VectorXd montarForcasInternasGlobais (const Estrutura& est, const Eigen::VectorXd& uGlobal)
    {
        Eigen::VectorXd FGlobal = Eigen::VectorXd::Zero(est.NumGDLs);

        for (const auto& elemento : est.Elementos)
        {
            Eigen::VectorXd FGlobalElem = elemento->getForcasInternasGlobais(uGlobal);
            std::vector<int> GDLs = elemento->getGDLsGlobais();

            for (size_t i = 0; i < GDLs.size(); ++i)
                FGlobal(GDLs[i]) += FGlobalElem(i);
        }

        return FGlobal;
    }
};

struct Resultado
{
    Eigen::VectorXd u;
    Eigen::VectorXd F;
    double FatorCarga;
};

class EstrategiaAnalise
{
public:
    virtual ~EstrategiaAnalise() = default;
    virtual std::vector<Resultado> executar(Estrutura& est) = 0;
};

class AnaliseLinear : public EstrategiaAnalise
{
public:
    std::vector<Resultado> executar(Estrutura& est) override
    {
        std::cout << "Iniciando solucao do sistema linear...\n";

        // Construir a matriz global da estrutura
        Eigen::VectorXd uZero = Eigen::VectorXd::Zero(est.NumGDLs);
        Eigen::MatrixXd KGlobalEst = Construtor::montarMatrizRigidezGlobal(est, uZero);

        // Vetor de forças externas
        Eigen::VectorXd FGlobal = est.ForcasExternas;
        std::cout << "Vetor de Forcas Externas (FGlobal):\n" << FGlobal << "\n";

        // Aplicar condições de contorno
        est.aplicarCondicoesContorno(KGlobalEst, FGlobal);
        std::cout << "Vetor após aplicar condições de contorno (FGlobal):\n" << FGlobal << "\n";

        // Resolver o sistema
        Eigen::LLT<Eigen::MatrixXd> llt(KGlobalEst);
        Eigen::VectorXd uGlobal;

        if (llt.info() == Eigen::Success)
        {
            uGlobal = llt.solve(FGlobal);
        }
        else
        {
            std::cout << "-> Aviso: LLT falhou. Tentando HouseholderQR...\n";
            uGlobal = KGlobalEst.colPivHouseholderQr().solve(FGlobal);
        }

        Resultado resultado;
        resultado.u = uGlobal;
        resultado.FatorCarga = 1.0;

        return {resultado};
    }
};

class AnaliseNaoLinearNR : public EstrategiaAnalise
{
private:
    int numPassos;
    int maxIter;
    double tol;

public:
    AnaliseNaoLinearNR(int passos = 10, int iteracoes = 50, double tolerancia = 1e-6)
    : numPassos(passos), maxIter(iteracoes), tol(tolerancia) {}

    std::vector<Resultado> executar(Estrutura& est) override
    {
        std::cout << "--- Iniciando Solver Newton-Raphson ---\n";
        std::vector<Resultado> historico;

        Eigen::VectorXd uAtual = Eigen::VectorXd::Zero(est.NumGDLs);
        historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

        // Loop incremental
        for (int passo = 1; passo <= numPassos; ++passo)
        {
            // Cálculo do fator de carga
            double lambda = (double)passo / numPassos;
            Eigen::VectorXd PassoCarga = lambda * est.ForcasExternas;

            int iter = 0;
            double erro = 1.0;

            // Loop iterativo
            while (iter < maxIter)
            {
                // Coletar matriz tangente e forças internas
                Eigen::MatrixXd KTangente = Construtor::montarMatrizRigidezGlobal(est, uAtual);
                Eigen::MatrixXd Fint = Construtor::montarForcasInternasGlobais(est, uAtual);

                // Calcular vetor de forças residuais
                Eigen::VectorXd g = PassoCarga - Fint;

                // Zerar o resíduo nos graus de liberdade fixos
                for (int dof : est.NosFixos)
                    g(dof) = 0.0;
                
                
                // Critérios de parada
                erro = g.norm();

                if (erro < tol) break;

                // Aplicar condições de contorno
                est.aplicarCondicoesContorno(KTangente, g);

                // Resolver o sistema 
                Eigen::VectorXd deltaU = KTangente.ldlt().solve(g);

                // Atualizar deslocamentos 
                uAtual += deltaU;
                iter++;
            }

            if (iter > maxIter) 
            {
                std::cout << "AVISO: Passo " << passo << " nao convergiu!\n";
            }
            else
            {
                std::cout << "Passo " << passo << " (Lambda=" << lambda
                          << ") convergiu em " << iter << " iteracoes. Erro final: " << erro << "\n";
            }

            historico.push_back(Resultado{uAtual, PassoCarga, lambda});
        }

        return historico;
    }
};

class AnaliseNaoLinearCompArco : public EstrategiaAnalise
{
private:
    int numPassos;
    int maxIter;
    double tol;
    double deltal0;
    double Nd;

public:
    AnaliseNaoLinearCompArco(
        int passos = 30,
        int iteracoes = 50,
        double tolerancia = 1e-6,
        double comprimentoArcoInicial = 0.025,
        double iterDesejadas = 5.0)
        : numPassos(passos),
          maxIter(iteracoes),
          tol(tolerancia),
          deltal0(comprimentoArcoInicial),
          Nd(iterDesejadas) {}

    std::vector<Resultado> executar(Estrutura& est) override
    {
        std::cout << "--- Iniciando Solver Arc-Length ---\n";

        std::vector<Resultado> historico;

        Eigen::VectorXd uAtual = Eigen::VectorXd::Zero(est.NumGDLs);
        Eigen::VectorXd DELTAU = Eigen::VectorXd::Zero(est.NumGDLs);

        double lambda = 0.0;
        double deltal = deltal0;

        historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

        for (int passo = 1; passo <= numPassos; ++passo)
        {
            // =========================
            // PREDITOR
            // =========================
            Eigen::MatrixXd Kt = Construtor::montarMatrizRigidezGlobal(est, uAtual);
            Eigen::VectorXd dummy = Eigen::VectorXd::Zero(est.NumGDLs);
            est.aplicarCondicoesContorno(Kt, dummy);

            Eigen::VectorXd deltaUr = Kt.colPivHouseholderQr().solve(est.ForcasExternas);

            double normDur = deltaUr.norm();
            if (normDur < 1e-12) normDur = 1e-12;

            double DeltaLambda = deltal / normDur;

            if (passo > 1)
            {
                if (DELTAU.dot(deltaUr) < 0.0)
                    DeltaLambda = -DeltaLambda;
            }

            Eigen::VectorXd deltaU = DeltaLambda * deltaUr;
            Eigen::VectorXd DELTAU0 = deltaU;
            DELTAU = deltaU;

            double lambdaTentativo = lambda + DeltaLambda;

            // =========================
            // CORRETOR
            // =========================
            int iter = 0;
            double erro = 1.0;

            while (iter < maxIter)
            {
                Eigen::VectorXd Fint =
                    Construtor::montarForcasInternasGlobais(est, uAtual + DELTAU);

                Eigen::VectorXd g = lambdaTentativo * est.ForcasExternas - Fint;

                for (int dof : est.NosFixos)
                    g(dof) = 0.0;

                erro = g.norm();
                if (erro <= tol)
                    break;

                Kt = Construtor::montarMatrizRigidezGlobal(est, uAtual + DELTAU);
                est.aplicarCondicoesContorno(Kt, dummy);

                Eigen::VectorXd deltaUg = Kt.colPivHouseholderQr().solve(g);
                deltaUr = Kt.colPivHouseholderQr().solve(est.ForcasExternas);

                double topo = DELTAU0.dot(deltaUg);
                double base = DELTAU0.dot(deltaUr);

                double dlambda = (std::abs(base) < 1e-12) ? 0.0 : -topo / base;

                DELTAU += deltaUg + dlambda * deltaUr;
                DeltaLambda += dlambda;
                lambdaTentativo = lambda + DeltaLambda;

                iter++;
            }

            if (iter >= maxIter)
            {
                std::cout << "AVISO: passo " << passo << " nao convergiu. Erro final = "
                          << erro << "\n";
            }
            else
            {
                std::cout << "Passo " << passo
                          << " convergiu em " << iter
                          << " iteracoes. Lambda = " << lambdaTentativo
                          << " | Erro = " << erro << "\n";
            }

            // Aceita o passo
            uAtual += DELTAU;
            lambda = lambdaTentativo;

            if (iter > 0)
                deltal = deltal0 * std::sqrt(Nd / static_cast<double>(iter));

            historico.push_back({uAtual, lambda * est.ForcasExternas, lambda});
        }

        return historico;
    }
};

struct ResultadoPassoUI
{
    Eigen::VectorXd udesl;   // deslocamentos globais completos
    double lambda;           // fator de carga
    double u_apex;           // deslocamento normalizado v/h
    double f_apex;           // força equivalente no gráfico
};

struct ArestaRender
{
    int n1;
    int n2;
};

std::vector<ResultadoPassoUI> prepararHistoricoUI(
    const std::vector<Resultado>& history,
    int dofApexY,
    double h_apex)
{
    std::vector<ResultadoPassoUI> out;
    out.reserve(history.size());

    for (const auto& step : history)
    {
        ResultadoPassoUI s;
        s.udesl = step.u;
        s.lambda = step.FatorCarga;

        double u_y_apex = step.u(dofApexY);

        // mesmo critério que você já usa
        s.u_apex = -u_y_apex / h_apex;
        s.f_apex = step.FatorCarga; // ou step.FatorCarga
        // Se quiser o mesmo sinal do seu print:
        // s.f_apex = -step.FatorCarga * (-1.0);

        out.push_back(s);
    }

    return out;
}

Vector2 WorldToScreen(double x, double y, float screenW, float screenH) {
    float scale = 35.0f; 
    float offsetX = 100.0f;
    float offsetY = screenH * 0.7f; 
    return { static_cast<float>(offsetX + x * scale), static_cast<float>(offsetY - y * scale) };
}

int main()
{
    std::cout << "--- TESTE ETAPA FINAL: PORTICO DE WILLIAMS (ARC-LENGTH) ---\n\n";

    Estrutura est;
    std::vector<ArestaRender> arestas;

    // Propriedades normalizadas clássicas do problema de Williams
    PropriedadesMaterial mat = {1.0, 1.885e6, 9.274e3};

    // Coordenadas dos 11 nós
    std::vector<std::pair<double, double>> coords = {
        {0.0, 0.0}, {2.5872, 0.0736}, {5.1744, 0.1472}, {7.7616, 0.2208},
        {10.3488, 0.2944}, {12.936, 0.368}, {15.5232, 0.2944}, {18.1104, 0.2208},
        {20.6976, 0.1472}, {23.2848, 0.0736}, {25.872, 0.0}
    };

    // Nós
    int dof_count = 0;
    for (int i = 0; i < 11; ++i) {
        auto no = std::make_shared<No>(
            i,
            coords[i].first,
            coords[i].second,
            std::vector<int>{dof_count, dof_count + 1, dof_count + 2}
        );
        est.adicionarNo(no);
        dof_count += 3;
    }

    // Elementos
    for (int i = 0; i < 10; ++i) {
        auto viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i + 1], mat);
        est.adicionarElemento(viga);
        arestas.push_back({i, i + 1});
    }

    // Apoios
    est.NosFixos = {0, 1, 2, 30, 31, 32};

    // Carga externa
    est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    est.ForcasExternas(16) = -1.0; // DOF Y do nó 5

    // Solver
    std::cout << "Calculando simulacao...\n";
    AnaliseNaoLinearCompArco solver(33, 50, 1e-6, 0.025, 5.0);
    std::vector<Resultado> history = solver.executar(est);
    std::cout << "Simulacao concluida. Abrindo UI...\n";

    const Resultado& ultimoPasso = history.back();
    EsforcosLocais esf = est.Elementos[0]->getEsforcosLocais(ultimoPasso.u);

    std::cout << "\nElemento 0 - Normal: " << esf.N 
          << ", Cortante: " << esf.V 
          << ", M1: " << esf.M1 
          << ", M2: " << esf.M2 << "\n";

    // Preparar dados para a UI
    double h_apex = 0.368;
    int dofApexY = 16;
    std::vector<ResultadoPassoUI> historyUI = prepararHistoricoUI(history, dofApexY, h_apex);

    std::vector<double> plot_u, plot_f;
    plot_u.reserve(historyUI.size());
    plot_f.reserve(historyUI.size());

    for (const auto& step : historyUI) {
        plot_u.push_back(step.u_apex);
        plot_f.push_back(step.f_apex);
    }

    // Setup Raylib
    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(1280, 720, "Williams Frame - Raylib + ImGui + ImPlot");
    SetTargetFPS(60);

    // Setup ImGui / ImPlot
    rlImGuiSetup(true);
    ImPlot::CreateContext();

    int current_step = 0;
    int max_steps = static_cast<int>(historyUI.size()) - 1;

    int tipoDiagrama = 0; // 0-Nenhum, 1-Normal, 2-Cortante, 3-Momento
    float escalaDiagrama = 0.002f;

    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(Color{40, 40, 45, 255});

        const ResultadoPassoUI& state = historyUI[current_step];

        // =========================================================
        // DESENHO DA ESTRUTURA
        // =========================================================

        // Estrutura original
        for (const auto& e : arestas)
        {
            const auto& no1 = est.Nos[e.n1];
            const auto& no2 = est.Nos[e.n2];

            Vector2 p1 = WorldToScreen(no1->x, no1->y, (float)GetScreenWidth(), (float)GetScreenHeight());
            Vector2 p2 = WorldToScreen(no2->x, no2->y, (float)GetScreenWidth(), (float)GetScreenHeight());

            DrawLineEx(p1, p2, 1.5f, Color{120, 120, 120, 120});
        }

        // Estrutura deformada
        for (const auto& e : arestas)
        {
            const auto& no1 = est.Nos[e.n1];
            const auto& no2 = est.Nos[e.n2];

            int gdl1x = no1->gdlGlobais[0];
            int gdl1y = no1->gdlGlobais[1];
            int gdl2x = no2->gdlGlobais[0];
            int gdl2y = no2->gdlGlobais[1];

            Vector2 p1 = WorldToScreen(
                no1->x + state.udesl(gdl1x),
                no1->y + state.udesl(gdl1y),
                (float)GetScreenWidth(),
                (float)GetScreenHeight()
            );

            Vector2 p2 = WorldToScreen(
                no2->x + state.udesl(gdl2x),
                no2->y + state.udesl(gdl2y),
                (float)GetScreenWidth(),
                (float)GetScreenHeight()
            );

            DrawLineEx(p1, p2, 3.5f, SKYBLUE);
            DrawCircleV(p1, 5.0f, WHITE);
            DrawCircleV(p2, 5.0f, WHITE);
        }

        if (tipoDiagrama != 0)
        {
            for (size_t i = 0; i < arestas.size(); ++i)
            {
                const auto& e = arestas[i];
                const auto& no1 = est.Nos[e.n1];
                const auto& no2 = est.Nos[e.n2];

                int gdl1x = no1->gdlGlobais[0]; int gdl1y = no1->gdlGlobais[1];
                int gdl2x = no2->gdlGlobais[0]; int gdl2y = no2->gdlGlobais[1];

                double x1_def = no1->x + state.udesl(gdl1x);
                double y1_def = no1->y + state.udesl(gdl1y);
                double x2_def = no2->x + state.udesl(gdl2x);
                double y2_def = no2->y + state.udesl(gdl2y);

                // Vetor normal (perpendicular) à barra
                double dx = x2_def - x1_def;
                double dy = y2_def - y1_def;
                double comp = std::hypot(dx, dy);
                double nx = -dy / comp; 
                double ny = dx / comp;  

                EsforcosLocais esf = est.Elementos[i]->getEsforcosLocais(state.udesl);

                // Variáveis para armazenar os valores do diagrama nos nós 1 e 2
                double val1 = 0.0, val2 = 0.0;
                Color corBorda = BLACK;
                Color corPreench = BLANK;

                // Define os valores e cores baseado no diagrama selecionado
                if (tipoDiagrama == 1) { // Normal (Constante)
                    val1 = esf.N * escalaDiagrama;
                    val2 = esf.N * escalaDiagrama;
                    corBorda = GREEN;
                    corPreench = ColorAlpha(GREEN, 0.4f);
                } 
                else if (tipoDiagrama == 2) { // Cortante (Constante)
                    val1 = esf.V * escalaDiagrama;
                    val2 = esf.V * escalaDiagrama;
                    corBorda = ORANGE;
                    corPreench = ColorAlpha(ORANGE, 0.4f);
                } 
                else if (tipoDiagrama == 3) { // Momento (Variavel)
                    val1 = esf.M1 * escalaDiagrama;
                    val2 = -esf.M2 * escalaDiagrama; 
                    corBorda = RED;
                    corPreench = ColorAlpha(RED, 0.4f);
                }

                // Pontos na base da barra (tela)
                Vector2 p1 = WorldToScreen(x1_def, y1_def, (float)GetScreenWidth(), (float)GetScreenHeight());
                Vector2 p2 = WorldToScreen(x2_def, y2_def, (float)GetScreenWidth(), (float)GetScreenHeight());

                // Pontos do topo do diagrama (tela)
                Vector2 d1 = WorldToScreen(x1_def + nx * val1, y1_def + ny * val1, (float)GetScreenWidth(), (float)GetScreenHeight());
                Vector2 d2 = WorldToScreen(x2_def + nx * val2, y2_def + ny * val2, (float)GetScreenWidth(), (float)GetScreenHeight());

                // Desenha o preenchimento e contorno do diagrama
                auto DrawTriangleRobust = [](Vector2 v1, Vector2 v2, Vector2 v3, Color col) {
                    DrawTriangle(v1, v2, v3, col);
                    DrawTriangle(v1, v3, v2, col); // Desenha com ordem invertida para garantir visibilidade (culling)
                };

                if ((val1 > 0 && val2 < 0) || (val1 < 0 && val2 > 0)) {
                    // Inversão de sinal: calcula o ponto onde o diagrama cruza a barra
                    float t = (float)(std::abs(val1) / (std::abs(val1) + std::abs(val2)));
                    Vector2 p_zero = { p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y) };

                    // Triângulo 1 (nó 1 até o zero)
                    DrawTriangleRobust(p1, d1, p_zero, corPreench);
                    DrawLineEx(p1, d1, 2.0f, corBorda);
                    DrawLineEx(d1, p_zero, 2.0f, corBorda);

                    // Triângulo 2 (zero até o nó 2)
                    DrawTriangleRobust(p_zero, d2, p2, corPreench);
                    DrawLineEx(p_zero, d2, 2.0f, corBorda);
                    DrawLineEx(d2, p2, 2.0f, corBorda);
                } else {
                    // Mesmo sinal ou um dos valores é zero: desenha como um trapézio (2 triângulos)
                    DrawTriangleRobust(p1, d1, d2, corPreench);
                    DrawTriangleRobust(p1, d2, p2, corPreench);

                    // Contorno
                    DrawLineEx(p1, d1, 2.0f, corBorda);
                    DrawLineEx(d1, d2, 2.0f, corBorda);
                    DrawLineEx(d2, p2, 2.0f, corBorda);
                    // Linha da barra (opcional, já que a barra principal é desenhada)
                }
            }
        }

        // =========================================================
        // UI
        // =========================================================
        rlImGuiBegin();

        ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(350, 170), ImGuiCond_Once);
        ImGui::Begin("Painel de Controle");

        ImGui::Text("Portico de Williams - Analise Nao Linear");
        ImGui::Separator();

        ImGui::SliderInt("Passo de Carga", &current_step, 0, max_steps);

        if (ImGui::Button("Anterior") && current_step > 0) current_step--;
        ImGui::SameLine();
        if (ImGui::Button("Proximo") && current_step < max_steps) current_step++;

        ImGui::Separator();
        
        // NOVO: Controles Unificados de Diagramas
        ImGui::Text("Diagramas de Esforcos:");
        ImGui::RadioButton("Nenhum", &tipoDiagrama, 0); ImGui::SameLine();
        ImGui::RadioButton("Normal (N)", &tipoDiagrama, 1); ImGui::SameLine();
        ImGui::RadioButton("Cortante (V)", &tipoDiagrama, 2); 
        ImGui::RadioButton("Momento Fletor (M)", &tipoDiagrama, 3);
        
        if (tipoDiagrama != 0) {
            ImGui::SliderFloat("Escala Visual", &escalaDiagrama, 0.0001f, 0.01f, "%.5f");
        }
        
        ImGui::Separator();

        ImGui::Separator();
        ImGui::Text("Fator de Carga (Lambda): %.6f", state.lambda);
        ImGui::Text("Desl. Normalizado (v/h): %.6f", state.u_apex);
        ImGui::Text("Forca P (normalizada): %.6f", state.f_apex);

        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2((float)GetScreenWidth() - 620.0f, 20.0f), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(600, 400), ImGuiCond_Once);
        ImGui::Begin("Grafico de Resposta P x v/h");

        if (ImPlot::BeginPlot("Curva de Equilibrio", ImVec2(-1, -1)))
        {
            ImPlot::SetupAxes("Deslocamento Normalizado (v/h)", "Forca P");

            // opcional: ajustar automaticamente
            // ImPlot::SetupAxesLimits(0, 3.0, 0, 1.5, ImGuiCond_Once);

            ImPlot::SetNextLineStyle(ImVec4(0.5f, 0.5f, 0.5f, 0.5f), 1.0f);
            ImPlot::PlotLine("Caminho Total", plot_u.data(), plot_f.data(), (int)plot_u.size());

            ImPlot::SetNextLineStyle(ImVec4(1.0f, 0.2f, 0.2f, 1.0f), 2.5f);
            ImPlot::PlotLine("Resposta Atual", plot_u.data(), plot_f.data(), current_step + 1);

            double curr_u = state.u_apex;
            double curr_f = state.f_apex;

            ImPlot::SetNextMarkerStyle(
                ImPlotMarker_Circle,
                5.0f,
                ImVec4(1, 1, 0, 1),
                1.0f,
                ImVec4(1, 1, 0, 1)
            );
            ImPlot::PlotScatter("Ponto de Equilibrio", &curr_u, &curr_f, 1);

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

    // std::cout << "--- TESTE ETAPA FINAL: PORTICO DE WILLIAMS (ARC-LENGTH) ---\n\n";

    // Estrutura est;
    
    // // Propriedades normalizadas clássicas do problema de Williams
    // PropriedadesMaterial mat = {1.0, 1.885e6, 9.274e3};

    // // Coordenadas dos 11 nós (formando o arco suave)
    // std::vector<std::pair<double, double>> coords = {
    //     {0.0, 0.0}, {2.5872, 0.0736}, {5.1744, 0.1472}, {7.7616, 0.2208},
    //     {10.3488, 0.2944}, {12.936, 0.368}, {15.5232, 0.2944}, {18.1104, 0.2208},
    //     {20.6976, 0.1472}, {23.2848, 0.0736}, {25.872, 0.0}
    // };

    // // 1. Criando os Nós Dinamicamente (Total de 33 graus de liberdade)
    // int dof_count = 0;
    // for (int i = 0; i < 11; ++i) {
    //     auto no = std::make_shared<No>(i, coords[i].first, coords[i].second, 
    //                                      std::vector<int>{dof_count, dof_count+1, dof_count+2});
    //     est.adicionarNo(no);
    //     dof_count += 3;
    // }

    // std::vector<ArestaRender> arestas;

    // // 2. Conectando 10 Vigas Corrotacionais
    // for (int i = 0; i < 10; ++i) {
    //     auto viga = std::make_shared<Viga2DCorrotacional>(est.Nos[i], est.Nos[i+1], mat);
    //     est.adicionarElemento(viga);
    // }

    // // 3. Condições de Contorno: Engaste total nas duas extremidades
    // // Nó 0 (dofs 0, 1, 2) e Nó 10 (dofs 30, 31, 32)
    // est.NosFixos = {0, 1, 2, 30, 31, 32};

    // // 4. Força Externa de Referência (Apertando o nó central para baixo)
    // // O Nó 5 é o ápice do arco. O DOF Y dele é o 16.
    // est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    // est.ForcasExternas(16) = -1.0;

    // // 5. Roda o Arc-Length
    // AnaliseNaoLinearCompArco solver(33, 50, 1e-6, 0.025, 5.0); // 30 passos na curva
    // std::vector<Resultado> history = solver.executar(est);

    // // 6. Impressão dos Resultados Normalizados
    // double h_apex = 0.368;
    
    // std::cout << "\n=== TRACADO DA CURVA DE FLAMBAGEM (SNAP-THROUGH) ===\n";
    // std::cout << " Passo | Desloc (v/h) | Forca Lambda \n";
    // std::cout << "--------------------------------------\n";
    
    // for (size_t i = 1; i < history.size(); ++i) {
    //     double lambda = history[i].FatorCarga;
    //     double u_y_apex = history[i].u(16); 
        
    //     // Normalização clássica do paper de Williams
    //     double u_normalizado = -u_y_apex / h_apex;
    //     double f_normalizada = -lambda * (-1.0); 

    //     std::cout << std::scientific << i 
    //               << " | " << std::scientific << u_normalizado 
    //               << " | " << std::scientific << f_normalizada 
    //               << "\n";
    // }

    // return 0;

    // std::cout << "--- TESTE ETAPA 4.1: GRANDES DEFORMACOES (CORROTACIONAL + NR) ---\n\n";

    // Estrutura est;
    // PropriedadesMaterial mat = {210E9, 0.01, 1e-5};

    // // Nós: Comprimento total de 10 metros, dividido em 2 elementos de 5m
    // auto no1 = std::make_shared<No>(1, 0.0, 0.0, std::vector<int>{0, 1, 2});
    // auto no2 = std::make_shared<No>(2, 5.0, 0.0, std::vector<int>{3, 4, 5});
    // auto no3 = std::make_shared<No>(3, 10.0, 0.0, std::vector<int>{6, 7, 8});

    // est.adicionarNo(no1);
    // est.adicionarNo(no2);
    // est.adicionarNo(no3);

    // // Elementos
    // auto barra1 = std::make_shared<Viga2DCorrotacional>(no1, no2, mat);
    // auto barra2 = std::make_shared<Viga2DCorrotacional>(no2, no3, mat);

    // est.adicionarElemento(barra1);
    // est.adicionarElemento(barra2);

    // // Engaste nó 1
    // est.NosFixos = {0, 1, 2};

    // // Carga alta no no3
    // est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    // est.ForcasExternas(7) = 1000.0;

    // AnaliseNaoLinearNR analiseEstrutural{10, 50, 1e-3};
    // std::vector<Resultado> historico = analiseEstrutural.executar(est);

    // Eigen::VectorXd uFinal = historico.back().u;

    // std::cout << "\n=== RESULTADO FINAL (Lambda 1.0) ===\n";
    // std::cout << "Deslocamento Y na ponta (Deve ser um valor positivo grande): " << uFinal(7) << " m\n";
    
    // std::cout << "\nO Pulo do Gato Não Linear:\n";
    // std::cout << "Deslocamento X na ponta (Teoria linear diz que eh ZERO): " << uFinal(6) << " m\n";

    // if (uFinal(6) < -0.1) {
    //     std::cout << "\n-> SUCESSO! O deslocamento X eh negativo. A barra recuou para compensar a curvatura geométrica (Grandes deformacoes comprovadas!).\n";
    // }

    // return 0;

    // std::cout << "--- TESTE ETAPA 3: SOLVER E CONDICOES DE CONTORNO ---\n\n";

    // Estrutura est;
    // PropriedadesMaterial mat = {210E9, 0.01, 0.0001};

    // // Nó 1: Origem (0, 0)
    // auto no1 = std::make_shared<No>(1, 0.0, 0.0, std::vector<int>{0, 1, 2});
    // // Nó 2: Na ponta, L = 5 metros
    // auto no2 = std::make_shared<No>(2, 5.0, 0.0, std::vector<int>{3, 4, 5});

    // est.adicionarNo(no1);
    // est.adicionarNo(no2);

    // auto barra = std::make_shared<Viga2DLinear>(no1, no2, mat);
    // est.adicionarElemento(barra);

    // // Engaste no nó 1
    // est.NosFixos = {0, 1, 2};

    // // Carga P = -1000 kN vertical no nó 2
    // est.ForcasExternas = Eigen::VectorXd::Zero(est.NumGDLs);
    // est.ForcasExternas(4) = -1000.0;
    
    // // Resolução
    // AnaliseNaoLinear analiseEstrutural;
    // std::vector<Resultado> historico = analiseEstrutural.executar(est);
    // Eigen::VectorXd uFinal = historico.back().u;

    // // Resultados 
    // std::cout << "\nDeslocamentos finais (Vetor u):\n";
    // std::cout << uFinal << "\n";

    // // Validação Teórica
    // double P = -1000.0;
    // double L = 5.0;
    // double E = 210E9;
    // double I = 0.0001;

    // double flechaTeorica = (P * std::pow(L, 3)) / (3.0 * E * I);
    // double rotacaoTeorica = (P * std::pow(L, 2)) / (2.0 * E * I);

    // std::cout << "--- VALIDACAO ANALITICA ---\n";
    // std::cout << "Deslocamento Y no No 2 (Algoritmo): " << std::scientific << uFinal(4) << " m\n";
    // std::cout << "Deslocamento Y no No 2 (Teorico)  : " << std::scientific << flechaTeorica << " m\n";
    // std::cout << "\nRotacao Z no No 2 (Algoritmo)     : " << std::scientific << uFinal(5) << " rad\n";
    // std::cout << "Rotacao Z no No 2 (Teorico)       : " << std::scientific << rotacaoTeorica << " rad\n";

    // return 0;

    // std::cout << "--- TESTE ETAPA 2: MONTADOR E MODELO DE ESTRUTURA ---\n\n";

    // Estrutura est;
    // PropriedadesMaterial mat = {210E9, 0.01, 0.0001};

    // // 1. Criando 3 Nós (Pórtico em "L" deitado)
    // // Nó 1 na origem (0, 0)
    // auto no1 = std::make_shared<No>(1, 0.0, 0.0, std::vector<int>{0, 1, 2});
    // // Nó 2 em X=3 (Viga 1 é horizontal, comprimento 3)
    // auto no2 = std::make_shared<No>(2, 3.0, 0.0, std::vector<int>{3, 4, 5});
    // // Nó 3 em X=3, Y=4 (Viga 2 é vertical, comprimento 4)
    // auto no3 = std::make_shared<No>(3, 3.0, 4.0, std::vector<int>{6, 7, 8});

    // est.adicionarNo(no1);
    // est.adicionarNo(no2);
    // est.adicionarNo(no3);

    // // 2. Criando 2 barras
    // auto barra1 = std::make_shared<Viga2DLinear>(no1, no2, mat);
    // auto barra2 = std::make_shared<Viga2DLinear>(no2, no3, mat);

    // est.adicionarElemento(barra1);
    // est.adicionarElemento(barra2);

    // Eigen::VectorXd uGlobal = Eigen::VectorXd::Zero(est.NumGDLs);
    // Eigen::MatrixXd KGlobalEst = Construtor::montarMatrizRigidezGlobal(est, uGlobal);

    // std::cout << "Matriz de rigidez global da Estrutura:\n";
    // std::cout << KGlobalEst << std::endl;

    // return 0;

    // std::cout << "--- TESTE ETAPA 1: ELEMENTO LINEAR ---\n\n";

    // // 1. Nó 1 na origem (0,0), Nó 2 em (3, 4). Comprimento = 5.
    // auto no1 = std::make_shared<Node>(1, 0.0, 0.0, std::vector<int>{0, 1, 2});
    // auto no2 = std::make_shared<Node>(2, 3.0, 4.0, std::vector<int>{3, 4, 5});

    // // 2. Material
    // PropriedadesMaterial material = {210E9, 0.01, 0.0001};

    // // 3. Instaciar a barra
    // std::shared_ptr<Element> barra = std::make_shared<Viga2DLinear>(no1, no2, material);

    // // 4. Testar BCN
    // std::cout << "Graus de liberdade (BCN) do elemento: ";
    // for (int gdl : barra->getGDLsGlobais())
    // {
    //     std::cout << gdl << " ";
    // }
    // std::cout << "\n\n";

    // // 5. Matriz de rigidez global, vetor de deslocamentos zerado
    // Eigen::VectorXd uGlobal = Eigen::VectorXd::Zero(6);
    // Eigen::MatrixXd KGlobal = barra.getMatrizRigidezGlobal(uGlobal);

    // std::cout << "Matriz de rigidez global do elemento:\n" << KGlobal << "\n";
    // std::cout << KGlobal << std::endl;
}


// #define _USE_MATH_DEFINES
// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <iomanip>

// // Raylib, ImGui e ImPlot
// #include <raylib.h>
// #include "imgui.h"
// #include "implot.h"
// #include "rlImGui.h" 

// // Eigen
// // #include <Eigen/Dense>
// #include "eigenpch.h"

// // --- ESTRUTURAS DE DADOS ---

// struct ModelData {
//     int NTEL, NTNOS, NTGL, NNOSCC;
//     Eigen::MatrixXd coord;
//     Eigen::MatrixXi inci;
//     Eigen::MatrixXi dofno;
//     Eigen::VectorXd E, A, I, Fr;
//     std::vector<int> NOCC;
// };

// struct StepResult {
//     Eigen::VectorXd udesl;     
//     double lambda;      
//     double u_apex;      
//     double f_apex;      
// };

// enum class CalcType { Force, Stiffness };

// // --- SUB-ROTINAS MATEMÁTICAS ---

// double normalize_angle(double angle) {
//     double a = std::fmod(angle + M_PI, 2.0 * M_PI);
//     if (a < 0.0) a += 2.0 * M_PI;
//     return a - M_PI;
// }

// ModelData krenk(double P_val) {
//     ModelData md;
//     md.NTNOS = 11;
//     md.NTEL = 10;
//     md.NTGL = md.NTNOS * 3;
//     md.NNOSCC = 6;

//     md.coord = Eigen::MatrixXd(md.NTNOS, 2);
//     md.coord << 0.0, 0.0,
//                 2.5872, 0.0736,
//                 5.1744, 0.1472,
//                 7.7616, 0.2208,
//                 10.3488, 0.2944,
//                 12.936, 0.368,   
//                 15.5232, 0.2944,
//                 18.1104, 0.2208,
//                 20.6976, 0.1472,
//                 23.2848, 0.0736,
//                 25.872, 0.0;

//     md.inci = Eigen::MatrixXi(md.NTEL, 2);
//     md.inci << 0, 1,   1, 2,   2, 3,   3, 4,   4, 5,
//                5, 6,   6, 7,   7, 8,   8, 9,   9, 10;

//     md.E = Eigen::VectorXd::Constant(md.NTEL, 1.0);
//     md.A = Eigen::VectorXd::Constant(md.NTEL, 1.885e6);
//     md.I = Eigen::VectorXd::Constant(md.NTEL, 9.274e3);

//     md.dofno = Eigen::MatrixXi(md.NTEL, 6);
//     for (int i = 0; i < md.NTEL; ++i) {
//         int no1 = md.inci(i, 0);
//         int no2 = md.inci(i, 1);
//         md.dofno(i, 0) = no1 * 3;     md.dofno(i, 1) = no1 * 3 + 1; md.dofno(i, 2) = no1 * 3 + 2;
//         md.dofno(i, 3) = no2 * 3;     md.dofno(i, 4) = no2 * 3 + 1; md.dofno(i, 5) = no2 * 3 + 2;
//     }

//     md.Fr = Eigen::VectorXd::Zero(md.NTGL);
//     md.Fr(16) = P_val; 

//     md.NOCC = {0, 1, 2, 30, 31, 32}; 

//     return md;
// }

// void element_routine(const Eigen::VectorXd& udesl, const Eigen::MatrixXd& coord0, const ModelData& md, 
//                      CalcType calc_type, Eigen::VectorXd& out_force, Eigen::MatrixXd& out_stiffness) {
    
//     if (calc_type == CalcType::Force) {
//         out_force = Eigen::VectorXd::Zero(md.NTGL);
//     } else {
//         out_stiffness = Eigen::MatrixXd::Zero(md.NTGL, md.NTGL);
//     }

//     for (int m = 0; m < md.NTEL; ++m) {
//         int no1 = md.inci(m, 0);
//         int no2 = md.inci(m, 1);
//         double X1 = coord0(no1, 0); double Y1 = coord0(no1, 1);
//         double X2 = coord0(no2, 0); double Y2 = coord0(no2, 1);

//         Eigen::Vector<double, 6> u_el;
//         for (int i = 0; i < 6; ++i) u_el(i) = udesl(md.dofno(m, i));

//         double L0 = std::sqrt(std::pow(X2 - X1, 2) + std::pow(Y2 - Y1, 2));

//         double dx = (X2 + u_el(3)) - (X1 + u_el(0));
//         double dy = (Y2 + u_el(4)) - (Y1 + u_el(1));
//         double L = std::sqrt(dx * dx + dy * dy);
//         if (L < 1e-12) L = 1e-12;

//         double C = dx / L;
//         double S = dy / L;

//         double beta0 = std::atan2(Y2 - Y1, X2 - X1);
//         double beta = std::atan2(dy, dx);

//         double theta1_g = u_el(2);
//         double theta2_g = u_el(5);

//         double teta1 = normalize_angle(theta1_g + beta0 - beta);
//         double teta2 = normalize_angle(theta2_g + beta0 - beta);

//         double ul = (L * L - L0 * L0) / (L + L0);

//         double Em = md.E(m); double Am = md.A(m); double Im = md.I(m);

//         double N = Em * Am * ul / L0;
//         double k_bend = 2.0 * Em * Im / L0;
//         double M1 = k_bend * (2.0 * teta1 + teta2);
//         double M2 = k_bend * (teta1 + 2.0 * teta2);

//         Eigen::Matrix<double, 3, 6> B;
//         B.row(0) << -C, -S, 0, C, S, 0;
//         B.row(1) << -S/L, C/L, 1.0, S/L, -C/L, 0.0;
//         B.row(2) << -S/L, C/L, 0.0, S/L, -C/L, 1.0;

//         if (calc_type == CalcType::Force) {
//             Eigen::Vector3d f_local(N, M1, M2);
//             Eigen::Vector<double, 6> Fel = B.transpose() * f_local;
//             for (int i = 0; i < 6; ++i) {
//                 out_force(md.dofno(m, i)) += Fel(i);
//             }
//         } else {
//             Eigen::Matrix3d D;
//             D << Em*Am/L0, 0, 0,
//                  0, 4.0*Em*Im/L0, 2.0*Em*Im/L0,
//                  0, 2.0*Em*Im/L0, 4.0*Em*Im/L0;

//             Eigen::Matrix<double, 6, 6> KM = B.transpose() * D * B;

//             Eigen::Vector<double, 6> r_vec, z_vec;
//             r_vec << -C, -S, 0, C, S, 0;
//             z_vec << S, -C, 0, -S, C, 0;

//             Eigen::Matrix<double, 6, 6> K1 = (N / L) * (z_vec * z_vec.transpose());
//             Eigen::Matrix<double, 6, 6> K2 = ((M1 + M2) / (L * L)) * (r_vec * z_vec.transpose() + z_vec * r_vec.transpose());

//             Eigen::Matrix<double, 6, 6> Kel = KM + K1 + K2;

//             for (int i = 0; i < 6; ++i) {
//                 for (int j = 0; j < 6; ++j) {
//                     out_stiffness(md.dofno(m, i), md.dofno(m, j)) += Kel(i, j);
//                 }
//             }
//         }
//     }
// }

// void apply_bc(Eigen::MatrixXd* K, Eigen::VectorXd* R, const std::vector<int>& NOCC) {
//     for (int dof : NOCC) {
//         if (K != nullptr) {
//             K->row(dof).setZero();
//             K->col(dof).setZero();
//             (*K)(dof, dof) = 1.0;
//         }
//         if (R != nullptr) {
//             (*R)(dof) = 0.0;
//         }
//     }
// }

// // --- SIMULAÇÃO ---

// std::vector<StepResult> run_simulation(ModelData& md, double P_inc, int nmax) {
//     std::vector<StepResult> history;
    
//     double tol = 1e-6; int kmax = 50;
//     double deltal0 = 0.025; double Nd = 5.0;
    
//     Eigen::VectorXd udesl = Eigen::VectorXd::Zero(md.NTGL);
//     Eigen::VectorXd DELTAU = Eigen::VectorXd::Zero(md.NTGL);
//     double lambda_val = 0.0;
//     double deltal = deltal0;
//     Eigen::MatrixXd coord0 = md.coord;
    
//     double h_apex = 0.368;
    
//     // Passo 0
//     history.push_back({udesl, 0.0, 0.0, 0.0});

//     Eigen::VectorXd out_force;
//     Eigen::MatrixXd K;

//     for (int np_step = 1; np_step <= nmax; ++np_step) {
        
//         // --- Preditor ---
//         element_routine(udesl, coord0, md, CalcType::Stiffness, out_force, K);
//         apply_bc(&K, nullptr, md.NOCC);

//         Eigen::VectorXd deltaur = K.colPivHouseholderQr().solve(md.Fr);

//         double norm_dur = deltaur.norm();
//         if (norm_dur == 0) norm_dur = 1e-10;
//         double Dlambda = deltal / norm_dur;

//         if (np_step > 1) {
//             if (DELTAU.dot(deltaur) < 0) {
//                 Dlambda = -Dlambda;
//             }
//         }

//         Eigen::VectorXd deltau = Dlambda * deltaur;
//         Eigen::VectorXd DELTAU0 = deltau;
//         DELTAU = deltau;

//         double lambda_temp = lambda_val + Dlambda;

//         // --- Corretor ---
//         element_routine(udesl + DELTAU, coord0, md, CalcType::Force, out_force, K);
//         Eigen::VectorXd g = (lambda_temp * md.Fr) - out_force;
//         apply_bc(nullptr, &g, md.NOCC);

//         int k = 0;
//         while (k < kmax) {
//             k++;

//             element_routine(udesl + DELTAU, coord0, md, CalcType::Stiffness, out_force, K);
//             apply_bc(&K, nullptr, md.NOCC);

//             Eigen::VectorXd deltaug = K.colPivHouseholderQr().solve(g);
//             deltaur = K.colPivHouseholderQr().solve(md.Fr);

//             double top = DELTAU0.dot(deltaug);
//             double bot = DELTAU0.dot(deltaur);

//             double dlambda = (std::abs(bot) < 1e-12) ? 0.0 : -top / bot;

//             DELTAU += (deltaug + dlambda * deltaur);
//             Dlambda += dlambda;
//             lambda_temp = lambda_val + Dlambda;

//             element_routine(udesl + DELTAU, coord0, md, CalcType::Force, out_force, K);
//             g = (lambda_temp * md.Fr) - out_force;
//             apply_bc(nullptr, &g, md.NOCC);

//             if (g.norm() <= tol) {
//                 break;
//             }
//         }

//         udesl += DELTAU;
//         lambda_val = lambda_temp;

//         if (k > 0) deltal = deltal0 * std::sqrt(Nd / k);

//         double u_y_apex = udesl(16); 
        
//         StepResult res;
//         res.udesl = udesl;
//         res.lambda = lambda_val;
//         res.u_apex = -u_y_apex / h_apex;
//         res.f_apex = -lambda_val * P_inc;
        
//         history.push_back(res);

//         // Imprime o passo, o deslocamento normalizado (v/h) e a força P no console
//         std::cout << "Passo " << std::setw(2) << np_step 
//                   << " | Desloc (v/h): " << std::fixed << std::setprecision(4) << std::setw(8) << res.u_apex 
//                   << " | Força P: "      << std::fixed << std::setprecision(4) << std::setw(8) << res.f_apex 
//                   << "\n";
//     }
//     return history;
// }

// // --- INTERFACE E RENDERIZAÇÃO ---

// Vector2 WorldToScreen(double x, double y, float screenW, float screenH) {
//     float scale = 35.0f; 
//     float offsetX = 100.0f;
//     float offsetY = screenH * 0.7f; 
//     return { static_cast<float>(offsetX + x * scale), static_cast<float>(offsetY - y * scale) };
// }

// // ==============================================================================
// // MAIN
// // ==============================================================================

// int main() {
//     std::cout << "Calculando simulação (NLG)...\n";
//     ModelData md = krenk(-1.0); 
//     std::vector<StepResult> history = run_simulation(md, -1.0, 33);
//     std::cout << "Simulação concluída. Abrindo UI...\n";

//     // Preparar dados para o ImPlot 
//     std::vector<double> plot_u, plot_f;
//     for(const auto& step : history) {
//         plot_u.push_back(step.u_apex);
//         plot_f.push_back(step.f_apex);
//     }

//     // Setup Raylib
//     SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
//     InitWindow(1280, 720, "Williams Frame - Raylib + ImGui + ImPlot");
//     SetTargetFPS(60);

//     // Setup ImGui & ImPlot
//     rlImGuiSetup(true);
//     ImPlot::CreateContext();

//     int current_step = 0;
//     int max_steps = history.size() - 1;

//     while (!WindowShouldClose()) {
//         BeginDrawing();
//         ClearBackground({ 40, 40, 45, 255 });

//         // --- DESENHO DA ESTRUTURA (RAYLIB) ---
//         StepResult& state = history[current_step];
//         for (int i = 0; i < md.NTEL; i++) {
//             int n1 = md.inci(i, 0); int n2 = md.inci(i, 1);
            
//             // Pega as posições deformadas
//             Vector2 p1 = WorldToScreen(md.coord(n1,0) + state.udesl(n1*3), md.coord(n1,1) + state.udesl(n1*3+1), GetScreenWidth(), GetScreenHeight());
//             Vector2 p2 = WorldToScreen(md.coord(n2,0) + state.udesl(n2*3), md.coord(n2,1) + state.udesl(n2*3+1), GetScreenWidth(), GetScreenHeight());
            
//             DrawLineEx(p1, p2, 3.5f, SKYBLUE);
//             DrawCircleV(p1, 5.0f, WHITE);
//             DrawCircleV(p2, 5.0f, WHITE);
//         }

//         // --- UI IMGUI & IMPLOT ---
//         rlImGuiBegin();

//         // Janela de Controle
//         ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiCond_Once);
//         ImGui::SetNextWindowSize(ImVec2(350, 150), ImGuiCond_Once);
//         ImGui::Begin("Painel de Controle");
//         ImGui::Text("Pórtico de Williams - Análise Não Linear");
//         ImGui::Separator();
        
//         ImGui::SliderInt("Passo de Carga", &current_step, 0, max_steps);
//         if (ImGui::Button("Anterior") && current_step > 0) current_step--;
//         ImGui::SameLine();
//         if (ImGui::Button("Próximo") && current_step < max_steps) current_step++;

//         ImGui::Separator();
//         ImGui::Text("Fator de Carga (Lambda): %.4f", state.lambda);
//         ImGui::Text("Desl. Normalizado (v/h): %.4f", state.u_apex);
//         ImGui::Text("Força P (lb): %.4f", state.f_apex);
//         ImGui::End();

//         // Janela do Gráfico com ImPlot
//         ImGui::SetNextWindowPos(ImVec2(GetScreenWidth() - 620.0f, 20.0f), ImGuiCond_Once);
//         ImGui::SetNextWindowSize(ImVec2(600, 400), ImGuiCond_Once);
//         ImGui::Begin("Gráfico de Resposta P x v/h");
        
//         if (ImPlot::BeginPlot("Curva de Equilíbrio", ImVec2(-1, -1))) {
//             ImPlot::SetupAxes("Deslocamento Normalizado (v/h)", "Força P (lb)");
//             ImPlot::SetupAxesLimits(0, 3.0, 0, 1.5);

//             // Desenha a curva completa em cinza transparente
//             ImPlot::SetNextLineStyle(ImVec4(0.5f, 0.5f, 0.5f, 0.5f), 1.0f);
//             ImPlot::PlotLine("Caminho Total", plot_u.data(), plot_f.data(), (int)plot_u.size());

//             // Desenha o caminho percorrido até o passo atual em vermelho
//             ImPlot::SetNextLineStyle(ImVec4(1.0f, 0.2f, 0.2f, 1.0f), 2.5f);
//             ImPlot::PlotLine("Resposta Atual", plot_u.data(), plot_f.data(), current_step + 1);

//             // Marcador no ponto atual
//             double curr_u = state.u_apex;
//             double curr_f = state.f_apex;
//             ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 5.0f, ImVec4(1, 1, 0, 1), 1.0f, ImVec4(1, 1, 0, 1));
//             ImPlot::PlotScatter("Ponto de Equilíbrio", &curr_u, &curr_f, 1);

//             ImPlot::EndPlot();
//         }
//         ImGui::End();

//         rlImGuiEnd();
//         EndDrawing();
//     }

//     // Limpeza
//     ImPlot::DestroyContext();
//     rlImGuiShutdown();
//     CloseWindow();

//     return 0;
// }