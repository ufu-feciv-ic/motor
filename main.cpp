#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>

#include "eigenpch.h"
#include <Eigen/Eigenvalues>

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
    virtual Eigen::MatrixXd getMatrizRigidezGeometrica(double N) const = 0;
    virtual void atualizarGeometria() = 0;
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
        atualizarGeometria();
    }

    void atualizarGeometria() override {
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

    Eigen::MatrixXd getMatrizRigidezGeometrica(double N) const override
    {
        double L = Linicial;
        Eigen::MatrixXd kgLocal = Eigen::MatrixXd::Zero(6, 6);
        double f = N / (30.0 * L);

        kgLocal << 
            0,  0,      0,      0,  0,      0,
            0,  36,     3*L,    0, -36,     3*L,
            0,  3*L,    4*L*L,  0, -3*L,   -L*L,
            0,  0,      0,      0,  0,      0,
            0, -36,    -3*L,    0,  36,    -3*L,
            0,  3*L,   -L*L,    0, -3*L,    4*L*L;
        
        kgLocal *= f;

        Eigen::MatrixXd T = calcularTransformacao();
        return T.transpose() * kgLocal * T;
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
        atualizarGeometria();
    }

    void atualizarGeometria() override {
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

    Eigen::MatrixXd getMatrizRigidezGeometrica(double N) const override
    {
        double L = L0;
        Eigen::MatrixXd kgLocal = Eigen::MatrixXd::Zero(6, 6);
        double f = N / (30.0 * L);

        kgLocal << 
            0,  0,      0,      0,  0,      0,
            0,  36,     3*L,    0, -36,     3*L,
            0,  3*L,    4*L*L,  0, -3*L,   -L*L,
            0,  0,      0,      0,  0,      0,
            0, -36,    -3*L,    0,  36,    -3*L,
            0,  3*L,   -L*L,    0, -3*L,    4*L*L;
        
        kgLocal *= f;

        // Transformação para global (usando geometria inicial)
        double dx = n2->x - n1->x;
        double dy = n2->y - n1->y;
        double C = dx / L;
        double S = dy / L;

        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);
        T << 
            C,  S, 0, 0, 0, 0,
           -S,  C, 0, 0, 0, 0,
            0,  0, 1, 0, 0, 0,
            0,  0, 0, C, S, 0,
            0,  0, 0,-S, C, 0,
            0,  0, 0, 0, 0, 1;

        return T.transpose() * kgLocal * T;
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

    void perturbarGeometria(const Eigen::VectorXd& modo, double amplitude)
    {
        double maxDesl = 0.0;
        for (const auto& no : Nos) {
            double dx = modo(no->gdlGlobais[0]);
            double dy = modo(no->gdlGlobais[1]);
            maxDesl = std::max(maxDesl, std::sqrt(dx*dx + dy*dy));
        }
        if (maxDesl < 1e-12) return;

        for (auto& no : Nos) {
            no->x += amplitude * (modo(no->gdlGlobais[0]) / maxDesl);
            no->y += amplitude * (modo(no->gdlGlobais[1]) / maxDesl);
        }
        
        for (auto& elem : Elementos) {
            elem->atualizarGeometria();
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

    static Eigen::MatrixXd montarMatrizRigidezGeometricaGlobal (const Estrutura& est, const Eigen::VectorXd& esforcosNormais)
    {
        Eigen::MatrixXd KGGlobalEst = Eigen::MatrixXd::Zero(est.NumGDLs, est.NumGDLs);

        for (size_t i = 0; i < est.Elementos.size(); ++i)
        {
            const auto& elemento = est.Elementos[i];
            Eigen::MatrixXd KGElem = elemento->getMatrizRigidezGeometrica(esforcosNormais(i));
            std::vector<int> GDLs = elemento->getGDLsGlobais();

            for (size_t row = 0; row < GDLs.size(); ++row)
                for (size_t col = 0; col < GDLs.size(); ++col)
                    KGGlobalEst(GDLs[row], GDLs[col]) += KGElem(row, col);
        }

        return KGGlobalEst;
    }
};

class AnaliseBuckling
{
public:
    struct ModoBuckling {
        Eigen::VectorXd modo;
        double lambda;
    };

    static std::vector<ModoBuckling> executar(Estrutura& est)
    {
        std::cout << "Iniciando Analise de Flambagem Linear...\n";

        // 1. Analise Linear de Referencia para obter esforcos normais
        Eigen::VectorXd uZero = Eigen::VectorXd::Zero(est.NumGDLs);
        Eigen::MatrixXd K = Construtor::montarMatrizRigidezGlobal(est, uZero);
        Eigen::VectorXd F = est.ForcasExternas;
        
        Eigen::MatrixXd K_solver = K;
        Eigen::VectorXd F_solver = F;
        est.aplicarCondicoesContorno(K_solver, F_solver);
        Eigen::VectorXd uRef = K_solver.colPivHouseholderQr().solve(F_solver);

        // 2. Calcular esforcos normais N_ref
        Eigen::VectorXd Nref = Eigen::VectorXd::Zero(est.Elementos.size());
        for (size_t i = 0; i < est.Elementos.size(); ++i) {
            EsforcosLocais esf = est.Elementos[i]->getEsforcosLocais(uRef);
            Nref(i) = esf.N;
        }

        // 3. Montar KG
        Eigen::MatrixXd KG = Construtor::montarMatrizRigidezGeometricaGlobal(est, Nref);

        // 4. Identificar GDLs livres
        std::vector<int> livres;
        std::vector<bool> eFixo(est.NumGDLs, false);
        for (int dof : est.NosFixos) eFixo[dof] = true;
        for (int i = 0; i < est.NumGDLs; ++i) if (!eFixo[i]) livres.push_back(i);

        int nf = livres.size();
        Eigen::MatrixXd Kff(nf, nf), KGff(nf, nf);
        for (int i = 0; i < nf; ++i) {
            for (int j = 0; j < nf; ++j) {
                Kff(i, j) = K(livres[i], livres[j]);
                KGff(i, j) = KG(livres[i], livres[j]);
            }
        }

        // 5. Resolver problema de autovalor: Kff * phi = -lambda * KGff * phi
        // Equivale a: (-Kff^-1 * KGff) * phi = (1/lambda) * phi
        Eigen::MatrixXd C = Kff.ldlt().solve(-KGff);
        Eigen::EigenSolver<Eigen::MatrixXd> es(C);

        std::vector<ModoBuckling> modos;
        for (int i = 0; i < es.eigenvalues().size(); ++i) {
            double mu = es.eigenvalues()(i).real();
            if (std::abs(es.eigenvalues()(i).imag()) < 1e-10 && mu > 1e-10) {
                double lambda = 1.0 / mu;
                Eigen::VectorXd modoF = es.eigenvectors().col(i).real();
                Eigen::VectorXd modoG = Eigen::VectorXd::Zero(est.NumGDLs);
                for (int j = 0; j < nf; ++j) modoG(livres[j]) = modoF(j);
                modos.push_back({modoG, lambda});
            }
        }

        // Ordenar por lambda crescente
        std::sort(modos.begin(), modos.end(), [](const ModoBuckling& a, const ModoBuckling& b) {
            return a.lambda < b.lambda;
        });

        return modos;
    }
};

struct Resultado
{
    Eigen::VectorXd u;
    Eigen::VectorXd F;
    Eigen::VectorXd reacoes; 
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

        Eigen::VectorXd uZero = Eigen::VectorXd::Zero(est.NumGDLs);
        Eigen::MatrixXd KGlobalOriginal = Construtor::montarMatrizRigidezGlobal(est, uZero);
        Eigen::MatrixXd KGlobalEst = KGlobalOriginal;
        Eigen::VectorXd FGlobal = est.ForcasExternas;

        est.aplicarCondicoesContorno(KGlobalEst, FGlobal);
        Eigen::VectorXd uGlobal = KGlobalEst.colPivHouseholderQr().solve(FGlobal);
        Eigen::VectorXd R = KGlobalOriginal * uGlobal - est.ForcasExternas;

        Resultado resultado;
        resultado.u = uGlobal;
        resultado.reacoes = R;
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
        historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

        for (int passo = 1; passo <= numPassos; ++passo)
        {
            double lambda = (double)passo / numPassos;
            Eigen::VectorXd PassoCarga = lambda * est.ForcasExternas;

            int iter = 0;
            double erro = 1.0;

            while (iter < maxIter)
            {
                Eigen::MatrixXd KTangente = Construtor::montarMatrizRigidezGlobal(est, uAtual);
                Eigen::VectorXd Fint = Construtor::montarForcasInternasGlobais(est, uAtual);
                Eigen::VectorXd g = PassoCarga - Fint;

                for (int dof : est.NosFixos) g(dof) = 0.0;
                erro = g.norm();
                if (erro < tol) break;

                est.aplicarCondicoesContorno(KTangente, g);
                Eigen::VectorXd deltaU = KTangente.ldlt().solve(g);
                uAtual += deltaU;
                iter++;
            }

            Eigen::VectorXd FintFinal = Construtor::montarForcasInternasGlobais(est, uAtual);
            Eigen::VectorXd R = FintFinal - PassoCarga;

            historico.push_back(Resultado{uAtual, PassoCarga, R, lambda});
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

        historico.push_back({uAtual, Eigen::VectorXd::Zero(est.NumGDLs), Eigen::VectorXd::Zero(est.NumGDLs), 0.0});

        for (int passo = 1; passo <= numPassos; ++passo)
        {
            Eigen::MatrixXd Kt = Construtor::montarMatrizRigidezGlobal(est, uAtual);
            Eigen::VectorXd dummy = Eigen::VectorXd::Zero(est.NumGDLs);
            est.aplicarCondicoesContorno(Kt, dummy);
            Eigen::VectorXd deltaUr = Kt.colPivHouseholderQr().solve(est.ForcasExternas);

            double normDur = deltaUr.norm();
            if (normDur < 1e-12) normDur = 1e-12;
            double DeltaLambda = deltal / normDur;
            if (passo > 1 && DELTAU.dot(deltaUr) < 0.0) DeltaLambda = -DeltaLambda;

            Eigen::VectorXd deltaU = DeltaLambda * deltaUr;
            Eigen::VectorXd DELTAU0 = deltaU;
            DELTAU = deltaU;
            double lambdaTentativo = lambda + DeltaLambda;

            int iter = 0;
            while (iter < maxIter)
            {
                Eigen::VectorXd Fint = Construtor::montarForcasInternasGlobais(est, uAtual + DELTAU);
                Eigen::VectorXd g = lambdaTentativo * est.ForcasExternas - Fint;
                for (int dof : est.NosFixos) g(dof) = 0.0;
                if (g.norm() <= tol) break;

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

            uAtual += DELTAU;
            lambda = lambdaTentativo;
            if (iter > 0) deltal = deltal0 * std::sqrt(Nd / static_cast<double>(iter));

            Eigen::VectorXd FintFinal = Construtor::montarForcasInternasGlobais(est, uAtual);
            Eigen::VectorXd R = FintFinal - (lambda * est.ForcasExternas);
            historico.push_back({uAtual, lambda * est.ForcasExternas, R, lambda});
        }
        return historico;
    }
};

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
