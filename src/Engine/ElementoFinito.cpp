#include "Engine/ElementoFinito.h"
#include <cmath>
#include <algorithm>

namespace Engine {

// --- Viga2DLinear ---

Viga2DLinear::Viga2DLinear(std::shared_ptr<No> n1, std::shared_ptr<No> n2, PropriedadesMaterial material)
    : n1(n1), n2(n2), material(material) {
    atualizarGeometria();
}

void Viga2DLinear::atualizarGeometria() {
    double dx = n2->x - n1->x;
    double dy = n2->y - n1->y;
    Linicial = std::sqrt(dx * dx + dy * dy);
    cosInicial = dx / Linicial;
    senInicial = dy / Linicial;
}

std::vector<int> Viga2DLinear::getGDLsGlobais() const {
    std::vector<int> GDLs;
    GDLs.insert(GDLs.end(), n1->gdlGlobais.begin(), n1->gdlGlobais.end());
    GDLs.insert(GDLs.end(), n2->gdlGlobais.begin(), n2->gdlGlobais.end());
    return GDLs;
}

Eigen::MatrixXd Viga2DLinear::calcularkLocal() const {
    double EAL = (material.E * material.A) / Linicial;
    double EIL = (material.E * material.I) / Linicial;
    double EIL2 = (material.E * material.I) / (Linicial * Linicial);
    double EIL3 = (material.E * material.I) / (Linicial * Linicial * Linicial);

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

Eigen::MatrixXd Viga2DLinear::calcularTransformacao() const {
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

Eigen::MatrixXd Viga2DLinear::getMatrizRigidezGlobal(const Eigen::VectorXd& uGlobal) const {
    Eigen::MatrixXd kLocal = calcularkLocal();
    Eigen::MatrixXd T = calcularTransformacao();
    return T.transpose() * kLocal * T;
}

Eigen::VectorXd Viga2DLinear::getForcasInternasGlobais(const Eigen::VectorXd& uGlobal) const {
    std::vector<int> GLDs = getGDLsGlobais();
    Eigen::VectorXd uGlobalElem = Eigen::VectorXd::Zero(6);
    for (int i = 0; i < 6; ++i) uGlobalElem(i) = uGlobal(GLDs[i]);

    Eigen::MatrixXd kLocal = calcularkLocal();
    Eigen::MatrixXd T = calcularTransformacao();

    Eigen::VectorXd uLocal = T * uGlobalElem;
    Eigen::VectorXd fLocal = kLocal * uLocal;
    return T.transpose() * fLocal;
}

EsforcosLocais Viga2DLinear::getEsforcosLocais(const Eigen::VectorXd& uGlobal) const {
    std::vector<int> GLDs = getGDLsGlobais();
    Eigen::VectorXd uGlobalElem = Eigen::VectorXd::Zero(6);
    for (int i = 0; i < 6; ++i) uGlobalElem(i) = uGlobal(GLDs[i]);

    Eigen::MatrixXd kLocal = calcularkLocal();
    Eigen::MatrixXd T = calcularTransformacao();

    Eigen::VectorXd uLocal = T * uGlobalElem;
    Eigen::VectorXd fLocal = kLocal * uLocal;

    double N = fLocal(3); 
    double V = fLocal(1);
    double M1 = fLocal(2);
    double M2 = fLocal(5);

    return {N, V, M1, M2};
}

Eigen::MatrixXd Viga2DLinear::getMatrizRigidezGeometrica(double N) const {
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

// --- Viga2DCorrotacional ---

Viga2DCorrotacional::Viga2DCorrotacional(std::shared_ptr<No> no1, std::shared_ptr<No> no2, PropriedadesMaterial m)
    : n1(no1), n2(no2), mat(m) {
    atualizarGeometria();
}

void Viga2DCorrotacional::atualizarGeometria() {
    double dx = n2->x - n1->x;
    double dy = n2->y - n1->y;
    L0 = std::sqrt(dx * dx + dy * dy);
}

std::vector<int> Viga2DCorrotacional::getGDLsGlobais() const {
    std::vector<int> GDLs;
    GDLs.insert(GDLs.end(), n1->gdlGlobais.begin(), n1->gdlGlobais.end());
    GDLs.insert(GDLs.end(), n2->gdlGlobais.begin(), n2->gdlGlobais.end());
    return GDLs;
}

Eigen::MatrixXd Viga2DCorrotacional::getMatrizRigidezGlobal(const Eigen::VectorXd& uGlobal) const {
    std::vector<int> GDLs = getGDLsGlobais();
    Eigen::VectorXd uGlobalElem = Eigen::VectorXd::Zero(6);
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

    double N = mat.E * mat.A * ul / L0;
    double constanteFlexao = 2.0 * mat.E * mat.I / L0;
    double M1 = constanteFlexao * (2.0 * teta1 + teta2);
    double M2 = constanteFlexao * (teta1 + 2.0 * teta2);

    Eigen::Matrix<double, 3, 6> B;
    B << -C,   -S,   0.0,  C,    S,   0.0,
         -S/L,  C/L, 1.0,  S/L, -C/L, 0.0,
         -S/L,  C/L, 0.0,  S/L, -C/L, 1.0;

    Eigen::Matrix3d D;
    D << mat.E*mat.A/L0, 0, 0,
         0, 4.0*mat.E*mat.I/L0, 2.0*mat.E*mat.I/L0,
         0, 2.0*mat.E*mat.I/L0, 4.0*mat.E*mat.I/L0;
    
    Eigen::Matrix<double, 6, 6> KM = B.transpose() * D * B;

    Eigen::Vector<double, 6> rVec, zVec;
    rVec << -C, -S, 0, C, S, 0;
    zVec << S, -C, 0, -S, C, 0;

    Eigen::Matrix<double, 6, 6> K1 = (N / L) * (zVec * zVec.transpose());
    Eigen::Matrix<double, 6, 6> K2 = ((M1 + M2) / (L * L)) * (rVec * zVec.transpose() + zVec * rVec.transpose());

    return KM + K1 + K2;
}

Eigen::VectorXd Viga2DCorrotacional::getForcasInternasGlobais(const Eigen::VectorXd& uGlobal) const {
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

    double N = mat.E * mat.A * ul / L0;
    double constanteFlexao = 2.0 * mat.E * mat.I / L0;
    double M1 = constanteFlexao * (2.0 * teta1 + teta2);
    double M2 = constanteFlexao * (teta1 + 2.0 * teta2);

    Eigen::Vector3d fLocal(N, M1, M2);

    Eigen::Matrix<double, 3, 6> B;
    B << -C,   -S,   0.0,  C,    S,   0.0,
         -S/L,  C/L, 1.0,  S/L, -C/L, 0.0,
         -S/L,  C/L, 0.0,  S/L, -C/L, 1.0;

    return B.transpose() * fLocal;
}

EsforcosLocais Viga2DCorrotacional::getEsforcosLocais(const Eigen::VectorXd& uGlobal) const {
    std::vector<int> GDLs = getGDLsGlobais();
    Eigen::Vector<double, 6> uGlobalElem;
    for (int i = 0; i < 6; ++i) uGlobalElem(i) = uGlobal(GDLs[i]);

    double X1 = n1->x; double Y1 = n1->y;
    double X2 = n2->x; double Y2 = n2->y;

    double dx = (X2 + uGlobalElem(3)) - (X1 + uGlobalElem(0));
    double dy = (Y2 + uGlobalElem(4)) - (Y1 + uGlobalElem(1));
    double L = std::sqrt(dx * dx + dy * dy);
    if (L < 1e-12) L = 1e-12;

    double beta0 = std::atan2(Y2 - Y1, X2 - X1);
    double beta = std::atan2(dy, dx);

    double teta1 = normalizaAngulo(uGlobalElem(2) + beta0 - beta);
    double teta2 = normalizaAngulo(uGlobalElem(5) + beta0 - beta);
    double ul = (L * L - L0 * L0) / (L + L0);

    double N = mat.E * mat.A * ul / L0;
    double constanteFlexao = 2.0 * mat.E * mat.I / L0;
    double M1 = constanteFlexao * (2.0 * teta1 + teta2);
    double M2 = constanteFlexao * (teta1 + 2.0 * teta2);
    
    double V = (M1 + M2) / L; 

    return {N, V, M1, M2};
}

Eigen::MatrixXd Viga2DCorrotacional::getMatrizRigidezGeometrica(double N) const {
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

} // namespace Engine
