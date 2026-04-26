#define _USE_MATH_DEFINES
#include "Elements.h"
#include <cmath>

// --- VIGA CORROTACIONAL ---

Viga2DCorrotacionalModerna::Viga2DCorrotacionalModerna(std::shared_ptr<NoModerno> node1, std::shared_ptr<NoModerno> node2, const MaterialData& mat)
    : n1(node1), n2(node2), E(mat.E), A(mat.A), I(mat.I) 
{
    double dx = n2->x - n1->x;
    double dy = n2->y - n1->y;
    L0 = std::sqrt(dx * dx + dy * dy);
}

double Viga2DCorrotacionalModerna::normalizarAngulo(double angulo) {
    double a = std::fmod(angulo + M_PI, 2.0 * M_PI);
    if (a < 0.0) a += 2.0 * M_PI;
    return a - M_PI;
}

std::vector<int> Viga2DCorrotacionalModerna::getIndicesGlobais(const DOFManager& dofManager) const {
    std::vector<int> i1 = dofManager.obterIndicesGlobais(n1->id);
    std::vector<int> i2 = dofManager.obterIndicesGlobais(n2->id);
    i1.insert(i1.end(), i2.begin(), i2.end());
    return i1;
}

Eigen::MatrixXd Viga2DCorrotacionalModerna::getMatrizRigidez(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const {
    auto indices = getIndicesGlobais(dofManager);
    Eigen::VectorXd uElem(6);
    for(int i=0; i<6; ++i) uElem(i) = uGlobal(indices[i]);

    double dx = (n2->x + uElem(3)) - (n1->x + uElem(0));
    double dy = (n2->y + uElem(4)) - (n1->y + uElem(1));
    double L = std::sqrt(dx * dx + dy * dy);
    if (L < 1e-12) L = 1e-12;

    double C = dx / L; double S = dy / L;
    double beta0 = std::atan2(n2->y - n1->y, n2->x - n1->x);
    double beta = std::atan2(dy, dx);

    double teta1 = normalizarAngulo(uElem(2) + beta0 - beta);
    double teta2 = normalizarAngulo(uElem(5) + beta0 - beta);
    double ul = (L * L - L0 * L0) / (L + L0);

    double N = E * A * ul / L0;
    double M1 = (2.0 * E * I / L0) * (2.0 * teta1 + teta2);
    double M2 = (2.0 * E * I / L0) * (teta1 + 2.0 * teta2);

    Eigen::Matrix<double, 3, 6> B;
    B << -C, -S, 0, C, S, 0,
         -S/L, C/L, 1, S/L, -C/L, 0,
         -S/L, C/L, 0, S/L, -C/L, 1;

    Eigen::Matrix3d D;
    D << E*A/L0, 0, 0, 0, 4.0*E*I/L0, 2.0*E*I/L0, 0, 2.0*E*I/L0, 4.0*E*I/L0;

    Eigen::Matrix<double, 6, 6> KM = B.transpose() * D * B;
    Eigen::Vector<double, 6> r, z;
    r << -C, -S, 0, C, S, 0; z << S, -C, 0, -S, C, 0;

    Eigen::MatrixXd K1 = (N/L)*(z*z.transpose());
    Eigen::MatrixXd K2 = ((M1+M2)/(L*L))*(r*z.transpose() + z*r.transpose());

    return KM + K1 + K2;
}

Eigen::VectorXd Viga2DCorrotacionalModerna::getForcasInternas(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const {
    auto indices = getIndicesGlobais(dofManager);
    Eigen::VectorXd uElem(6);
    for(int i=0; i<6; ++i) uElem(i) = uGlobal(indices[i]);

    double dx = (n2->x + uElem(3)) - (n1->x + uElem(0));
    double dy = (n2->y + uElem(4)) - (n1->y + uElem(1));
    double L = std::sqrt(dx * dx + dy * dy);
    if (L < 1e-12) L = 1e-12;

    double C = dx / L; double S = dy / L;
    double beta0 = std::atan2(n2->y - n1->y, n2->x - n1->x);
    double beta = std::atan2(dy, dx);

    double teta1 = normalizarAngulo(uElem(2) + beta0 - beta);
    double teta2 = normalizarAngulo(uElem(5) + beta0 - beta);
    double ul = (L * L - L0 * L0) / (L + L0);

    double N = E * A * ul / L0;
    double M1 = (2.0 * E * I / L0) * (2.0 * teta1 + teta2);
    double M2 = (2.0 * E * I / L0) * (teta1 + 2.0 * teta2);

    Eigen::Matrix<double, 3, 6> B;
    B << -C, -S, 0, C, S, 0,
         -S/L, C/L, 1, S/L, -C/L, 0,
         -S/L, C/L, 0, S/L, -C/L, 1;

    return B.transpose() * Eigen::Vector3d(N, M1, M2);
}

EsforcosLocaisModernos Viga2DCorrotacionalModerna::getEsforcosLocais(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const {
    auto indices = getIndicesGlobais(dofManager);
    Eigen::VectorXd uElem(6);
    for(int i=0; i<6; ++i) uElem(i) = uGlobal(indices[i]);

    double dx = (n2->x + uElem(3)) - (n1->x + uElem(0));
    double dy = (n2->y + uElem(4)) - (n1->y + uElem(1));
    double L = std::sqrt(dx * dx + dy * dy);
    if (L < 1e-12) L = 1e-12;

    double beta0 = std::atan2(n2->y - n1->y, n2->x - n1->x);
    double beta = std::atan2(dy, dx);

    double teta1 = normalizarAngulo(uElem(2) + beta0 - beta);
    double teta2 = normalizarAngulo(uElem(5) + beta0 - beta);
    double ul = (L * L - L0 * L0) / (L + L0);

    double N = E * A * ul / L0;
    double M1 = (2.0 * E * I / L0) * (2.0 * teta1 + teta2);
    double M2 = (2.0 * E * I / L0) * (teta1 + 2.0 * teta2);
    double V = (M1 + M2) / L;

    return {N, V, M1, M2};
}

// --- VIGA LINEAR ---

Viga2DLinearModerna::Viga2DLinearModerna(std::shared_ptr<NoModerno> node1, std::shared_ptr<NoModerno> node2, const MaterialData& mat)
    : n1(node1), n2(node2), E(mat.E), A(mat.A), I(mat.I) 
{
    double dx = n2->x - n1->x;
    double dy = n2->y - n1->y;
    L0 = std::sqrt(dx * dx + dy * dy);
    cos0 = dx / L0;
    sin0 = dy / L0;
}

std::vector<int> Viga2DLinearModerna::getIndicesGlobais(const DOFManager& dofManager) const {
    std::vector<int> i1 = dofManager.obterIndicesGlobais(n1->id);
    std::vector<int> i2 = dofManager.obterIndicesGlobais(n2->id);
    i1.insert(i1.end(), i2.begin(), i2.end());
    return i1;
}

Eigen::MatrixXd Viga2DLinearModerna::getMatrizRigidez(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const {
    double EAL = E * A / L0;
    double EIL = E * I / L0;
    double EIL2 = E * I / (L0 * L0);
    double EIL3 = E * I / (L0 * L0 * L0);

    Eigen::MatrixXd kLoc = Eigen::MatrixXd::Zero(6, 6);
    kLoc << EAL, 0, 0, -EAL, 0, 0,
            0, 12*EIL3, 6*EIL2, 0, -12*EIL3, 6*EIL2,
            0, 6*EIL2, 4*EIL, 0, -6*EIL2, 2*EIL,
            -EAL, 0, 0, EAL, 0, 0,
            0, -12*EIL3, -6*EIL2, 0, 12*EIL3, -6*EIL2,
            0, 6*EIL2, 2*EIL, 0, -6*EIL2, 4*EIL;

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);
    T << cos0, sin0, 0, 0, 0, 0,
         -sin0, cos0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0,
         0, 0, 0, cos0, sin0, 0,
         0, 0, 0, -sin0, cos0, 0,
         0, 0, 0, 0, 0, 1;

    return T.transpose() * kLoc * T;
}

Eigen::VectorXd Viga2DLinearModerna::getForcasInternas(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const {
    auto indices = getIndicesGlobais(dofManager);
    Eigen::VectorXd uElem(6);
    for(int i=0; i<6; ++i) uElem(i) = uGlobal(indices[i]);

    Eigen::MatrixXd K = getMatrizRigidez(uGlobal, dofManager);
    return K * uElem;
}

EsforcosLocaisModernos Viga2DLinearModerna::getEsforcosLocais(const Eigen::VectorXd& uGlobal, const DOFManager& dofManager) const {
    auto indices = getIndicesGlobais(dofManager);
    Eigen::VectorXd uElem(6);
    for(int i=0; i<6; ++i) uElem(i) = uGlobal(indices[i]);

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);
    T << cos0, sin0, 0, 0, 0, 0,
         -sin0, cos0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0,
         0, 0, 0, cos0, sin0, 0,
         0, 0, 0, -sin0, cos0, 0,
         0, 0, 0, 0, 0, 1;

    Eigen::VectorXd uLoc = T * uElem;
    
    double EAL = E * A / L0;
    double EIL = E * I / L0;
    double EIL2 = E * I / (L0 * L0);
    double EIL3 = E * I / (L0 * L0 * L0);

    Eigen::MatrixXd kLoc = Eigen::MatrixXd::Zero(6, 6);
    kLoc << EAL, 0, 0, -EAL, 0, 0,
            0, 12*EIL3, 6*EIL2, 0, -12*EIL3, 6*EIL2,
            0, 6*EIL2, 4*EIL, 0, -6*EIL2, 2*EIL,
            -EAL, 0, 0, EAL, 0, 0,
            0, -12*EIL3, -6*EIL2, 0, 12*EIL3, -6*EIL2,
            0, 6*EIL2, 2*EIL, 0, -6*EIL2, 4*EIL;

    Eigen::VectorXd fLoc = kLoc * uLoc;
    return { fLoc(3), fLoc(4), fLoc(2), fLoc(5) };
}
