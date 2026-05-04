#ifndef ENGINE_COMMON_H
#define ENGINE_COMMON_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "eigenpch.h"

namespace Engine {

double normalizaAngulo(double angulo);

struct PropriedadesMaterial {
    double E;
    double A;
    double I;
};

struct CargaDistribuida {
    double qx; // Carga axial (local x)
    double qy; // Carga transversal (local y)
};

struct EsforcosLocais {
    double N;  // Força normal
    double V;  // força cortante
    double M1; // Momento fletor no nó 1
    double M2; // Momento fletor no nó 2
};

struct Resultado {
    Eigen::VectorXd u;
    Eigen::VectorXd F;
    Eigen::VectorXd reacoes;
    double FatorCarga;
};

} // namespace Engine

#endif // ENGINE_COMMON_H
