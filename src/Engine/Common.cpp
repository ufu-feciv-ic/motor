#include "Engine/Common.h"

namespace Engine {

double normalizaAngulo(double angulo) {
    double a = std::fmod(angulo + M_PI, 2.0 * M_PI);
    if (a < 0.0) a += 2 * M_PI;
    return a - M_PI;
}

} // namespace Engine
