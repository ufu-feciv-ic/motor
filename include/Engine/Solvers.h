#ifndef ENGINE_SOLVERS_H
#define ENGINE_SOLVERS_H

#include <vector>
#include <iostream>
#include "Engine/Estrutura.h"
#include "Engine/Common.h"

namespace Engine {

class AnaliseBuckling {
public:
    struct ModoBuckling {
        Eigen::VectorXd modo;
        double lambda;
    };

    static std::vector<ModoBuckling> executar(Estrutura& est);
};

class EstrategiaAnalise {
public:
    virtual ~EstrategiaAnalise() = default;
    virtual std::vector<Resultado> executar(Estrutura& est) = 0;
};

class AnaliseLinear : public EstrategiaAnalise {
public:
    std::vector<Resultado> executar(Estrutura& est) override;
};

class AnaliseNaoLinearNR : public EstrategiaAnalise {
private:
    int numPassos;
    int maxIter;
    double tol;

public:
    AnaliseNaoLinearNR(int passos = 10, int iteracoes = 50, double tolerancia = 1e-6);
    std::vector<Resultado> executar(Estrutura& est) override;
};

class AnaliseNaoLinearCompArco : public EstrategiaAnalise {
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
        double iterDesejadas = 5.0);
    std::vector<Resultado> executar(Estrutura& est) override;
};

} // namespace Engine

#endif // ENGINE_SOLVERS_H
