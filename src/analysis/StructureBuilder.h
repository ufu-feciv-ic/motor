#pragma once
#include <memory>
#include <Eigen/Core>
#include "StructuralModel.h"
#include "DOFManager.h"
#include "FiniteElementCore.h"

class StructureBuilder {
public:
    static void construir(const StructuralModel& model, DOFManager& outDofManager, EstruturaModerna& outEst);
    static Eigen::VectorXd montarVetorCargas(int totalDofs, const StructuralModel& model, const DOFManager& dofManager);
};
