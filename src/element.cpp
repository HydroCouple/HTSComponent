#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "htsmodel.h"

#include <math.h>

using namespace std;

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  HTSModel *model)
  : id(id),
    mainChannelTemperature(0),
    groundTemperature(0),
    numSolutes(0),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    mainChannelSoluteConcs(nullptr),
    groundSoluteConcs(nullptr),
    sedSoluteDiffCoefficients(nullptr),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    length(0.0),
    depth(0.001),
    xSectionArea(0.0),
    width(0.001),
    mainChannelConductionHeat(0.0),
    groundConductionHeat(0.0),
    advectionHeat(0.0),
    externalHeatFluxes(0.0),
    radiationFluxes(0.0),
    externalSoluteFluxes(nullptr),
    totalHeatBalance(0.0),
    totalRadiationFluxesHeatBalance(0.0),
    totalExternalHeatFluxesBalance(0.0),
    totalSoluteMassBalance(nullptr),
    totalExternalSoluteFluxesMassBalance(nullptr),
    groundConductionDepth(0.01),
    model(model)
{

  initializeSolutes();

  upstream->outgoingElements.insert(this);
  downstream->incomingElements.insert(this);

  x = (upstream->x +  downstream->x) / 2.0;
  y = (upstream->y +  downstream->y) / 2.0;
  z = (upstream->z +  downstream->z) / 2.0;

}

Element::~Element()
{
  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
    delete[] externalSoluteFluxes;
    delete[] mainChannelSoluteConcs; mainChannelSoluteConcs = nullptr;
    delete[] groundSoluteConcs; groundSoluteConcs = nullptr;
    delete[] sedSoluteDiffCoefficients;  sedSoluteDiffCoefficients = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;
    delete[] totalExternalSoluteFluxesMassBalance; totalExternalSoluteFluxesMassBalance = nullptr;
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

}

void Element::initialize()
{
  //set upstream and downstream elements
  mainChannelConductionHeat = 0.0;
  groundConductionHeat = 0.0;
  advectionHeat = 0.0;

  totalHeatBalance =  totalRadiationFluxesHeatBalance = 0.0;

  for(int i = 0; i < numSolutes; i++)
  {
    totalSoluteMassBalance[i] = 0.0;
    totalExternalSoluteFluxesMassBalance[i] = 0.0;
  }
}

void Element::initializeSolutes()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] mainChannelSoluteConcs; mainChannelSoluteConcs = nullptr;
    delete[] groundSoluteConcs; groundSoluteConcs = nullptr;
    delete[] sedSoluteDiffCoefficients;  sedSoluteDiffCoefficients = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;
    delete[] totalExternalSoluteFluxesMassBalance; totalExternalSoluteFluxesMassBalance = nullptr;
  }

  if(model->m_solutes.size() > 0)
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
    mainChannelSoluteConcs = new double[numSolutes]();
    groundSoluteConcs = new double[numSolutes];
    sedSoluteDiffCoefficients = new double[numSolutes]();
    totalSoluteMassBalance = new double[numSolutes]();
    totalExternalSoluteFluxesMassBalance = new double[numSolutes]();
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0;
  xSectionArea = width * depth;
  double volume = xSectionArea * length;
  double rho_cp_vol = model->m_sedDensity * model->m_sedCp * volume;
  double currentTemp = T[index];

  if(volume)
  {
    //conduction main channel
    mainChannelConductionHeat = sedThermalDiffCoefficient * width * length * (mainChannelTemperature - currentTemp) *
                                model->m_sedDensity * model->m_sedCp / depth;

    //conduction ground
    groundConductionHeat = sedThermalDiffCoefficient * width * length * (groundTemperature - currentTemp) * model->m_sedDensity * model->m_sedCp / groundConductionDepth;

    //advection
    advectionHeat = model->m_waterDensity * model->m_cp * mainChannelAdvectionCoeff * (mainChannelTemperature - currentTemp);

    //Sum all heat
    DTDt = (mainChannelConductionHeat + groundConductionHeat + advectionHeat + externalHeatFluxes + (radiationFluxes * length * width)) / rho_cp_vol;
  }

  return DTDt;
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;
  xSectionArea = width * depth;
  double volume = xSectionArea * length;
  double rho_vol = model->m_sedDensity * volume;

  //conduction main channel
  DSoluteDt += sedSoluteDiffCoefficients[soluteIndex] * width * length * (mainChannelSoluteConcs[soluteIndex] - S[index]) / (depth * volume);

  //conduction ground
  DSoluteDt += sedSoluteDiffCoefficients[soluteIndex] * width * length * (groundSoluteConcs[soluteIndex] - S[index]) / (groundConductionDepth * volume);

  //advection
  DSoluteDt += model->m_waterDensity * mainChannelAdvectionCoeff * (mainChannelSoluteConcs[soluteIndex] - S[index]) / rho_vol;

  //External heat sources
  DSoluteDt += externalSoluteFluxes[soluteIndex] / rho_vol;

  return DSoluteDt;
}

double Element::computeDispersionFactor() const
{
  double maxFactor = length * xSectionArea / mainChannelAdvectionCoeff;

  if(depth)
  {
    maxFactor   = max(maxFactor,  2.0 * sedThermalDiffCoefficient / (depth * depth));

    for(int i = 0; i < numSolutes; i++)
      maxFactor   = max(maxFactor,  2.0 * sedSoluteDiffCoefficients[i] / (depth * depth));
  }

  if(groundConductionDepth)
  {
    maxFactor   = max(maxFactor,  2.0 * sedThermalDiffCoefficient / (groundConductionDepth * groundConductionDepth));

    for(int i = 0; i < numSolutes; i++)
    {
      maxFactor   = max(maxFactor,  2.0 * sedSoluteDiffCoefficients[i] / (groundConductionDepth * groundConductionDepth));
    }
  }

  return maxFactor;
}

void Element::computeHeatBalance(double timeStep)
{
  double radiationEnergy = radiationFluxes * length * width * timeStep / 1000.0;
  totalRadiationFluxesHeatBalance += radiationEnergy;

  double externalEnergy = externalHeatFluxes * xSectionArea * length * timeStep / 1000.0;
  totalExternalHeatFluxesBalance += externalEnergy;

  double totalHeatEnergy = model->m_waterDensity * model->m_cp * xSectionArea * length * (temperature.value - prevTemperature.value) / 1000.0;
  totalHeatBalance +=  totalHeatEnergy;

}

void Element::computeSoluteBalance(double timeStep, int soluteIndex)
{
  double externalMass = externalSoluteFluxes[soluteIndex] * xSectionArea * length * timeStep;
  totalExternalSoluteFluxesMassBalance[soluteIndex] += externalMass;

  double totalMass = model->m_waterDensity * xSectionArea * length * (soluteConcs[soluteIndex].value - prevSoluteConcs[soluteIndex].value) ;
  totalSoluteMassBalance[soluteIndex] +=  totalMass;
}
