#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "htsmodel.h"

#include <math.h>

using namespace std;

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  HTSModel *model)
  : id(id),
    groundTemperature(0.0),
    channelTemperature(0),
    channelSoluteConcs(nullptr),
    sedimentDensity(1.0),
    sedimentSpecificHeatCapacity(4187),
    coefficientAdvectiveTransport(0.0),
    sedimentCoefficientofThermalDiffusivity(0.0),
    groundConductionDepth(0.0),
    numSolutes(0),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    length(0.0),
    depth(0.0),
    width(0.0),
    beta(0.0),
    bank(Bank::Lumped),
    externalHeatFluxes(0.0),
    externalSoluteFluxes(nullptr),
    model(model)
{
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
    delete[] channelSoluteConcs; channelSoluteConcs = nullptr;
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

}

void Element::initialize()
{
  //set upstream and downstream elements
}

void Element::initializeSolutes()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] channelSoluteConcs; channelSoluteConcs = nullptr;
  }

  if(model->m_solutes.size() > 0)
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
    channelSoluteConcs = new double[numSolutes];
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0;
  double volume = xSectionArea * length;
  double rho_cp_vol = sedimentDensity * sedimentSpecificHeatCapacity * volume;

  //conduction main channel
  if(depth)
    DTDt += sedimentCoefficientofThermalDiffusivity * width * beta * length * (channelTemperature - T[index]) / depth / volume;

  //conduction ground
  if(groundConductionDepth)
    DTDt += sedimentCoefficientofThermalDiffusivity * width * beta * length * (groundTemperature - T[index]) / groundConductionDepth / volume;

  //advection
  DTDt += model->m_waterDensity * model->m_cp * coefficientAdvectiveTransport * (channelTemperature - T[index]) / rho_cp_vol;

  //External heat sources
  DTDt += externalHeatFluxes / rho_cp_vol;

  return DTDt;
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;
  double volume = xSectionArea * length;
  double rho_vol = sedimentDensity * volume;

  DSoluteDt += model->m_waterDensity * model->m_cp * coefficientAdvectiveTransport * (channelSoluteConcs[soluteIndex] - S[index]) / rho_vol;

  //Add external sources
  {
    DSoluteDt += externalSoluteFluxes[soluteIndex] / rho_vol;
  }

  return DSoluteDt;
}

double Element::computeDispersionFactor() const
{
  double maxFactor = length * xSectionArea / coefficientAdvectiveTransport;

  if(depth)
    maxFactor   = max(maxFactor,  2.0 * sedimentCoefficientofThermalDiffusivity / (depth * depth));

  if(groundConductionDepth)
    maxFactor   = max(maxFactor,  2.0 * sedimentCoefficientofThermalDiffusivity / (groundConductionDepth * groundConductionDepth));



  return maxFactor ;
}
