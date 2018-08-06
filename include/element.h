#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>

#include "variable.h"

struct Element;
struct ElementJunction;
class HTSModel;

/*!
 * \brief This struct represents the channel control volume
 */
struct Element
{

    /*!
    * \brief Element - Creates an instance of the control volume element used to represent a computational
    * element in a reach.
    * \param numSolutes - Number of solutes that are going to be transported in addition to temperature.
    * \param from - The upstream junction of this element.
    * \param to - The downstream junction of this element.
    * \param project
    */
   Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  HTSModel *model);

   /*!
    * \brief ~Element - Destructor for this class.
    */
   ~Element();

   /*!
    * \brief index unique identifier for element
    */
   int index;

   /*!
    * \brief id
    */
   std::string id;

   /*!
    * \brief x
    */
   double x;

   /*!
    * \brief y
    */
   double y;

   /*!
    * \brief z
    */
   double z;

   /*!
    * \brief temperature  (°C)
    */
   Variable temperature;

   /*!
    * \brief prevTemperature  (°C)
    */
   Variable prevTemperature;

   /*!
    * \brief mainChannelTemperature  (°C)
    */
   double mainChannelTemperature;

   /*!
    * \brief groundTemperature  (°C)
    */
   double groundTemperature;

   /*!
    * \brief numSolutes - Number of solutes
    */
   int numSolutes = 0;

   /*!
    * \brief soluteConcs [kg/m^3]
    */
   Variable *soluteConcs;

   /*!
    * \brief prevSoluteConcs [kg/m^3]
    */
   Variable *prevSoluteConcs;

   /*!
    * \brief mainChannelSoluteConcs [kg/m^3]
    */
   double *mainChannelSoluteConcs;

   /*!
    * \brief groundSoluteConcs [kg/m^3]
    */
   double *groundSoluteConcs;

   /*!
    * \brief temperatureExchangeCoefficient (m^2/s)
    */
   double sedThermalDiffCoefficient;

   /*!
    * \brief soluteExchangeCoefficients (m^2/s)
    */
   double *sedSoluteDiffCoefficients;

   /*!
    * \brief mainChannelAdvectionCoeff (m^3/s)
    */
   double mainChannelAdvectionCoeff;

   /*!
    * \brief fromJunction
    */
   ElementJunction *upstreamJunction;

   /*!
    * \brief toJunction
    */
   ElementJunction *downstreamJunction;

   /*!
    * \brief length (m)
    */
   double length;

   /*!
    * \brief depth  (m)
    */
   double depth;

   /*!
    * \brief xSectionArea (m^2)
    */
   double xSectionArea;

   /*!
    * \brief width  (m)
    */
   double width;

   /*!
    * \brief mainChannelConductionHeat (J/s)
    */
   double mainChannelConductionHeat;

   /*!
    * \brief groundConductionHeat (J/s)
    */
   double groundConductionHeat;

   /*!
    * \brief advectionHeat (J/s)
    */
   double advectionHeat;

   /*!
    * \brief externalHeatFluxes of (J/s)
    */
   double externalHeatFluxes;

   /*!
    * \brief externalSoluteFluxes (W/m^2)
    */
   double radiationFluxes;

   /*!
    * \brief externalSoluteFluxes of the form (kg/s)
    */
   double *externalSoluteFluxes;

   /*!
    * \brief heatBalance (KJ)
    */
   double totalHeatBalance;

   /*!
    * \brief totalRadiationHeatBalance (KJ)
    */
   double totalRadiationFluxesHeatBalance;

   /*!
    * \brief totalExternalHeatFluxesBalance (KJ)
    */
   double totalExternalHeatFluxesBalance;

   /*!
    * \brief soluteMassBalance (kg)
    */
   double *totalSoluteMassBalance;

   /*!
    * \brief totalExternalSoluteFluxesMassBalance  (kg)
    */
   double *totalExternalSoluteFluxesMassBalance;

   /*!
    * \brief groundConductionDepth (m)
    */
   double groundConductionDepth;

   /*!
    * \brief distanceFromUpStreamJunction
    */
   double distanceFromUpStreamJunction;

   /*!
    * \brief model
    */
   HTSModel *model;

   /*!
    * \brief initializeSolutes
    * \param numSolutes
    */
   void initialize();

   /*!
    * \brief initializeSolutes
    */
   void initializeSolutes();

   /*!
    * \brief computeDTDt - Computes the time derivative of temperature based on data generated by the ODE solver.
    * \param dt - The timestep over which to compute the solute gradient.
    * \param T - The temperature array for all elements.
    * \return
    */
   double computeDTDt(double dt, double T[]);

   /*!
    * \brief computeDSoluteDt
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDt(double dt, double S[], int soluteIndex);

   /*!
    * \brief commputeDispersionFactor 1/s
    * \return
    */
   double computeDispersionFactor() const;

   /*!
    * \brief computeHeatBalance
    */
   void computeHeatBalance(double timeStep);

   /*!
    * \brief computeSoluteBalance
    * \param soluteIndex
    */
   void computeSoluteBalance(double timeStep, int soluteIndex);
};

#endif // ELEMENT_H
