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

    enum Bank
    {
      Left,
      Right,
      Lumped
    };

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
    * \brief temperature
    */
   Variable temperature;

   /*!
    * \brief prevTemperature
    */
   Variable prevTemperature;

   /*!
    * \brief groundTemperature
    */
   double groundTemperature;

   /*!
    * \brief mainChannelTemperature
    */
   double channelTemperature;

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
    * \brief mainChannelSoluteConcs
    */
   double *channelSoluteConcs;

   /*!
    * \brief sedimentDensity
    */
   double sedimentDensity;

   /*!
    * \brief sedimentSpecificHeatCapacity
    */
   double sedimentSpecificHeatCapacity;

   /*!
    * \brief coefficientAdvectiveTransport
    */
   double coefficientAdvectiveTransport;

   /*!
    * \brief sedimentCoefficientofThermalDiffusivity
    */
   double sedimentCoefficientofThermalDiffusivity;

   /*!
    * \brief groundConductionDepth
    */
   double groundConductionDepth;

   /*!
    * \brief fromJunction
    */
   ElementJunction *upstreamJunction;

   /*!
    * \brief toJunction
    */
   ElementJunction *downstreamJunction;

   /*!
    * \brief length
    */
   double length;

   /*!
    * \brief depth
    */
   double depth;

   /*!
    * \brief xSectionArea
    */
   double xSectionArea;

   /*!
    * \brief width
    */
   double width;

   /*!
    * \brief beta- Fraction of the channel width that acts as a surface transient zone.
    */
   double beta;

   Bank bank;

   /*!
    * \brief externalHeatFluxes of J / s
    */
   double externalHeatFluxes;

   /*!
    * \brief externalSoluteFluxes of the form m^3. C / s
    */
   double *externalSoluteFluxes;

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

};

#endif // ELEMENT_H
