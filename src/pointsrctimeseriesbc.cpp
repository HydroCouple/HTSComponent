/*!
*  \file    pointsrctimeseriesbc.cpp
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
*  This file and its associated files and libraries are free software;
*  you can redistribute it and/or modify it under the terms of the
*  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
*  either version 3 of the License, or (at your option) any later version.
*  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
*  \date 2018
*  \pre
*  \bug
*  \todo Test transport on branching networks
*  \warning
*/

#include "pointsrctimeseriesbc.h"
#include "element.h"

PointSrcTimeSeriesBC::PointSrcTimeSeriesBC(Element *element, int variableIndex, HTSModel *model)
  :AbstractTimeSeriesBC(model),
   m_element(element),
   m_variableIndex(variableIndex)
{

}

PointSrcTimeSeriesBC::~PointSrcTimeSeriesBC()
{

}

void PointSrcTimeSeriesBC::findAssociatedGeometries()
{

}

void PointSrcTimeSeriesBC::prepare()
{

}

void PointSrcTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableIndex)
    {
      case -2:
        m_element->externalHeatFluxes += value;
        break;
      default:
        {
          m_element->externalSoluteFluxes[m_variableIndex]  += value;
        }
        break;
    }
  }
}

Element *PointSrcTimeSeriesBC::element() const
{
  return m_element;
}
