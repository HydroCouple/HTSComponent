/*!
*  \file    boundarycondition.cpp
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


#include "boundarycondition.h"
#include "element.h"
#include "htsmodel.h"

BoundaryCondition::BoundaryCondition(Element *element, int sourceIndex, int variableIndex, HTSModel *model)
  : AbstractTimeSeriesBC(model),
    m_element(element),
    m_sourceIndex(sourceIndex),
    m_variableIndex(variableIndex)
{

}

BoundaryCondition::~BoundaryCondition()
{

}

void BoundaryCondition::findAssociatedGeometries()
{

}

void BoundaryCondition::prepare()
{

}

void BoundaryCondition::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_sourceIndex)
    {
      case 1:
        {
          switch (m_variableIndex)
          {
            case -1:
              {
                m_element->mainChannelTemperature = value;
              }
              break;
            default:
              {
                m_element->mainChannelSoluteConcs[m_variableIndex] = value;
              }
              break;
          }
        }
        break;
      case 2:
        {
          switch (m_variableIndex)
          {
            case -1:
              {
                m_element->groundTemperature = value;
              }
              break;
            default:
              {
                m_element->groundSoluteConcs[m_variableIndex] = value;
              }
              break;
          }
        }
        break;
    }
  }
}

Element *BoundaryCondition::element() const
{
  return m_element;
}

void BoundaryCondition::setElement(Element *element)
{
  m_element = element;
}



UniformBoundaryCondition::UniformBoundaryCondition(Element *startElement, Element *endElement, int sourceIndex, int variableIndex, HTSModel *model)
  :AbstractTimeSeriesBC(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_sourceIndex(sourceIndex),
    m_variableIndex(variableIndex)
{

}

UniformBoundaryCondition::~UniformBoundaryCondition()
{

}

void UniformBoundaryCondition::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void UniformBoundaryCondition::prepare()
{

}

void UniformBoundaryCondition::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);


  if(found)
  {
    switch (m_sourceIndex)
    {
      case 1:
        {
          switch (m_variableIndex)
          {
            case -1:
              {
                for(Element *element : m_profile)
                {
                  element->mainChannelTemperature = value;
                }
              }
              break;
            default:
              {
                for(Element *element : m_profile)
                {
                  element->mainChannelSoluteConcs[m_variableIndex] = value;
                }
              }
              break;
          }
        }
        break;
      case 2:
        {
          switch (m_variableIndex)
          {
            case -1:
              {
                for(Element *element : m_profile)
                {
                  element->groundTemperature = value;
                }
              }
              break;
            default:
              {
                for(Element *element : m_profile)
                {
                  element->groundSoluteConcs[m_variableIndex] = value;
                }
              }
              break;
          }
        }
        break;
    }
  }
}

Element *UniformBoundaryCondition::startElement() const
{
  return m_startElement;
}

void UniformBoundaryCondition::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *UniformBoundaryCondition::endElement() const
{
  return m_endElement;
}

void UniformBoundaryCondition::setEndElement(Element *element)
{
  m_endElement = element;
}

