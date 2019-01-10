
#include "stdafx.h"
#include "element.h"
#include "htsmodel.h"
#include "hydraulicsbc.h"

HydraulicsBC::HydraulicsBC(Element *element, int variableIndex, HTSModel *model)
  :AbstractTimeSeriesBC(model),
   m_element(element),
   m_variableIndex(variableIndex)
{

}

HydraulicsBC::~HydraulicsBC()
{

}

void HydraulicsBC::findAssociatedGeometries()
{

}

void HydraulicsBC::prepare()
{

}

void HydraulicsBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableIndex)
    {
      case 1:
        m_element->depth = value;
        break;
      case 2:
        m_element->width = value;
        break;
      case 3:
        m_element->mainChannelAdvectionCoeff = value;
        break;
      case 4:
        m_element->groundConductionDepth = value;
        break;
    }
  }
}

Element *HydraulicsBC::element() const
{
  return m_element;
}

void HydraulicsBC::setElement(Element *element)
{
  m_element = element;
}


HydraulicsBC::HydraulicsBC(Element *startElement, Element *endElement, int variableIndex, HTSModel *model)
  :AbstractTimeSeriesBC(model),
   m_startElement(startElement),
   m_endElement(endElement),
   m_variableIndex(variableIndex)
{

}

HydraulicsBC::~HydraulicsBC()
{

}

void HydraulicsBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void HydraulicsBC::prepare()
{

}

void HydraulicsBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    for(Element *element : m_profile)
    {
      switch (m_variableIndex)
      {
        case 1:
          element->depth = value;
          break;
        case 2:
          element->width = value;
          break;
        case 3:
          element->mainChannelAdvectionCoeff = value;
          break;
        case 4:
          element->groundConductionDepth = value;
          break;
      }
    }
  }
}

Element *HydraulicsBC::startElement() const
{
  return m_startElement;
}

void HydraulicsBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *HydraulicsBC::endElement() const
{
  return m_endElement;
}

void HydraulicsBC::setEndElement(Element *element)
{
  m_endElement = element;
}
