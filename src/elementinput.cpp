#include "stdafx.h"
#include "htscomponent.h"
#include "elementinput.h"
#include "spatial/point.h"
#include "spatial/linestring.h"
#include "htsmodel.h"
#include "element.h"
#include "temporal/timedata.h"

#include <QDebug>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace HydroCouple::Temporal;
using namespace HydroCouple::SpatioTemporal;


ElementInput::ElementInput(const QString &id,
                           Dimension *timeDimension,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           VariableType varType,
                           HTSComponent *modelComponent)
  : TimeGeometryMultiInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
                                 valueDefinition, modelComponent),
    m_component(modelComponent),
    m_varType(varType),
    m_soluteIndex(0)
{

}

bool ElementInput::addProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::addProvider(provider))
  {
    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;
    ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;
    IIdBasedComponentDataItem *idBasedDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
       timeGeometryDataItem->geometryCount())
    {
      std::vector<bool> mapped(timeGeometryDataItem->geometryCount(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

        if(lineString->pointCount())
        {
          HCPoint *p1 = lineString->pointInternal(0);
          HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

          for(int j = 0; j < timeGeometryDataItem->geometryCount() ; j++)
          {
            if(!mapped[j])
            {
              ILineString *lineStringProvider = dynamic_cast<ILineString*>(timeGeometryDataItem->geometry(j));

              IPoint *pp1 = lineStringProvider->point(0);
              IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);

              double deltap1p1 = hypot(p1->x() - pp1->x() , p1->y() - pp1->y());
              double deltap2p2 = hypot(p2->x() - pp2->x() , p2->y() - pp2->y());

              double deltap1p2 = hypot(p1->x() - pp2->x() , p1->y() - pp2->y());
              double deltap2p1 = hypot(p2->x() - pp1->x() , p2->y() - pp1->y());

              if( deltap1p1 < 1e-3 && deltap2p2 < 1e-3)
              {
                m_geometryMapping[provider][i] = j;
                m_geometryMappingOrientation[provider][i] = 1.0;
                mapped[j] = true;
                break;
              }
              else if(deltap1p2 < 1e-3 &&  deltap2p1 < 1e-3)
              {
                m_geometryMapping[provider][i] = j;
                m_geometryMappingOrientation[provider][i] = -1.0;
                mapped[j] = true;
                break;
              }
            }
          }
        }
      }
    }
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)) &&
            geometryDataItem->geometryCount())
    {
      std::vector<bool> mapped(geometryDataItem->geometryCount(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

        if(lineString->pointCount())
        {
          HCPoint *p1 = lineString->pointInternal(0);
          HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

          for(int j = 0; j < geometryDataItem->geometryCount() ; j++)
          {
            if(!mapped[j])
            {
              ILineString *lineStringProvider = dynamic_cast<ILineString*>(geometryDataItem->geometry(j));

              IPoint *pp1 = lineStringProvider->point(0);
              IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);

              double deltap1p1 = hypot(p1->x() - pp1->x() , p1->y() - pp1->y());
              double deltap2p2 = hypot(p2->x() - pp2->x() , p2->y() - pp2->y());

              double deltap1p2 = hypot(p1->x() - pp2->x() , p1->y() - pp2->y());
              double deltap2p1 = hypot(p2->x() - pp1->x() , p2->y() - pp1->y());

              if( deltap1p1 < 1e-3 && deltap2p2 < 1e-3)
              {
                m_geometryMapping[provider][i] = j;
                m_geometryMappingOrientation[provider][i] = 1.0;
                mapped[j] = true;
                break;
              }
              else if(deltap1p2 < 1e-3 &&  deltap2p1 < 1e-3)
              {
                m_geometryMapping[provider][i] = j;
                m_geometryMappingOrientation[provider][i] = -1.0;
                mapped[j] = true;
                break;
              }
            }
          }
        }
      }
    }
    else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)))
    {
      QStringList identifiers = timeIdBasedDataItem->identifiers();

      std::vector<bool> mapped(identifiers.size(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        Element *element = m_component->modelInstance()->getElement(i);

        for(int j = 0; j < identifiers.size() ; j++)
        {
          if(!mapped[j])
          {
            QString providerId = identifiers[j];

            if(element->id == providerId.toStdString())
            {
              m_geometryMapping[provider][i] = j;
              m_geometryMappingOrientation[provider][i] = 1.0;
              mapped[j] = true;
              break;
            }
          }
        }
      }
    }
    else if((idBasedDataItem = dynamic_cast<IIdBasedComponentDataItem*>(provider)))
    {
      QStringList identifiers = idBasedDataItem->identifiers();

      std::vector<bool> mapped(identifiers.size(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        Element *element = m_component->modelInstance()->getElement(i);

        for(int j = 0; j < identifiers.size() ; j++)
        {
          if(!mapped[j])
          {
            QString providerId = identifiers[j];

            if(element->id == providerId.toStdString())
            {
              m_geometryMapping[provider][i] = j;
              m_geometryMappingOrientation[provider][i] = 1.0;
              mapped[j] = true;
              break;
            }
          }
        }
      }
    }

    return true;
  }

  return false;
}

bool ElementInput::removeProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::removeProvider(provider))
  {
    m_geometryMapping.erase(provider);
    m_geometryMappingOrientation.erase(provider);
    return true;
  }

  return false;
}

bool ElementInput::canConsume(HydroCouple::IOutput *provider, QString &message) const
{
  ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
  IGeometryComponentDataItem *geometryDataItem = nullptr;
  ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;
  IIdBasedComponentDataItem *idBasedDataItem = nullptr;

  if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
     (timeGeometryDataItem->geometryType() == IGeometry::LineString ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZ ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZM) &&
     (provider->valueDefinition()->type() == QVariant::Double ||
      provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)) &&
          (geometryDataItem->geometryType() == IGeometry::LineString ||
           geometryDataItem->geometryType() == IGeometry::LineStringZ ||
           geometryDataItem->geometryType() == IGeometry::LineStringZM) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((idBasedDataItem = dynamic_cast<IIdBasedComponentDataItem*>(provider)) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }

  message = "Provider must be a LineString";

  return false;
}

void ElementInput::retrieveValuesFromProvider()
{
  moveDataToPrevTime();
  int currentTimeIndex = m_times.size() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_component->modelInstance()->currentDateTime());

  for(HydroCouple::IOutput *provider : providers())
  {
    provider->updateValues(this);
  }
}

void ElementInput::applyData()
{
  double currentTime = m_component->modelInstance()->currentDateTime();

  for(HydroCouple::IOutput *provider : providers())
  {
    std::unordered_map<int,int>& geomap = m_geometryMapping[provider];
    std::unordered_map<int,double>& geoorient = m_geometryMappingOrientation[provider];

    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;
    ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;
    IIdBasedComponentDataItem *idBasedDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)))
    {
      int currentTimeIndex = timeGeometryDataItem->timeCount() - 1;
      int previousTimeIndex = std::max(0 , timeGeometryDataItem->timeCount() - 2);

      double providerCurrentTime = timeGeometryDataItem->time(currentTimeIndex)->julianDay();
      double providerPreviousTime = timeGeometryDataItem->time(previousTimeIndex)->julianDay();

      if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
      {

        double factor = 0.0;

        if(providerCurrentTime > providerPreviousTime)
        {
          double denom = providerCurrentTime - providerPreviousTime;
          double numer = currentTime - providerPreviousTime;
          factor = numer / denom;
        }

        switch (m_varType)
        {
          case MainChannelTemperature:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                double orientation = geoorient[it.first];
                value1 *= orientation;
                value2 *= orientation;

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelTemperature = value2 + factor *(value1 - value2);

              }
            }
            break;
          case MainChannelSoluteConc:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                double orientation = geoorient[it.first];
                value1 *= orientation;
                value2 *= orientation;

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelSoluteConcs[m_soluteIndex] = value2 + factor *(value1 - value2);

              }
            }
            break;
          case Width:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->width = value2 + factor *(value1 - value2);
              }
            }
            break;
          case QHTS:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelAdvectionCoeff = value2 + factor *(value1 - value2);
              }
            }
            break;
          case YHTS:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->depth = value2 + factor *(value1 - value2);
              }
            }
            break;
          case AlphaSed:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->sedThermalDiffCoefficient = value2 + factor *(value1 - value2);
              }
            }
            break;
          case GroundTemperature:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundTemperature = value2 + factor *(value1 - value2);
              }
            }
            break;
          case GroundConductionDepth:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundConductionDepth = value2 + factor *(value1 - value2);
              }
            }
            break;
        }
      }
      else
      {
        switch (m_varType)
        {
          case MainChannelTemperature:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                value *= geoorient[it.first];
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelTemperature = value;
              }
            }
            break;
          case MainChannelSoluteConc:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                value *= geoorient[it.first];
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelSoluteConcs[m_soluteIndex] = value;
              }
            }
            break;
          case Width:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->width = value;
              }
            }
            break;
          case QHTS:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelAdvectionCoeff = value;
              }
            }
            break;
          case YHTS:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->depth = value;
              }
            }
            break;
          case AlphaSed:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->sedThermalDiffCoefficient = value;
              }
            }
            break;
          case GroundTemperature:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundTemperature = value;
              }
            }
            break;
          case GroundConductionDepth:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundConductionDepth = value;
              }
            }
            break;
        }
      }
    }
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)))
    {
      switch (m_varType)
      {
        case MainChannelTemperature:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              value *= geoorient[it.first];
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->mainChannelTemperature = value;
            }
          }
          break;
        case MainChannelSoluteConc:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              value *= geoorient[it.first];
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->mainChannelSoluteConcs[m_soluteIndex] = value;
            }
          }
          break;
        case Width:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->width = value;
            }
          }
          break;
        case QHTS:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->mainChannelAdvectionCoeff = value;
            }
          }
          break;
        case YHTS:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->depth = value;
            }
          }
          break;
        case AlphaSed:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->sedThermalDiffCoefficient = value;
            }
          }
          break;
        case GroundTemperature:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->groundTemperature = value;
            }
          }
          break;
        case GroundConductionDepth:
          {
            for(auto it : geomap)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->groundConductionDepth = value;
            }
          }
          break;
      }
    }
    else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)))
    {
      int currentTimeIndex = timeIdBasedDataItem->timeCount() - 1;
      int previousTimeIndex = std::max(0 , timeIdBasedDataItem->timeCount() - 2);

      double providerCurrentTime = timeIdBasedDataItem->time(currentTimeIndex)->julianDay();
      double providerPreviousTime = timeIdBasedDataItem->time(previousTimeIndex)->julianDay();

      if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
      {

        double factor = 0.0;

        if(providerCurrentTime > providerPreviousTime)
        {
          double denom = providerCurrentTime - providerPreviousTime;
          double numer = currentTime - providerPreviousTime;
          factor = numer / denom;
        }

        switch (m_varType)
        {
          case MainChannelTemperature:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                double orientation = geoorient[it.first];
                value1 *= orientation;
                value2 *= orientation;

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelTemperature = value2 + factor *(value1 - value2);

              }
            }
            break;
          case MainChannelSoluteConc:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                double orientation = geoorient[it.first];
                value1 *= orientation;
                value2 *= orientation;

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelSoluteConcs[m_soluteIndex] = value2 + factor *(value1 - value2);

              }
            }
            break;
          case Width:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->width = value2 + factor *(value1 - value2);
              }
            }
            break;
          case QHTS:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelAdvectionCoeff = value2 + factor *(value1 - value2);
              }
            }
            break;
          case YHTS:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->depth = value2 + factor *(value1 - value2);
              }
            }
            break;
          case AlphaSed:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->sedThermalDiffCoefficient = value2 + factor *(value1 - value2);
              }
            }
            break;
          case GroundTemperature:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundTemperature = value2 + factor *(value1 - value2);
              }
            }
            break;
          case GroundConductionDepth:
            {
              for(auto it : geomap)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundConductionDepth = value2 + factor *(value1 - value2);
              }
            }
            break;
        }
      }
      else
      {
        switch (m_varType)
        {
          case MainChannelTemperature:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                value *= geoorient[it.first];
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelTemperature = value;
              }
            }
            break;
          case MainChannelSoluteConc:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                value *= geoorient[it.first];
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelSoluteConcs[m_soluteIndex] = value;
              }
            }
            break;
          case Width:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->width = value;
              }
            }
            break;
          case QHTS:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->mainChannelAdvectionCoeff = value;
              }
            }
            break;
          case YHTS:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->depth = value;
              }
            }
            break;
          case AlphaSed:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->sedThermalDiffCoefficient = value;
              }
            }
            break;
          case GroundTemperature:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundTemperature = value;
              }
            }
            break;
          case GroundConductionDepth:
            {
              for(auto it : geomap)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->groundConductionDepth = value;
              }
            }
            break;
        }
      }
    }
    else if((idBasedDataItem = dynamic_cast<IIdBasedComponentDataItem*>(provider)))
    {
      switch (m_varType)
      {
        case MainChannelTemperature:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              value *= geoorient[it.first];
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->mainChannelTemperature = value;
            }
          }
          break;
        case MainChannelSoluteConc:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              value *= geoorient[it.first];
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->mainChannelSoluteConcs[m_soluteIndex] = value;
            }
          }
          break;
        case Width:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->width = value;
            }
          }
          break;
        case QHTS:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->mainChannelAdvectionCoeff = value;
            }
          }
          break;
        case YHTS:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->depth = value;
            }
          }
          break;
        case AlphaSed:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->sedThermalDiffCoefficient = value;
            }
          }
          break;
        case GroundTemperature:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->groundTemperature = value;
            }
          }
          break;
        case GroundConductionDepth:
          {
            for(auto it : geomap)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->groundConductionDepth = value;
            }
          }
          break;
      }
    }
  }
}

ElementInput::VariableType ElementInput::variableType() const
{
  return m_varType;
}

void ElementInput::setVariableType(VariableType variableType)
{
  m_varType = variableType;
}

int ElementInput::soluteIndex() const
{
  return m_soluteIndex;
}

void ElementInput::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}


ElementHeatSourceInput::ElementHeatSourceInput(const QString &id,
                                               Dimension *timeDimension,
                                               Dimension *geometryDimension,
                                               ValueDefinition *valueDefinition,
                                               SourceType srcType,
                                               HTSComponent *modelComponent)
  : TimeGeometryMultiInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
                                 valueDefinition, modelComponent),
    m_component(modelComponent),
    m_srcType(srcType)
{

}

bool ElementHeatSourceInput::addProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::addProvider(provider))
  {
    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)))
    {
      std::unordered_map<int, int> geometryMapping;

      if(timeGeometryDataItem->geometryCount())
      {
        std::vector<bool> mapped(timeGeometryDataItem->geometryCount(), false);

        for(int i = 0; i < geometryCount() ; i++)
        {
          HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

          if(lineString->pointCount())
          {
            HCPoint *p1 = lineString->pointInternal(0);
            HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

            for(int j = 0; j < timeGeometryDataItem->geometryCount() ; j++)
            {
              if(!mapped[j])
              {
                ILineString *lineStringProvider = dynamic_cast<ILineString*>(timeGeometryDataItem->geometry(j));

                IPoint *pp1 = lineStringProvider->point(0);
                IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);


                if(hypot(p1->x() - pp1->x() , p1->y() - pp1->y()) < 1e-3 && hypot(p2->x() - pp2->x() , p2->y() - pp2->y()) < 1e-3)
                {
                  geometryMapping[i] = j;
                  mapped[j] = true;
                  break;
                }
                else if(hypot(p1->x() - pp2->x() , p1->y() - pp2->y()) < 1e-3 && hypot(p2->x() - pp1->x() , p2->y() - pp1->y()) < 1e-3)
                {
                  geometryMapping[i] = j;
                  mapped[j] = true;
                  break;
                }
              }
            }
          }
        }
      }

      m_geometryMapping[provider] = geometryMapping;

      return true;
    }
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)))
    {
      std::unordered_map<int, int> geometryMapping;

      if(geometryDataItem->geometryCount())
      {
        std::vector<bool> mapped(geometryDataItem->geometryCount(), false);

        for(int i = 0; i < geometryCount() ; i++)
        {
          HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

          if(lineString->pointCount())
          {
            HCPoint *p1 = lineString->pointInternal(0);
            HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

            for(int j = 0; j < geometryDataItem->geometryCount() ; j++)
            {
              if(!mapped[j])
              {
                ILineString *lineStringProvider = dynamic_cast<ILineString*>(geometryDataItem->geometry(j));

                IPoint *pp1 = lineStringProvider->point(0);
                IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);


                if(hypot(p1->x() - pp1->x() , p1->y() - pp1->y()) < 1e-3 && hypot(p2->x() - pp2->x() , p2->y() - pp2->y()) < 1e-3)
                {
                  geometryMapping[i] = j;
                  mapped[j] = true;
                  break;
                }
                else if(hypot(p1->x() - pp2->x() , p1->y() - pp2->y()) < 1e-3 && hypot(p2->x() - pp1->x() , p2->y() - pp1->y()) < 1e-3)
                {
                  geometryMapping[i] = j;
                  mapped[j] = true;
                  break;
                }
              }
            }
          }
        }
      }

      m_geometryMapping[provider] = geometryMapping;

      return true;
    }
  }

  return false;
}

bool ElementHeatSourceInput::removeProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::removeProvider(provider))
  {
    m_geometryMapping.erase(provider);
    return true;
  }

  return false;
}

bool ElementHeatSourceInput::canConsume(HydroCouple::IOutput *provider, QString &message) const
{
  ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
  IGeometryComponentDataItem *geometryDataItem = nullptr;

  if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
     (timeGeometryDataItem->geometryType() == IGeometry::LineString ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZ ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZM) &&
     (provider->valueDefinition()->type() == QVariant::Double ||
      provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)) &&
          (geometryDataItem->geometryType() == IGeometry::LineString ||
           geometryDataItem->geometryType() == IGeometry::LineStringZ ||
           geometryDataItem->geometryType() == IGeometry::LineStringZM) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }

  message = "Provider must be a LineString";

  return false;
}

void ElementHeatSourceInput::retrieveValuesFromProvider()
{
  moveDataToPrevTime();
  int currentTimeIndex = m_times.size() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_component->modelInstance()->currentDateTime());

  for(QList<IOutput*>::iterator it = m_providers.begin(); it !=  m_providers.end() ; it++)
  {
    (*it)->updateValues(this);
  }
}

void ElementHeatSourceInput::applyData()
{
  double currentTime = m_component->modelInstance()->currentDateTime();

  for(IOutput *provider : m_providers)
  {

    std::unordered_map<int,int> &geometryMapping = m_geometryMapping[provider];

    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)))
    {
      int currentTimeIndex = timeGeometryDataItem->timeCount() - 1;
      int previousTimeIndex = std::max(0 , timeGeometryDataItem->timeCount() - 2);

      double providerCurrentTime = timeGeometryDataItem->time(currentTimeIndex)->julianDay();
      double providerPreviousTime = timeGeometryDataItem->time(previousTimeIndex)->julianDay();

      if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
      {
        double factor = 0.0;

        if(providerCurrentTime > providerPreviousTime)
        {
          double denom = providerCurrentTime - providerPreviousTime;
          double numer =currentTime - providerPreviousTime;
          factor = numer / denom;
        }

        switch (m_srcType)
        {
          case RadiativeFlux:
            {
              for(auto it : geometryMapping)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->radiationFluxes += value2 + factor *(value1 - value2);

              }
            }
            break;
          case HeatFlux:
            {
              for(auto it : geometryMapping)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->externalHeatFluxes += value2 + factor *(value1 - value2);
              }
            }
            break;
        }
      }
      else
      {
        switch (m_srcType)
        {
          case RadiativeFlux:
            {
              for(auto it : geometryMapping)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->radiationFluxes += value;

              }
            }
            break;
          case HeatFlux:
            {
              for(auto it : geometryMapping)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->externalHeatFluxes += value;
              }
            }
            break;
        }
      }
    }
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)))
    {
      switch (m_srcType)
      {
        case RadiativeFlux:
          {
            for(auto it : geometryMapping)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, &value);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->radiationFluxes += value;

            }
          }
          break;
        case HeatFlux:
          {
            for(auto it : geometryMapping)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, &value);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->externalHeatFluxes += value;
            }
          }
          break;
      }
    }
  }
}

ElementHeatSourceInput::SourceType ElementHeatSourceInput::sourceType() const
{
  return m_srcType;
}

void ElementHeatSourceInput::setSourceType(SourceType srcType)
{
  m_srcType = srcType;
}
