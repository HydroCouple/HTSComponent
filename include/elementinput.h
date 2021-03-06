/*!
 * \file elementinput.h
 * \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 * \version 1.0.0
 * \description
 * \license
 * This file and its associated files, and libraries are free software.
 * You can redistribute it and/or modify it under the terms of the
 * Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 * either version 3 of the License, or (at your option) any later version.
 * This file and its associated files is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 * \copyright Copyright 2014-2018, Caleb Buahin, All rights reserved.
 * \date 2014-2018
 * \pre
 * \bug
 * \warning
 * \todo
 */

#ifndef ELEMENTINPUT_H
#define ELEMENTINPUT_H

#include "htscomponent_global.h"
#include "spatiotemporal/timegeometrymultiinput.h"

#include <unordered_map>

class HTSComponent;


class HTSCOMPONENT_EXPORT ElementInput : public TimeGeometryMultiInputDouble
{
    Q_OBJECT

  public:

    enum VariableType
    {
      MainChannelTemperature,
      MainChannelSoluteConc,
      Width,
      QHTS,
      YHTS,
      AlphaSed,
      GroundTemperature,
      GroundConductionDepth,
    };

    ElementInput(const QString &id,
                 Dimension *timeDimension,
                 Dimension *geometryDimension,
                 ValueDefinition *valueDefinition,
                 VariableType varType,
                 HTSComponent *modelComponent);

    virtual ~ElementInput() override
    {

    }

    /*!
     * \brief setProvider
     * \param provider
     */
    bool addProvider(HydroCouple::IOutput *provider) override;

    /*!
     * \brief removeProvider
     * \param provider
     * \return
     */
    bool removeProvider(HydroCouple::IOutput *provider) override;

    /*!
     * \brief canConsume
     * \param provider
     * \param message
     * \return
     */
    bool canConsume(HydroCouple::IOutput *provider, QString &message) const override;

    /*!
     * \brief retrieveValuesFromProvider
     */
    void retrieveValuesFromProvider() override;

    /*!
     * \brief applyData
     */
    void applyData() override;

    /*!
     * \brief variableType
     * \return
     */
    VariableType variableType() const;

    /*!
     * \brief setVariableType
     * \param variableType
     */
    void setVariableType(VariableType variableType);

    /*!
     * \brief soluteIndex
     * \return
     */
    int soluteIndex() const;

    /*!
     * \brief setSoluteIndex
     * \param soluteIndex
     */
    void setSoluteIndex(int soluteIndex);

  private:

    std::unordered_map<HydroCouple::IOutput*, std::unordered_map<int,int>> m_geometryMapping;
    std::unordered_map<HydroCouple::IOutput*, std::unordered_map<int,double>> m_geometryMappingOrientation;
    HTSComponent *m_component;
    VariableType m_varType;
    int m_soluteIndex;

};

class HTSCOMPONENT_EXPORT ElementHeatSourceInput : public  TimeGeometryMultiInputDouble
{
    Q_OBJECT

  public:

    enum SourceType
    {
      RadiativeFlux,
      HeatFlux,
    };

    ElementHeatSourceInput(const QString &id,
                           Dimension *timeDimension,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           SourceType srcType,
                           HTSComponent *modelComponent);

    bool addProvider(HydroCouple::IOutput *provider) override;

    bool removeProvider(HydroCouple::IOutput *provider) override;

    bool canConsume(HydroCouple::IOutput *provider, QString &message) const override;

    void retrieveValuesFromProvider() override;

    void applyData() override;

    SourceType sourceType() const;

    void setSourceType(SourceType srcType);

  private:

    HTSComponent *m_component;
    SourceType m_srcType;
    std::unordered_map<HydroCouple::IOutput*, std::unordered_map<int,int>> m_geometryMapping;
};


#endif // ELEMENTINPUT_H
