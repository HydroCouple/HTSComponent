/*!
 *  \file    htscomponent.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 *  This file and its associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2018
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */

#ifndef HTSCOMPONENT_H
#define HTSCOMPONENT_H

#include "htscomponent_global.h"
#include "htscomponentinfo.h"
#include "temporal/abstracttimemodelcomponent.h"

class Dimension;
class HTSModel;
class HCGeometry;
class ElementInput;
class ElementOutput;
class Unit;
class ElementHeatSourceInput;

class HTSCOMPONENT_EXPORT HTSComponent : public AbstractTimeModelComponent,
    public virtual HydroCouple::ICloneableModelComponent
{

    Q_OBJECT

    Q_INTERFACES(HydroCouple::ICloneableModelComponent)

  public:

    /*!
     * \brief HTSComponent constructor
     * \param id Unique identifier for this component instance.
     * \param modelComponentInfo the parent ModelComponentInfo that generated this component instance.
     */
    HTSComponent(const QString &id, HTSComponentInfo* modelComponentInfo = nullptr);

    /*!
     * \brief ~HTSComponent destructor
     */
    virtual ~HTSComponent();

    /*!
     * \brief validate validates this component model instance
     * \return Returns a list of error messages.
     */
    QList<QString> validate() override;

    /*!
     * \brief prepare Prepares the model component instance.
     */
    void prepare() override;

    /*!
     * \brief update
     * \param requiredOutputs
     */
    void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

    /*!
     * \brief finish
     */
    void finish() override;

    /*!
     * \brief modelInstance
     * \return
     */
    HTSModel *modelInstance() const;

    /*!
     * \brief parent
     * \return
     */
    HydroCouple::ICloneableModelComponent* parent() const override;

    /*!
     * \brief clone
     * \return
     */
    HydroCouple::ICloneableModelComponent* clone() override;

    /*!
     * \brief clones
     * \return
     */
    QList<HydroCouple::ICloneableModelComponent*> clones() const override;

  protected:

    bool removeClone(HTSComponent *component);


    /*!
     * \brief intializeFailureCleanUp
     */
    void initializeFailureCleanUp() override;

  private:

    /*!
     * \brief createArguments
     */
    void createArguments() override;

    /*!
     * \brief createInputFileArguments
     */
    void createInputFileArguments();

    /*!
     * \brief initializeArguments
     * \param message
     * \return
     */
    bool initializeArguments(QString &message) override;

    /*!
     * \brief initializeInputFilesArguments
     * \param message
     * \return
     */
    bool initializeInputFilesArguments(QString &message);

    /*!
     * \brief createGeometriesMap
     */
    void createGeometries();

    /*!
     * \brief createInputs
     */
    void createInputs() override;

    /*!
     * \brief createMainChannelTemperatureInput
     */
    void createMainChannelTemperatureInput();

    /*!
     * \brief createWidthInput
     */
    void createWidthInput();

    void createQHTSInput();

    void createYHTSInput();

    void createAlphaSedInput();

    void createGroundTempInput();

    void createGroundDepthInput();

    /*!
     * \brief createRadiationFluxInput
     */
    void createExternalRadiationFluxInput();

    /*!
     * \brief createExternalHeatFluxInput
     */
    void createExternalHeatFluxInput();

    /*!
     * \brief createOutputs
     */
    void createOutputs() override;

    /*!
     * \brief createMainChannelConductionHeatOutput
     */
    void createMainChannelConductionHeatOutput();

    /*!
     * \brief createMainChannelAdvectionHeatOutput
     */
    void createMainChannelAdvectionHeatOutput();

  private:

    IdBasedArgumentString *m_inputFilesArgument;

    ElementInput *m_mainChannelTempInput,
                 *m_widthInput,
                 *m_QHTSInput,
                 *m_YHTSInput,
                 *m_alphaSediment,
                 *m_groundTemp,
                 *m_groundDepth;

    Dimension *m_timeDimension,
              *m_geometryDimension;

    Unit *m_radiationFluxUnit,
         *m_heatFluxUnit,
         *m_temperatureUnit,
         *m_diffCoeffUnit;

    ElementOutput *m_mainChannelConductionHeat,
                  *m_mainChannelAdvectionHeat;

    ElementHeatSourceInput *m_externalRadiationFluxInput,
                           *m_externalHeatFluxInput;

    std::vector<QSharedPointer<HCGeometry>> m_elementGeometries;
    std::vector<QSharedPointer<HCGeometry>> m_elementJunctionGeometries;
    HTSModel *m_modelInstance;
    HTSComponentInfo *m_HTSComponentInfo;

    HTSComponent *m_parent;
    QList<HydroCouple::ICloneableModelComponent*> m_clones;
};

#endif //HTSCOMPONENT_H
