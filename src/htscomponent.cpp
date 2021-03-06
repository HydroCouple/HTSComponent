/*!
 *  \file    htscomponent.cpp
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
 *  \todo
 *  \warning
 */

#include "stdafx.h"
#include "htscomponent.h"
#include "core/dimension.h"
#include "core/unit.h"
#include "core/unitdimensions.h"
#include "core/abstractoutput.h"
#include "htsmodel.h"
#include "progresschecker.h"
#include "core/valuedefinition.h"
#include "core/idbasedargument.h"
#include "element.h"
#include "elementjunction.h"
#include "spatial/linestring.h"
#include "spatial/point.h"
#include "elementinput.h"
#include "elementoutput.h"
#include "temporal/timedata.h"

using namespace HydroCouple;

HTSComponent::HTSComponent(const QString &id, HTSComponentInfo *modelComponentInfo)
  : AbstractTimeModelComponent(id, modelComponentInfo),
    m_inputFilesArgument(nullptr),
    m_mainChannelTempInput(nullptr),
    m_widthInput(nullptr),
    m_timeDimension(nullptr),
    m_geometryDimension(nullptr),
    m_radiationFluxUnit(nullptr),
    m_heatFluxUnit(nullptr),
    m_temperatureUnit(nullptr),
    m_mainChannelConductionHeat(nullptr),
    m_mainChannelAdvectionHeat(nullptr),
    m_externalRadiationFluxInput(nullptr),
    m_externalHeatFluxInput(nullptr),
    m_modelInstance(nullptr),
    m_HTSComponentInfo(modelComponentInfo),
    m_parent(nullptr)
{
  m_timeDimension = new Dimension("TimeDimension",this);
  m_geometryDimension = new Dimension("ElementGeometryDimension", this);

  m_radiationFluxUnit = new Unit(this);
  m_radiationFluxUnit->setCaption("Radiation Flux (W/m^2)");
  m_radiationFluxUnit->setConversionFactorToSI(1.0);
  m_radiationFluxUnit->setOffsetToSI(0.0);
  m_radiationFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_radiationFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -3.0);

  m_heatFluxUnit = new Unit(this);
  m_heatFluxUnit->setCaption("Heat Source (W or J/s)");
  m_heatFluxUnit->setConversionFactorToSI(1.0);
  m_heatFluxUnit->setOffsetToSI(0.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Length, 2.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -3.0);

  m_temperatureUnit = new Unit(this);
  m_temperatureUnit->setCaption("Temperature (°C)");
  m_temperatureUnit->setConversionFactorToSI(1.0);
  m_temperatureUnit->setOffsetToSI(273.15);
  m_temperatureUnit->dimensionsInternal()->setPower(HydroCouple::Temperature, 1.0);

  m_diffCoeffUnit = new Unit(this);
  m_diffCoeffUnit->setCaption("Diffusion Coefficient (m^2/s)");
  m_diffCoeffUnit->setConversionFactorToSI(1.0);
  m_diffCoeffUnit->setOffsetToSI(0.0);
  m_diffCoeffUnit->dimensionsInternal()->setPower(HydroCouple::Length, 2);
  m_diffCoeffUnit->dimensionsInternal()->setPower(HydroCouple::Time, -1);

  m_soluteFluxUnit = new Unit(this);
  m_soluteFluxUnit->setCaption("Mass Flux (kg/s)");
  m_soluteFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_soluteFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -1.0);

  m_soluteUnit = new Unit(this);
  m_soluteUnit->setCaption("Concentration (kg/m^3)");
  m_soluteUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_soluteUnit->dimensionsInternal()->setPower(HydroCouple::Length, -3.0);

  m_soluteConcQuantity = new Quantity(QVariant::Double, m_soluteUnit, this);
  m_soluteConcFluxQuantity = new Quantity(QVariant::Double, m_soluteFluxUnit, this);

  createArguments();
}

HTSComponent::~HTSComponent()
{
  initializeFailureCleanUp();

  while (m_clones.size())
  {
    HTSComponent *clone =  dynamic_cast<HTSComponent*>(m_clones.first());
    removeClone(clone);
    delete clone;
  }

  if(m_parent)
  {
    m_parent->removeClone(this);
    m_parent = nullptr;
  }
}

QList<QString> HTSComponent::validate()
{
  if(isInitialized())
  {
    setStatus(IModelComponent::Validating,"Validating...");

    //check connections

    setStatus(IModelComponent::Valid,"");
  }
  else
  {
    //throw has not been initialized yet.
  }

  return QList<QString>();
}

void HTSComponent::prepare()
{
  if(!isPrepared() && isInitialized() && m_modelInstance)
  {
    for(auto output :  outputsInternal())
    {
      for(auto adaptedOutput : output->adaptedOutputs())
      {
        adaptedOutput->initialize();
      }
    }

    updateOutputValues(QList<HydroCouple::IOutput*>());

    setStatus(IModelComponent::Updated ,"Finished preparing model");
    setPrepared(true);
  }
  else
  {
    setPrepared(false);
    setStatus(IModelComponent::Failed ,"Error occured when preparing model");
  }
}

void HTSComponent::update(const QList<HydroCouple::IOutput *> &requiredOutputs)
{
  if(status() == IModelComponent::Updated)
  {
    setStatus(IModelComponent::Updating);

    double minConsumerTime = std::max(m_modelInstance->currentDateTime(), getMinimumConsumerTime());

    while (m_modelInstance->currentDateTime() <= minConsumerTime &&
           m_modelInstance->currentDateTime() < m_modelInstance->endDateTime())
    {
      m_modelInstance->update();

      if(progressChecker()->performStep(m_modelInstance->currentDateTime()))
      {
        setStatus(IModelComponent::Updated , "Simulation performed time-step | DateTime: " + QString::number(m_modelInstance->currentDateTime(), 'f') , progressChecker()->progress());
      }
    }

    updateOutputValues(requiredOutputs);

    currentDateTimeInternal()->setJulianDay(m_modelInstance->currentDateTime());

    if(m_modelInstance->currentDateTime() >=  m_modelInstance->endDateTime())
    {
      setStatus(IModelComponent::Done , "Simulation finished successfully", 100);
    }
    else
    {
      if(progressChecker()->performStep(m_modelInstance->currentDateTime()))
      {
        setStatus(IModelComponent::Updated , "Simulation performed time-step | DateTime: " + QString::number(m_modelInstance->currentDateTime(), 'f') , progressChecker()->progress());
      }
      else
      {
        setStatus(IModelComponent::Updated);
      }
    }
  }
}

void HTSComponent::finish()
{
  if(isPrepared())
  {
    setStatus(IModelComponent::Finishing , "STMComponent with id " + id() + " is being disposed" , 100);

    std::list<std::string> errors;
    m_modelInstance->finalize(errors);
    initializeFailureCleanUp();

    setPrepared(false);
    setInitialized(false);

    setStatus(IModelComponent::Finished , "STMComponent with id " + id() + " has been disposed" , 100);
    setStatus(IModelComponent::Created , "STMComponent with id " + id() + " ran successfully and has been re-created" , 100);
  }
}

HTSModel *HTSComponent::modelInstance() const
{
  return m_modelInstance;
}

ICloneableModelComponent *HTSComponent::parent() const
{
  return m_parent;
}

ICloneableModelComponent *HTSComponent::clone()
{
  if(isInitialized())
  {
    HTSComponent *cloneComponent = dynamic_cast<HTSComponent*>(m_HTSComponentInfo->createComponentInstance());
    cloneComponent->setReferenceDirectory(referenceDirectory());

    IdBasedArgumentString *identifierArg = identifierArgument();

    IdBasedArgumentString *cloneIndentifierArg = cloneComponent->identifierArgument();
    (*cloneIndentifierArg)["Id"] = QString((*identifierArg)["Id"]);
    (*cloneIndentifierArg)["Caption"] = QString((*identifierArg)["Caption"]);
    (*cloneIndentifierArg)["Description"] = QString((*identifierArg)["Description"]);

    QString appendName = "_clone_" + QString::number(m_clones.size()) + "_" + QUuid::createUuid().toString().replace("{","").replace("}","");

    //(*cloneComponent->m_inputFilesArgument)["Input File"] = QString((*m_inputFilesArgument)["Input File"]);

    QString inputFilePath = QString((*m_inputFilesArgument)["Input File"]);
    QFileInfo inputFile = getAbsoluteFilePath(inputFilePath);

    if(inputFile.absoluteDir().exists())
    {
      QString suffix = "." + inputFile.completeSuffix();
      inputFilePath = inputFile.absoluteFilePath().replace(suffix,"") + appendName + suffix;
      QFile::copy(inputFile.absoluteFilePath(), inputFilePath);
      (*cloneComponent->m_inputFilesArgument)["Input File"] = inputFilePath;
    }

    QString outputNetCDFFilePath = QString((*m_inputFilesArgument)["Output NetCDF File"]);
    QFileInfo outputNetCDFFile = getAbsoluteFilePath(outputNetCDFFilePath);

    if(!outputNetCDFFilePath.isEmpty() && outputNetCDFFile.absoluteDir().exists())
    {
      QString suffix = "." + outputNetCDFFile.completeSuffix();
      outputNetCDFFilePath = outputNetCDFFile.absoluteFilePath().replace(suffix,"") + appendName + suffix;
      (*cloneComponent->m_inputFilesArgument)["Output NetCDF File"] = outputNetCDFFilePath;
    }

    QString  outputCSVFilePath = QString((*m_inputFilesArgument)["Output CSV File"]);
    QFileInfo outputCSVFile = getAbsoluteFilePath(outputCSVFilePath);

    if(!outputCSVFilePath.isEmpty() && outputCSVFile.absoluteDir().exists())
    {
      QString suffix = "." + outputCSVFile.completeSuffix();
      outputCSVFilePath = outputCSVFile.absoluteFilePath().replace(suffix,"") + appendName + suffix;
      (*cloneComponent->m_inputFilesArgument)["Output CSV File"] = outputCSVFilePath;
    }

    cloneComponent->m_parent = this;
    m_clones.append(cloneComponent);

    emit propertyChanged("Clones");

    cloneComponent->initialize();

    return cloneComponent;
  }


  return nullptr;
}

QList<ICloneableModelComponent*> HTSComponent::clones() const
{
  return m_clones;
}

bool HTSComponent::removeClone(HTSComponent *component)
{
  int removed;

#ifdef USE_OPENMP
#pragma omp critical (HTSComponent)
#endif
  {
    removed = m_clones.removeAll(component);
  }


  if(removed)
  {
    component->m_parent = nullptr;
    emit propertyChanged("Clones");
  }

  return removed;
}

void HTSComponent::initializeFailureCleanUp()
{
  if(m_modelInstance)
  {
    delete m_modelInstance;
    m_modelInstance = nullptr;
  }
}

void HTSComponent::createArguments()
{
  createInputFileArguments();
}

void HTSComponent::createInputFileArguments()
{
  QStringList fidentifiers;
  fidentifiers.append("Input File");
  fidentifiers.append("Output NetCDF File");
  fidentifiers.append("Output CSV File");

  Quantity *fquantity = Quantity::unitLessValues("InputFilesQuantity", QVariant::String, this);
  fquantity->setDefaultValue("");
  fquantity->setMissingValue("");

  Dimension *dimension = new Dimension("IdDimension","Dimension for identifiers",this);

  m_inputFilesArgument = new IdBasedArgumentString("InputFiles", fidentifiers, dimension, fquantity, this);
  m_inputFilesArgument->setCaption("Model Input Files");
  m_inputFilesArgument->addFileFilter("Input File (*.inp)");
  m_inputFilesArgument->setMatchIdentifiersWhenReading(true);

  addArgument(m_inputFilesArgument);
}

bool HTSComponent::initializeArguments(QString &message)
{
  bool initialized = initializeInputFilesArguments(message);

  if(initialized)
  {
    createGeometries();
  }

  return initialized;
}

bool HTSComponent::initializeInputFilesArguments(QString &message)
{

  QString inputFilePath = (*m_inputFilesArgument)["Input File"];
  QFileInfo inputFile = getAbsoluteFilePath(inputFilePath);

  if(inputFile.exists())
  {
    initializeFailureCleanUp();

    m_modelInstance = new HTSModel(this);
    m_modelInstance->setInputFile(inputFile);

    QString netCDFOutput = QString((*m_inputFilesArgument)["Output NetCDF File"]);
    if(!netCDFOutput.isEmpty() && !netCDFOutput.isNull())
      m_modelInstance->setOutputNetCDFFile(QFileInfo(netCDFOutput));

    QString csvOutput = QString((*m_inputFilesArgument)["Output CSV File"]);
    if(!csvOutput.isEmpty() && !csvOutput.isNull())
      m_modelInstance->setOutputCSVFile(QFileInfo(csvOutput));

    std::list<std::string> errors;
    bool initialized = m_modelInstance->initialize(errors);

    for (std::string errorMsg : errors)
    {
      message += "/n" + QString::fromStdString(errorMsg);
    }

    if(initialized)
    {
      timeHorizonInternal()->setJulianDay(m_modelInstance->startDateTime());
      timeHorizonInternal()->setDuration(m_modelInstance->endDateTime() - m_modelInstance->startDateTime());
      currentDateTimeInternal()->setJulianDay(m_modelInstance->startDateTime());
      progressChecker()->reset(m_modelInstance->startDateTime(), m_modelInstance->endDateTime());
    }

    return initialized;
  }
  else
  {
    message = "Input file does not exist: " + inputFile.absoluteFilePath();
    return false;
  }

  return true;
}

void HTSComponent::createGeometries()
{

  m_elementGeometries.clear();
  m_elementJunctionGeometries.clear();

  for(int i = 0; i < m_modelInstance->numElements() ; i++)
  {
    Element *element = m_modelInstance->getElement(i);
    ElementJunction *from = element->upstreamJunction;
    ElementJunction *to   = element->downstreamJunction;

    HCLineString *lineString = new HCLineString(QString::fromStdString(element->id));
    lineString->setMarker(i);
    HCPoint *p1 = new HCPoint(from->x , from->y, QString::fromStdString(from->id), lineString);
    HCPoint *p2 = new HCPoint(to->x , to->y, QString::fromStdString(to->id), lineString);
    lineString->addPoint(p1);
    lineString->addPoint(p2);

    m_elementJunctionGeometries.push_back(QSharedPointer<HCPoint>(new HCPoint(from->x , from->y, from->z, QString::fromStdString(from->id), nullptr)));
    m_elementJunctionGeometries.push_back(QSharedPointer<HCPoint>(new HCPoint(to->x , to->y, to->z, QString::fromStdString(to->id), nullptr)));

    m_elementGeometries.push_back(QSharedPointer<HCLineString>(lineString));
  }
}

void HTSComponent::createInputs()
{
  createMainChannelTemperatureInput();
  createWidthInput();
  createQHTSInput();
  createYHTSInput();
  createAlphaSedInput();
  createGroundTempInput();
  createGroundDepthInput();
  createExternalRadiationFluxInput();
  createExternalHeatFluxInput();

  for(int i = 0 ; i <  m_modelInstance->numSolutes(); i++)
  {
    createMainChannelSoluteConcInput(i);
  }
}

void HTSComponent::createMainChannelTemperatureInput()
{
  Quantity *temperatureQuantity = new Quantity(QVariant::Double, m_temperatureUnit, this);

  m_mainChannelTempInput = new ElementInput("ChannelTemperatureInput",
                                            m_timeDimension,
                                            m_geometryDimension,
                                            temperatureQuantity,
                                            ElementInput::MainChannelTemperature,
                                            this);

  m_mainChannelTempInput->setCaption("Channel Temperature Input (°C)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_mainChannelTempInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_mainChannelTempInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_mainChannelTempInput);

  m_mainChannelTempInput->addTime(dt1);
  m_mainChannelTempInput->addTime(dt2);

  addInput(m_mainChannelTempInput);
}

void HTSComponent::createWidthInput()
{
  Quantity *widthQuantity = Quantity::lengthInMeters(this);

  m_widthInput = new ElementInput("WidthInput",
                                  m_timeDimension,
                                  m_geometryDimension,
                                  widthQuantity,
                                  ElementInput::Width,
                                  this);

  m_widthInput->setCaption("Element Width (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_widthInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_widthInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_widthInput);

  m_widthInput->addTime(dt1);
  m_widthInput->addTime(dt2);

  addInput(m_widthInput);
}

void HTSComponent::createQHTSInput()
{
  Quantity *flowQuantity = Quantity::flowInCMS(this);

  m_QHTSInput = new ElementInput("AdvectiveTransportCoeffInput",
                                 m_timeDimension,
                                 m_geometryDimension,
                                 flowQuantity,
                                 ElementInput::QHTS,
                                 this);

  m_QHTSInput->setCaption("Advective Transport Coefficient (m^3/s)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_QHTSInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_QHTSInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_QHTSInput);

  m_QHTSInput->addTime(dt1);
  m_QHTSInput->addTime(dt2);

  addInput(m_QHTSInput);
}

void HTSComponent::createYHTSInput()
{
  Quantity *depthQuantity = Quantity::lengthInMeters(this);

  m_YHTSInput = new ElementInput("HTSDepthInput",
                                 m_timeDimension,
                                 m_geometryDimension,
                                 depthQuantity,
                                 ElementInput::YHTS,
                                 this);

  m_YHTSInput->setCaption("HTS Depth (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_YHTSInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_YHTSInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_YHTSInput);

  m_YHTSInput->addTime(dt1);
  m_YHTSInput->addTime(dt2);

  addInput(m_YHTSInput);
}

void HTSComponent::createAlphaSedInput()
{
  Quantity *diffCoeffQuantity = new Quantity(QVariant::Double, m_diffCoeffUnit, this);

  m_alphaSediment = new ElementInput("SedThermaDiffusivityCoeffInput",
                                     m_timeDimension,
                                     m_geometryDimension,
                                     diffCoeffQuantity,
                                     ElementInput::AlphaSed,
                                     this);

  m_alphaSediment->setCaption("Sediment Thermal Diffusivity Coefficient (m^2/s)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_alphaSediment->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_alphaSediment);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_alphaSediment);

  m_alphaSediment->addTime(dt1);
  m_alphaSediment->addTime(dt2);

  addInput(m_alphaSediment);
}

void HTSComponent::createGroundTempInput()
{
  Quantity *temperatureQuantity = new Quantity(QVariant::Double, m_temperatureUnit, this);

  m_groundTemp = new ElementInput("GroundTemperatureInput",
                                  m_timeDimension,
                                  m_geometryDimension,
                                  temperatureQuantity,
                                  ElementInput::GroundTemperature,
                                  this);

  m_groundTemp->setCaption("Ground Temperature Input (°C)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_groundTemp->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_groundTemp);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_groundTemp);

  m_groundTemp->addTime(dt1);
  m_groundTemp->addTime(dt2);

  addInput(m_groundTemp);
}

void HTSComponent::createGroundDepthInput()
{
  Quantity *depthQuantity = Quantity::lengthInMeters(this);

  m_groundDepth = new ElementInput("GroundDepthInput",
                                   m_timeDimension,
                                   m_geometryDimension,
                                   depthQuantity,
                                   ElementInput::GroundConductionDepth,
                                   this);

  m_groundDepth->setCaption("Ground Depth (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_groundDepth->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_groundDepth);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_groundDepth);

  m_groundDepth->addTime(dt1);
  m_groundDepth->addTime(dt2);

  addInput(m_groundDepth);
}

void HTSComponent::createExternalRadiationFluxInput()
{
  Quantity *radiationQuantity = new Quantity(QVariant::Double, m_radiationFluxUnit, this);

  m_externalRadiationFluxInput = new ElementHeatSourceInput("RadiationFluxInput",
                                                            m_timeDimension,
                                                            m_geometryDimension,
                                                            radiationQuantity,
                                                            ElementHeatSourceInput::RadiativeFlux,
                                                            this);
  m_externalRadiationFluxInput->setCaption("External Radiation Flux (W/m^2)");


  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_externalRadiationFluxInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_externalRadiationFluxInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_externalRadiationFluxInput);

  m_externalRadiationFluxInput->addTime(dt1);
  m_externalRadiationFluxInput->addTime(dt2);

  addInput(m_externalRadiationFluxInput);
}

void HTSComponent::createExternalHeatFluxInput()
{
  Quantity *heatFluxQuantity = new Quantity(QVariant::Double, m_heatFluxUnit, this);

  m_externalHeatFluxInput = new ElementHeatSourceInput("HeatFluxInput",
                                                       m_timeDimension,
                                                       m_geometryDimension,
                                                       heatFluxQuantity,
                                                       ElementHeatSourceInput::HeatFlux,
                                                       this);

  m_externalHeatFluxInput->setCaption("External Heat Flux (J/s)");


  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_externalHeatFluxInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_externalHeatFluxInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_externalHeatFluxInput);

  m_externalHeatFluxInput->addTime(dt1);
  m_externalHeatFluxInput->addTime(dt2);

  addInput(m_externalHeatFluxInput);
}

void HTSComponent::createMainChannelSoluteConcInput(int soluteIndex)
{
  QString soluteName = QString::fromStdString(m_modelInstance->solute(soluteIndex));

  ElementInput *soluteConcInput  = new ElementInput(soluteName + "Input",
                                                    m_timeDimension,
                                                    m_geometryDimension,
                                                    m_soluteConcQuantity,
                                                    ElementInput::MainChannelSoluteConc,
                                                    this);

  soluteConcInput->setCaption("Element " + soluteName + " Concentration (kg/m^3)");
  soluteConcInput->setSoluteIndex(soluteIndex);

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  soluteConcInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteConcInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteConcInput);

  soluteConcInput->addTime(dt1);
  soluteConcInput->addTime(dt2);

  addInput(soluteConcInput);
}

void HTSComponent::createOutputs()
{
  createMainChannelConductionHeatOutput();
  createMainChannelAdvectionHeatOutput();

  for(int i = 0 ; i <  m_modelInstance->numSolutes(); i++)
  {
    createMainChannelSoluteDiffFluxOutput(i);
    createMainChannelSoluteAdvFluxOutput(i);
  }
}

void HTSComponent::createMainChannelConductionHeatOutput()
{
  Quantity *heatFluxQuantity = new Quantity(QVariant::Double, m_heatFluxUnit, this);

  m_mainChannelConductionHeat = new ElementOutput("MainChannelConductionHeatFluxOutput",
                                                  m_timeDimension,
                                                  m_geometryDimension,
                                                  heatFluxQuantity,
                                                  ElementOutput::ChannelConductionHeat,
                                                  this);

  m_mainChannelConductionHeat->setCaption("Main Channel Conduction Heat Flux (J/s)");


  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_mainChannelConductionHeat->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_mainChannelConductionHeat);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_mainChannelConductionHeat);

  m_mainChannelConductionHeat->addTime(dt1);
  m_mainChannelConductionHeat->addTime(dt2);

  addOutput(m_mainChannelConductionHeat);
}

void HTSComponent::createMainChannelAdvectionHeatOutput()
{
  Quantity *heatFluxQuantity = new Quantity(QVariant::Double, m_heatFluxUnit, this);

  m_mainChannelAdvectionHeat = new ElementOutput("MainChannelAdvectionHeatFluxOutput",
                                                 m_timeDimension,
                                                 m_geometryDimension,
                                                 m_soluteConcFluxQuantity,
                                                 ElementOutput::ChannelAdvectionHeat,
                                                 this);

  m_mainChannelAdvectionHeat->setCaption("Main Channel Advection Heat Flux (J/s)");


  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_mainChannelAdvectionHeat->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_mainChannelAdvectionHeat);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_mainChannelAdvectionHeat);

  m_mainChannelAdvectionHeat->addTime(dt1);
  m_mainChannelAdvectionHeat->addTime(dt2);

  addOutput(m_mainChannelAdvectionHeat);
}

void HTSComponent::createMainChannelSoluteDiffFluxOutput(int soluteIndex)
{
  QString soluteName = QString::fromStdString(m_modelInstance->solute(soluteIndex));

  ElementOutput *soluteFluxOutput  = new ElementOutput(soluteName + "DiffOutput",
                                                            m_timeDimension,
                                                            m_geometryDimension,
                                                            m_soluteConcFluxQuantity,
                                                            ElementOutput::ChannelSoluteDiffusionFlux,
                                                            this);
  soluteFluxOutput->setCaption("Element " + soluteName + " Diffusive Flux (kg/s)");
  soluteFluxOutput->setSoluteIndex(soluteIndex);

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  soluteFluxOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteFluxOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteFluxOutput);

  soluteFluxOutput->addTime(dt1);
  soluteFluxOutput->addTime(dt2);

  addOutput(soluteFluxOutput);
}

void HTSComponent::createMainChannelSoluteAdvFluxOutput(int soluteIndex)
{
  QString soluteName = QString::fromStdString(m_modelInstance->solute(soluteIndex));

  ElementOutput *soluteFluxOutput  = new ElementOutput(soluteName + "AdvOutput",
                                                            m_timeDimension,
                                                            m_geometryDimension,
                                                            m_soluteConcFluxQuantity,
                                                            ElementOutput::ChannelSoluteAdvectionFlux,
                                                            this);

  soluteFluxOutput->setCaption("Element " + soluteName + " Advective Flux (kg/s)");
  soluteFluxOutput->setSoluteIndex(soluteIndex);

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  soluteFluxOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteFluxOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteFluxOutput);

  soluteFluxOutput->addTime(dt1);
  soluteFluxOutput->addTime(dt2);

  addOutput(soluteFluxOutput);
}
