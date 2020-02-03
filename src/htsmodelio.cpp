#include "stdafx.h"
#include "htsmodel.h"
#include "element.h"
#include "elementjunction.h"
#include "sourcebc.h"
#include "hydraulicsbc.h"
#include "radiativefluxbc.h"
#include "elementbc.h"
#include "temporal/timedata.h"
#include "threadsafenetcdf/threadsafencfile.h"
#include "threadsafenetcdf/threadsafencdim.h"
#include "threadsafenetcdf/threadsafencvar.h"
#include "threadsafenetcdf/threadsafencatt.h"
#include "temporal/timeseries.h"

#include <QDir>
#include <QDate>
#include <cstdlib>
#include <errno.h>

#ifdef USE_NETCDF

using namespace netCDF;
using namespace netCDF::exceptions;
using namespace SDKTemporal;

#endif

using namespace std;

bool HTSModel::verbose() const
{
  return m_verbose;
}

void HTSModel::setVerbose(bool verbose)
{
  m_verbose = verbose;
}

int HTSModel::printFrequency() const
{
  return m_printFrequency;
}

void HTSModel::setPrintFrequency(int printFreq)
{
  m_printFrequency = printFreq;
}

int HTSModel::flushToDiskFrequency() const
{
  return m_flushToDiskFrequency;
}

void HTSModel::setFlushToDiskFrequency(int diskFlushFrequency)
{
  m_flushToDiskFrequency = diskFlushFrequency;
}

QFileInfo HTSModel::inputFile() const
{
  return m_inputFile;
}

void HTSModel::setInputFile(const QFileInfo &inputFile)
{
  m_inputFile = inputFile;
}

QFileInfo HTSModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void HTSModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

QFileInfo HTSModel::outputNetCDFFile() const
{
  return m_outputNetCDFFileInfo;
}

void HTSModel::setOutputNetCDFFile(const QFileInfo &outputNetCDFFile)
{
  m_outputNetCDFFileInfo = outputNetCDFFile;
}

void HTSModel::printStatus()
{
  m_currentPrintCount++;

  if (m_currentPrintCount >= m_printFrequency)
  {

    printf("HTS TimeStep (s): %f\tDateTime: %f\tTemp (°C) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalHeatBalance: %g (KJ)}", m_timeStep, m_currentDateTime,
           m_heatSolver->getIterations(), m_heatSolver->maxIterations(), m_minTemp, m_maxTemp, m_totalHeatBalance);

    for (size_t j = 0; j < m_solutes.size(); j++)
    {
      std::string &solute = m_solutes[j];
      ODESolver *solver = m_soluteSolvers[j];
      printf("\t%s (kg/m) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalMassBalance: %g (kg)}", solute.c_str(), solver->getIterations(), solver->maxIterations(),
             m_minSolute[j], m_maxSolute[j], m_totalSoluteMassBalance[j]);
    }

    printf("\n");

    m_currentPrintCount = 0;
  }
}

bool HTSModel::initializeInputFiles(list<string> &errors)
{

  if(QFile::exists(m_inputFile.absoluteFilePath()))
  {
    QFile file(m_inputFile.absoluteFilePath());

    if(file.open(QIODevice::ReadOnly))
    {
      m_outNetCDFVariablesOnOff.clear();
      m_timeSeries.clear();

      m_delimiters = QRegExp("(\\,|\\t|\\;|\\s+)");
      int currentFlag = -1;
      m_addedSoluteCount = 0;

      QTextStream streamReader(&file);
      int lineCount = 0;
      while (!streamReader.atEnd())
      {
        QString line = streamReader.readLine().trimmed();
        lineCount++;

        if (!line.isEmpty() && !line.isNull())
        {
          bool readSuccess = true;
          QString error = "";

          auto it = m_inputFileFlags.find(line.toStdString());

          if (it != m_inputFileFlags.cend())
          {
            currentFlag = it->second;
          }
          else if (!QStringRef::compare(QStringRef(&line, 0, 2), ";;"))
          {
            //commment do nothing
          }
          else
          {
            switch (currentFlag)
            {
              case 1:
                readSuccess = readInputFileOptionTag(line, error);
                break;
              case 2:
                readSuccess = readInputFileOutputTag(line, error);
                break;
              case 3:
                readSuccess = readInputFileSolutesTag(line, error);
                break;
              case 4:
                readSuccess = readInputFileElementJunctionsTag(line, error);
                break;
              case 5:
                readSuccess = readInputFileElementsTag(line, error);
                break;
              case 6:
                readSuccess = readInputFileSourcesTag(line, error);
                break;
              case 7:
                readSuccess = readInputFileHydraulicsTag(line, error);
                break;
              case 8:
                readSuccess = readInputFileRadiativeFluxesTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileElementBCTag(line, error);
                break;
              case 10:
                readSuccess = readInputFileTimeSeriesTag(line, error);
                break;
              case 11:
                readSuccess = readOutputVariableOnOff(line, error);
                break;
            }
          }

          if (!readSuccess)
          {
            errors.push_back("Line " + std::to_string(lineCount) + " : " + error.toStdString());
            file.close();
            return false;
          }
        }
      }

      file.close();
    }
  }

  return true;
}

bool HTSModel::initializeOutputFiles(list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool HTSModel::initializeCSVOutputFile(list<string> &errors)
{

  if (m_outputCSVFileInfo.isRelative())
  {
    m_outputCSVFileInfo = relativePathToAbsolute(m_outputCSVFileInfo);
  }

  QString file = m_outputCSVFileInfo.absoluteFilePath();

  if (!file.isEmpty() && !file.isNull() && !m_outputCSVFileInfo.absoluteDir().exists())
  {
    errors.push_back("Output shapefile directory does not exist: " + file.toStdString());
    return false;
  }

  if (!file.isEmpty() && !file.isNull())
  {
    if(m_outputCSVFileInfo.isDir())
      return true;

    if (m_outputCSVStream.device() == nullptr)
    {
      QFile *device = new QFile(file, this);
      m_outputCSVStream.setDevice(device);
    }

    if (m_outputCSVStream.device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_outputCSVStream.setRealNumberPrecision(10);
      m_outputCSVStream.setRealNumberNotation(QTextStream::SmartNotation);
      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Depth, Width, XSectionArea, Dispersion, "
                           "Temperature, TotalExternalHeatFluxesBalance, TotalRadFluxesHeatBalance, "
                           "TotalHeatBalance";

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        QString soluteName = QString::fromStdString(m_solutes[i]);
        m_outputCSVStream << ", " << soluteName << ", "
                          << "TotalAdvDispMassBalance_" + soluteName << ", "
                          << "TotalExternalFluxMassBalance_" + soluteName << ","
                          << "TotalMassBalance_" + soluteName ;
      }

      m_outputCSVStream << endl;
      m_outputCSVStream.flush();
    }


    return true;
  }

  return false;
}

bool HTSModel::initializeNetCDFOutputFile(list<string> &errors)
{
#ifdef USE_NETCDF

  if (m_outputNetCDFFileInfo.isRelative())
  {
    m_outputNetCDFFileInfo = relativePathToAbsolute(m_outputNetCDFFileInfo);
  }

  if (m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() || m_outputNetCDFFileInfo.isDir())
  {
    return true;
  }
  else if (!m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() &&
           !m_outputNetCDFFileInfo.absoluteFilePath().isNull() &&
           !m_outputNetCDFFileInfo.absoluteDir().exists())
  {
    std::string message = "NetCDF output file directory does not exist: " + m_outputNetCDFFileInfo.absoluteFilePath().toStdString();
    errors.push_back(message);
    return false;
  }


  bool returnValue = false;


  closeOutputNetCDFFile();

  try
  {

    m_outNetCDFVariables.clear();

    m_outputNetCDF = new ThreadSafeNcFile(m_outputNetCDFFileInfo.absoluteFilePath().toStdString(), NcFile::replace);

    //time variable
    ThreadSafeNcDim timeDim =  m_outputNetCDF->addDim("time");
    ThreadSafeNcVar timeVar =  m_outputNetCDF->addVar("time", NcType::nc_DOUBLE, timeDim);
    timeVar.putAtt("long_name", "Time");
    timeVar.putAtt("standard_name", "time");
    timeVar.putAtt("calendar", "julian");
    m_outNetCDFVariables["time"] = timeVar;

    //Add Solutes
    ThreadSafeNcDim solutesDim =  m_outputNetCDF->addDim("solutes", m_solutes.size());
    ThreadSafeNcVar solutes =  m_outputNetCDF->addVar("solute_names", NcType::nc_STRING, solutesDim);
    solutes.putAtt("long_name", "Solutes");
    m_outNetCDFVariables["solute_names"] = solutes;

    if (m_solutes.size())
    {
      char **soluteNames = new char *[m_solutes.size()];

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        string soluteName = m_solutes[i];
        soluteNames[i] = new char[soluteName.size() + 1];
        strcpy(soluteNames[i], soluteName.c_str());
      }

      solutes.putVar(soluteNames);

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        delete[] soluteNames[i];
      }

      delete[] soluteNames;
    }

    //Add element junctions
    ThreadSafeNcDim junctionDim =  m_outputNetCDF->addDim("element_junctions", m_elementJunctions.size());

    ThreadSafeNcVar junctionIdentifiers =  m_outputNetCDF->addVar("element_junction_id", NcType::nc_STRING, junctionDim);
    junctionIdentifiers.putAtt("long_name", "Element Junction Identifier");
    m_outNetCDFVariables["element_junction_id"] = junctionIdentifiers;

    ThreadSafeNcVar junctionX =  m_outputNetCDF->addVar("x", NcType::nc_FLOAT, junctionDim);
    junctionX.putAtt("long_name", "Junction X-Coordinate");
    junctionX.putAtt("units", "m");
    m_outNetCDFVariables["x"] = junctionX;

    ThreadSafeNcVar junctionY =  m_outputNetCDF->addVar("y", NcType::nc_FLOAT, junctionDim);
    junctionY.putAtt("long_name", "Junction Y-coordinate");
    junctionY.putAtt("units", "m");
    m_outNetCDFVariables["y"] = junctionY;

    ThreadSafeNcVar junctionZ =  m_outputNetCDF->addVar("z", NcType::nc_FLOAT, junctionDim);
    junctionZ.putAtt("long_name", "Junction Z-Coordinate");
    junctionZ.putAtt("units", "m");
    m_outNetCDFVariables["z"] = junctionZ;

    float *vertx = new float[m_elementJunctions.size()];
    float *verty = new float[m_elementJunctions.size()];
    float *vertz = new float[m_elementJunctions.size()];
    char **junctionIds = new char *[m_elementJunctions.size()];

    //write other relevant junction attributes here.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elementJunctions.size(); i++)
    {
      ElementJunction *junction = m_elementJunctions[i];

      junctionIds[i] = new char[junction->id.size() + 1];
      strcpy(junctionIds[i], junction->id.c_str());

      vertx[i] = junction->x;
      verty[i] = junction->y;
      vertz[i] = junction->z;
    }

    junctionX.putVar(vertx);
    junctionY.putVar(verty);
    junctionZ.putVar(vertz);
    junctionIdentifiers.putVar(junctionIds);

    delete[] vertx;
    delete[] verty;
    delete[] vertz;

    for (size_t i = 0; i < m_elementJunctions.size(); i++)
    {
      delete[] junctionIds[i];
    }

    delete[] junctionIds;

    //Add Elements
    ThreadSafeNcDim elementsDim =  m_outputNetCDF->addDim("elements", m_elements.size());

    ThreadSafeNcVar elementIdentifiers =  m_outputNetCDF->addVar("element_id", NcType::nc_STRING, elementsDim);
    elementIdentifiers.putAtt("long_name", "Element Identifier");
    m_outNetCDFVariables["element_id"] = elementIdentifiers;

    ThreadSafeNcVar elementFromJunction =  m_outputNetCDF->addVar("from_junction", NcType::nc_INT64, elementsDim);
    elementFromJunction.putAtt("long_name", "Upstream Junction");
    m_outNetCDFVariables["from_junction"] = elementFromJunction;

    ThreadSafeNcVar elementToJunction =  m_outputNetCDF->addVar("to_junction", NcType::nc_INT64, elementsDim);
    elementToJunction.putAtt("long_name", "Downstream Junction");
    m_outNetCDFVariables["to_junction"] = elementToJunction;


    //    ThreadSafeNcVar elementsVar =  m_outputNetCDF->addVar("elements", NcType::nc_FLOAT, elementsDim);
    //    elementsVar.putAtt("long_name", "Distance");
    //    elementsVar.putAtt("units", "m");
    //    m_outNetCDFVariables["elements"] = elementsVar;


    int *fromJunctions = new int[m_elements.size()];
    int *toJunctions = new int[m_elements.size()];
    char **elementIds = new char *[m_elements.size()];
    //float *els = new float[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int) m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      elementIds[i] = new char[element->id.size() + 1];
      strcpy(elementIds[i], element->id.c_str());

      fromJunctions[i] = element->upstreamJunction->index;
      toJunctions[i] = element->downstreamJunction->index;
      //els[i] = element->distanceFromUpStreamJunction;
    }

    elementIdentifiers.putVar(elementIds);
    elementFromJunction.putVar(fromJunctions);
    elementToJunction.putVar(toJunctions);
    //elementsVar.putVar(els);

    delete[] fromJunctions;
    delete[] toJunctions;
    //delete[] els;

    for (size_t i = 0; i < m_elements.size(); i++)
    {
      delete[] elementIds[i];
    }

    delete[] elementIds;

    auto varOnOff = [this](const std::string& name) -> bool
    {
      return m_outNetCDFVariablesOnOff.find(name) != m_outNetCDFVariablesOnOff.end() ? m_outNetCDFVariablesOnOff[name] : true;
    };


    if((m_outNetCDFVariablesOnOff["depth"] = varOnOff("depth"))){

      ThreadSafeNcVar depthVar =  m_outputNetCDF->addVar("depth", "float",
                                                         std::vector<std::string>({"time", "elements"}));
      depthVar.putAtt("long_name", "Flow Depth");
      depthVar.putAtt("units", "m");
      m_outNetCDFVariables["depth"] = depthVar;
      m_outNetCDFVariablesIOFunctions["depth"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *depth = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          depth[i] = static_cast<float>(elements[i]->depth);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), depth);
        delete[] depth;
      };
    }

    if((m_outNetCDFVariablesOnOff["width"] = varOnOff("width"))){

      ThreadSafeNcVar widthVar =  m_outputNetCDF->addVar("width", "float",
                                                         std::vector<std::string>({"time", "elements"}));
      widthVar.putAtt("long_name", "Flow Top Width");
      widthVar.putAtt("units", "m");
      m_outNetCDFVariables["width"] = widthVar;
      m_outNetCDFVariablesIOFunctions["width"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *width = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          width[i] = static_cast<float>(elements[i]->depth);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), width);
        delete[] width;
      };
    }

    if((m_outNetCDFVariablesOnOff["xsection_area"] = varOnOff("xsection_area"))){

      ThreadSafeNcVar xsectAreaVar =  m_outputNetCDF->addVar("xsection_area", "float",
                                                             std::vector<std::string>({"time", "elements"}));
      xsectAreaVar.putAtt("long_name", "Flow Cross-Sectional Area");
      xsectAreaVar.putAtt("units", "m^2");
      m_outNetCDFVariables["xsection_area"] = xsectAreaVar;
      m_outNetCDFVariablesIOFunctions["xsection_area"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *xsection_area = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          xsection_area[i] = static_cast<float>(elements[i]->xSectionArea);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), xsection_area);
        delete[] xsection_area;
      };
    }

    if((m_outNetCDFVariablesOnOff["hts_temperature"] = varOnOff("hts_temperature"))){

      ThreadSafeNcVar temperatureVar =  m_outputNetCDF->addVar("hts_temperature", "float",
                                                               std::vector<std::string>({"time", "elements"}));
      temperatureVar.putAtt("long_name", "HTS Temperature");
      temperatureVar.putAtt("units", "°C");
      m_outNetCDFVariables["hts_temperature"] = temperatureVar;
      m_outNetCDFVariablesIOFunctions["hts_temperature"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *hts_temperature = new float[elements.size()];

        for (size_t i = 0; i < elements.size(); i++)
        {
          hts_temperature[i] = static_cast<float>(elements[i]->temperature.value);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), hts_temperature);
        delete[] hts_temperature;
      };
    }

    if((m_outNetCDFVariablesOnOff["channel_conduction_heat_flux"] = varOnOff("channel_conduction_heat_flux"))){

      ThreadSafeNcVar channelConductionFluxVar =  m_outputNetCDF->addVar("channel_conduction_heat_flux", "float",
                                                                         std::vector<std::string>({"time", "elements"}));
      channelConductionFluxVar.putAtt("long_name", "Channel Conduction Heat Flux");
      channelConductionFluxVar.putAtt("units", "W/m2");
      m_outNetCDFVariables["channel_conduction_heat_flux"] = channelConductionFluxVar;
      m_outNetCDFVariablesIOFunctions["channel_conduction_heat_flux"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *channel_conduction_heat_flux = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          channel_conduction_heat_flux[i] = static_cast<float>(element->mainChannelConductionHeat / (element->width * element->length));
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), channel_conduction_heat_flux);
        delete[] channel_conduction_heat_flux;
      };
    }

    if((m_outNetCDFVariablesOnOff["ground_conduction_heat_flux"] = varOnOff("ground_conduction_heat_flux"))){

      ThreadSafeNcVar groundConductionFluxVar =  m_outputNetCDF->addVar("ground_conduction_heat_flux", "float",
                                                                        std::vector<std::string>({"time", "elements"}));
      groundConductionFluxVar.putAtt("long_name", "Ground Conduction Heat Flux");
      groundConductionFluxVar.putAtt("units", "W/m2");
      m_outNetCDFVariables["ground_conduction_heat_flux"] = groundConductionFluxVar;
      m_outNetCDFVariablesIOFunctions["ground_conduction_heat_flux"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *ground_conduction_heat_flux = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          ground_conduction_heat_flux[i] = static_cast<float>(element->groundConductionHeat / (element->width * element->length));
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), ground_conduction_heat_flux);
        delete[] ground_conduction_heat_flux;
      };
    }

    if((m_outNetCDFVariablesOnOff["channel_advection_heat_flux"] = varOnOff("channel_advection_heat_flux"))){

      ThreadSafeNcVar channelAdvectionHeatFluxVar =  m_outputNetCDF->addVar("channel_advection_heat_flux", "float",
                                                                            std::vector<std::string>({"time", "elements"}));
      channelAdvectionHeatFluxVar.putAtt("long_name", "Channel Advection Heat Flux");
      channelAdvectionHeatFluxVar.putAtt("units", "W/m2");
      m_outNetCDFVariables["channel_advection_heat_flux"] = channelAdvectionHeatFluxVar;
      m_outNetCDFVariablesIOFunctions["channel_advection_heat_flux"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *channel_advection_heat_flux = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          channel_advection_heat_flux[i] = static_cast<float>(element->mainChannelAdvectionHeat / (element->width * element->length));
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), channel_advection_heat_flux);
        delete[] channel_advection_heat_flux;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_heat_balance"] = varOnOff("total_heat_balance"))){

      ThreadSafeNcVar totalHeatBalanceVar =  m_outputNetCDF->addVar("total_heat_balance", "float",
                                                                    std::vector<std::string>({"time"}));
      totalHeatBalanceVar.putAtt("long_name", "Total Heat Balance");
      totalHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_heat_balance"] = totalHeatBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["total_radiation_flux_heat_balance"] = varOnOff("total_radiation_flux_heat_balance"))){

      ThreadSafeNcVar totalRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_radiation_flux_heat_balance", "float",
                                                                                 std::vector<std::string>({"time"}));
      totalRadiationFluxHeatBalanceVar.putAtt("long_name", "Total Radiation Flux Heat Balance");
      totalRadiationFluxHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_radiation_flux_heat_balance"] = totalRadiationFluxHeatBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["total_external_heat_flux_balance"] = varOnOff("total_external_heat_flux_balance"))){

      ThreadSafeNcVar totalExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_external_heat_flux_balance", "float",
                                                                                std::vector<std::string>({"time"}));
      totalExternalHeatFluxBalanceVar.putAtt("long_name", "Total External Heat Flux Balance");
      totalExternalHeatFluxBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_external_heat_flux_balance"] = totalExternalHeatFluxBalanceVar;
    }

    if((m_outNetCDFVariablesOnOff["total_element_heat_balance"] = varOnOff("total_element_heat_balance"))){

      ThreadSafeNcVar totalElementHeatBalanceVar =  m_outputNetCDF->addVar("total_element_heat_balance", "float",
                                                                           std::vector<std::string>({"time", "elements"}));
      totalElementHeatBalanceVar.putAtt("long_name", "Total Heat Balance");
      totalElementHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_heat_balance"] = totalElementHeatBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_heat_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_heat_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_heat_balance[i] = static_cast<float>(element->totalHeatBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_heat_balance);
        delete[] total_element_heat_balance;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_element_radiation_flux_heat_balance"] = varOnOff("total_element_radiation_flux_heat_balance"))){

      ThreadSafeNcVar totalElementRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_element_radiation_flux_heat_balance", "float",
                                                                                        std::vector<std::string>({"time", "elements"}));
      totalElementRadiationFluxHeatBalanceVar.putAtt("long_name", "Total Element Radiation Flux Heat Balance");
      totalElementRadiationFluxHeatBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_radiation_flux_heat_balance"] = totalElementRadiationFluxHeatBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_radiation_flux_heat_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_radiation_flux_heat_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_radiation_flux_heat_balance[i] = static_cast<float>(element->totalRadiationFluxesHeatBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_radiation_flux_heat_balance);
        delete[] total_element_radiation_flux_heat_balance;
      };
    }

    if((m_outNetCDFVariablesOnOff["total_element_external_heat_flux_balance"] = varOnOff("total_element_external_heat_flux_balance"))){

      ThreadSafeNcVar totalElementExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_element_external_heat_flux_balance", "float",
                                                                                       std::vector<std::string>({"time", "elements"}));
      totalElementExternalHeatFluxBalanceVar.putAtt("long_name", "Total Element External Heat Flux Balance");
      totalElementExternalHeatFluxBalanceVar.putAtt("units", "KJ");
      m_outNetCDFVariables["total_element_external_heat_flux_balance"] = totalElementExternalHeatFluxBalanceVar;
      m_outNetCDFVariablesIOFunctions["total_element_external_heat_flux_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
      {
        float *total_element_external_heat_flux_balance = new float[elements.size()];
        for (size_t i = 0; i < elements.size(); i++)
        {
          Element *element = elements[i];
          total_element_external_heat_flux_balance[i] = static_cast<float>(element->totalExternalHeatFluxesBalance);
        }
        variable.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, elements.size()}), total_element_external_heat_flux_balance);
        delete[] total_element_external_heat_flux_balance;
      };
    }

    if(m_solutes.size())
    {
      if((m_outNetCDFVariablesOnOff["solute_concentration"] = varOnOff("solute_concentration"))){

        ThreadSafeNcVar solutesVar =  m_outputNetCDF->addVar("solute_concentration", "float",
                                                             std::vector<std::string>({"time", "solutes", "elements"}));
        solutesVar.putAtt("long_name", "Solute Concentration");
        solutesVar.putAtt("units", "kg/m^3");
        m_outNetCDFVariables["solute_concentration"] = solutesVar;
        m_outNetCDFVariablesIOFunctions["solute_concentration"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
        {
          int numSolutes = elements[0]->numSolutes;
          float *solute_concentration = new float[elements.size() * numSolutes];

          for (size_t i = 0; i < elements.size(); i++)
          {
            Element *element = elements[i];

            for(int j = 0 ; j < numSolutes; j++)
            {
              solute_concentration[i + j * elements.size()] = static_cast<float>(element->soluteConcs[j].value);
            }
          }
          variable.putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)numSolutes, elements.size()}), solute_concentration);
          delete[] solute_concentration;
        };
      }

      if((m_outNetCDFVariablesOnOff["total_element_solute_mass_balance"] = varOnOff("total_element_solute_mass_balance"))){

        ThreadSafeNcVar totalElementSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_element_solute_mass_balance", "float",
                                                                                   std::vector<std::string>({"time", "solutes", "elements"}));
        totalElementSoluteMassBalanceVar.putAtt("long_name", "Total Element Solute Mass Balance");
        totalElementSoluteMassBalanceVar.putAtt("units", "kg");
        m_outNetCDFVariables["total_element_solute_mass_balance"] = totalElementSoluteMassBalanceVar;
        m_outNetCDFVariablesIOFunctions["total_element_solute_mass_balance"] = [](size_t currentTime, ThreadSafeNcVar &variable, const std::vector<Element*>& elements)
        {
          int numSolutes = elements[0]->numSolutes;
          float *total_element_solute_mass_balance = new float[elements.size() * numSolutes];

          for (size_t i = 0; i < elements.size(); i++)
          {
            Element *element = elements[i];

            for(int j = 0 ; j < numSolutes; j++)
            {
              total_element_solute_mass_balance[i + j * elements.size()] = static_cast<float>(element->totalSoluteMassBalance[j]);
            }
          }
          variable.putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)numSolutes, elements.size()}), total_element_solute_mass_balance);
          delete[] total_element_solute_mass_balance;
        };
      }

      if((m_outNetCDFVariablesOnOff["total_solute_mass_balance"] = varOnOff("total_solute_mass_balance"))){

        ThreadSafeNcVar totalSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_solute_mass_balance", "float",
                                                                            std::vector<std::string>({"time", "solutes"}));
        totalSoluteMassBalanceVar.putAtt("long_name", "Total Solute Mass Balance");
        totalSoluteMassBalanceVar.putAtt("units", "kg");
        m_outNetCDFVariables["total_solute_mass_balance"] = totalSoluteMassBalanceVar;
      }

      if((m_outNetCDFVariablesOnOff["total_external_solute_flux_mass_balance"] = varOnOff("total_external_solute_flux_mass_balance"))){

        ThreadSafeNcVar totalExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->addVar("total_external_solute_flux_mass_balance", "float",
                                                                                        std::vector<std::string>({"time", "solutes"}));
        totalExternalSoluteFluxMassBalanceVar.putAtt("long_name", "Total External Solute Flux Mass Balance");
        totalExternalSoluteFluxMassBalanceVar.putAtt("units", "kg");
        m_outNetCDFVariables["total_external_solute_flux_mass_balance"] = totalExternalSoluteFluxMassBalanceVar;
      }

      if((m_outNetCDFVariablesOnOff["total_element_external_solute_flux_mass_balance"] = varOnOff("total_element_external_solute_flux_mass_balance"))){

        ThreadSafeNcVar totalElementExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->addVar("total_element_external_solute_flux_mass_balance", "float",
                                                                                               std::vector<std::string>({"time", "solutes"}));
        totalElementExternalSoluteFluxMassBalanceVar.putAtt("long_name", "Total External Solute Flux Mass Balance");
        totalElementExternalSoluteFluxMassBalanceVar.putAtt("units", "kg");
        m_outNetCDFVariables["total_element_external_solute_flux_mass_balance"] = totalElementExternalSoluteFluxMassBalanceVar;
      }
    }


    m_optionalOutputVariables.clear();
    m_optionalOutputVariables.reserve(m_outNetCDFVariablesOnOff.size());

    for (const auto& pair : m_outNetCDFVariablesOnOff)
    {
      if(pair.second && (m_outNetCDFVariablesIOFunctions.find(pair.first) != m_outNetCDFVariablesIOFunctions.end()))
        m_optionalOutputVariables.push_back(pair.first);
    }

    m_outputNetCDF->sync();

    returnValue = true;
  }
  catch (NcException &e)
  {
    std::string message = std::string(e.what());
    printf("%s\n", e.what());
    errors.push_back(message);
    returnValue = false;
  }


#endif

  return returnValue;
}

bool HTSModel::readInputFileOptionTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  std::string optionsFlag = options[0].toStdString();
  auto it = m_optionsFlags.find(optionsFlag);

  if (it != m_optionsFlags.end())
  {
    int optionsIndex = it->second;

    switch (optionsIndex)
    {
      case 1:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_startDateTime = DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 2:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;

            if (DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_endDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 3:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_outputInterval = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Report interval error";
            return false;
          }
        }
        break;
      case 4:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_maxTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Max time step error";
            return false;
          }
        }
        break;
      case 5:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_minTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Min time step error";
            return false;
          }
        }
        break;
      case 6:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_numInitFixedTimeSteps = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of initial time step error";
            return false;
          }
        }
        break;
      case 7:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useAdaptiveTimeStep = QString::compare(options[1], "No", Qt::CaseInsensitive) && QString::compare(options[1], "False", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use adaptive time step error";
            return false;
          }
        }
        break;
      case 8:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_timeStepRelaxationFactor = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Time step relaxation factor error";
            return false;
          }
        }
        break;
      case 9:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_solverTypeFlags.find(code);

            int heatSolverMode = -1;

            if (it != m_solverTypeFlags.end())
              heatSolverMode = it->second;

            switch (heatSolverMode)
            {
              case 1:
                m_heatSolver->setSolverType(ODESolver::RK4);
                break;
              case 2:
                m_heatSolver->setSolverType(ODESolver::RKQS);
                break;
              case 3:
                m_heatSolver->setSolverType(ODESolver::CVODE_ADAMS);
                break;
              case 4:
                m_heatSolver->setSolverType(ODESolver::CVODE_BDF);
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver type error";
            return false;
          }
        }
        break;
      case 10:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            double abs_tol = options[1].toDouble(&ok);

            if (ok)
              m_heatSolver->setAbsoluteTolerance(abs_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver absolute tolerance error";
            return false;
          }
        }
        break;
      case 11:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            double rel_tol = options[1].toDouble(&ok);

            if (ok)
              m_heatSolver->setRelativeTolerance(rel_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver relative tolerance error";
            return false;
          }
        }
        break;
      case 12:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_waterDensity = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Water density error";
            return false;
          }
        }
        break;
      case 13:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_cp = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Specific heat capacity of water error";
            return false;
          }
        }
        break;
      case 14:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            int numSolutes = options[1].toInt(&parsed);

            if (parsed)
              setNumSolutes(numSolutes);

            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of solutes error";
            return false;
          }
        }
        break;
      case 15:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_verbose = QString::compare(options[1], "No", Qt::CaseInsensitive) && QString::compare(options[1], "False", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Verbose error";
            return false;
          }
        }
        break;
      case 16:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_flushToDiskFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Flush to disk frequency error";
            return false;
          }
        }
        break;
      case 17:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_printFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Print frequency error";
            return false;
          }
        }
        break;
      case 18:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_sedDensity = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Sediment density error";
            return false;
          }
        }
        break;
      case 19:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_sedCp = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Sediment specific heat capacity error";
            return false;
          }
        }
        break;
    }
  }

  return true;
}

bool HTSModel::readInputFileOutputTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  QString optionsFlag = options[0];

  if (options.size() == 2)
  {
    if (!QString::compare(optionsFlag, "csv", Qt::CaseInsensitive))
    {
      m_outputCSVFileInfo = QFileInfo(options[1]);
    }
    else if (!QString::compare(optionsFlag, "netcdf", Qt::CaseInsensitive))
    {
      m_outputNetCDFFileInfo = QFileInfo(options[1]);
    }
  }
  else
  {
    errorMessage = "Output file error";
    return false;
  }

  return true;
}

bool HTSModel::readInputFileSolutesTag(const QString &line, QString &errorMessage)
{
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    bool foundError = false;

    if (m_addedSoluteCount < (int)m_solutes.size())
    {
      m_solutes[m_addedSoluteCount] = columns[0].toStdString();

      std::string solverType = columns[2].toStdString();
      auto it = m_solverTypeFlags.find(solverType);

      if (it != m_solverTypeFlags.end())
      {
        int solverTypeCode = it->second;

        switch (solverTypeCode)
        {
          case 1:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::RK4);
            break;
          case 2:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::RKQS);
            break;
          case 3:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::CVODE_ADAMS);
            break;
          case 4:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::CVODE_BDF);
            break;
          default:
            foundError = true;
            break;
        }

        if (foundError)
        {
          errorMessage = "Solute error";
          return false;
        }

        bool parsed;

        double first_order_k = columns[1].toDouble(&parsed);

        if(parsed)
        {
          m_solute_first_order_k[m_addedSoluteCount] = first_order_k;
        }
        else
        {
          errorMessage = "Invalid solute first order reaction rate";
          return false;
        }

        double abs_tol = columns[3].toDouble(&parsed);

        if (parsed)
        {
          m_soluteSolvers[m_addedSoluteCount]->setAbsoluteTolerance(abs_tol);
        }
        else
        {
          errorMessage = "Solute absolute tolerance error";
          return false;
        }

        double rel_tol = columns[4].toDouble(&parsed);

        if (parsed)
        {
          m_soluteSolvers[m_addedSoluteCount]->setRelativeTolerance(rel_tol);
        }
        else
        {
          errorMessage = "Solute relative tolerance error";
          return false;
        }
      }
      else
      {
        errorMessage = "Solute error";
        return false;
      }

      m_addedSoluteCount++;
    }
    else
    {
      errorMessage = "Solute error count";
      return false;
    }
  }
  else
  {
    errorMessage = "Solute error";
    return false;
  }

  return true;
}

bool HTSModel::readInputFileElementJunctionsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];

    bool workedX;
    bool workedY;
    bool workedZ;

    double x = columns[1].toDouble(&workedX);

    double y = columns[2].toDouble(&workedY);

    double z = columns[3].toDouble(&workedZ);

    if (workedX && workedY && workedZ)
    {
      addElementJunction(id.toStdString(), x, y, z);
    }
    else
    {
      errorMessage = "Junctions error";
      return false;
    }
  }
  else
  {
    errorMessage = "Junctions error";
    return false;
  }

  return true;
}

bool HTSModel::readInputFileElementsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() > 11)
  {
    QString id = columns[0];
    QString fromId = columns[1];
    QString toId = columns[2];

    auto fromIt = m_elementJunctionsById.find(fromId.toStdString());
    auto toIt = m_elementJunctionsById.find(toId.toStdString());

    if (fromIt != m_elementJunctionsById.end() &&
        toIt != m_elementJunctionsById.end())
    {
      ElementJunction *ej1 = m_elementJunctionsById[fromId.toStdString()];
      ElementJunction *ej2 = m_elementJunctionsById[toId.toStdString()];


      bool lengthOk;
      double length = columns[3].toDouble(&lengthOk);

      bool depthOk ;
      double depth = columns[4].toDouble(&depthOk);

      bool widthOk ;
      double width = columns[5].toDouble(&widthOk);

      bool tempOk ;
      double temp = columns[6].toDouble(&tempOk);

      bool advCoeffOk ;
      double advCoeff = columns[7].toDouble(&advCoeffOk);

      bool mcTempOk ;
      double mcTemp = columns[8].toDouble(&mcTempOk);

      bool grndDepthOk;
      double grndDepth = columns[9].toDouble(&grndDepthOk);

      bool grTempOk ;
      double grTemp = columns[10].toDouble(&grTempOk);

      bool thermDiffOk ;
      double thermDiff = columns[11].toDouble(&thermDiffOk);


      if (lengthOk && depthOk && advCoeffOk &&
          tempOk && mcTempOk && grTempOk &&
          thermDiffOk && widthOk && grndDepthOk)
      {
        Element *element = addElement(id.toStdString(), ej1, ej2);
        element->length = length;
        element->depth = depth;
        element->mainChannelAdvectionCoeff = advCoeff;
        element->width = width;
        element->temperature.value = temp;
        element->groundConductionDepth = grndDepth;
        element->mainChannelTemperature = mcTemp;
        element->groundTemperature = grTemp;
        element->sedThermalDiffCoefficient = thermDiff;

        if(m_solutes.size() && columns.size() > 11)
        {
          int soluteIndex = 0;

          for (int i = 12; i < columns.size(); i = i+3)
          {
            bool soluteOk ;
            double solute = columns[i].toDouble(&soluteOk);

            bool mainChannelSoluteOk;
            double mainChannelSolute = columns[i+1].toDouble(&mainChannelSoluteOk);

            bool soluteDiffCoeffOk;
            double soluteDiffCoeff = columns[i+2].toDouble(&soluteDiffCoeffOk);


            if (soluteOk && mainChannelSoluteOk &&  soluteDiffCoeffOk && soluteIndex < (int)m_solutes.size())
            {
              element->soluteConcs[soluteIndex].value = solute;
              element->mainChannelSoluteConcs[soluteIndex] = mainChannelSolute;
              element->sedSoluteDiffCoefficients[soluteIndex] = soluteDiffCoeff;
            }
            else
            {
              errorMessage = "Wrong initial solute parameters";
              return false;
            }

            soluteIndex++;
          }
        }
      }
      else
      {
        errorMessage = "";
        return false;
      }
    }
    else
    {
      errorMessage = "Wrong upstream or downstream junction";
      return false;
    }
  }

  return true;
}

bool HTSModel::readInputFileSourcesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 7)
  {
    QString idFrom = columns[0];
    QString idTo = columns[2];

    bool okStart;
    double startFactor = columns[1].toDouble(&okStart);

    bool okEnd;
    double endFactor = columns[3].toDouble(&okEnd);

    auto itFrom = m_elementsById.find(idFrom.toStdString());
    auto itTo = m_elementsById.find(idTo.toStdString());

    if (okStart && okEnd && itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *elementFrom = itFrom->second;
      Element *elementTo = itTo->second;

      QString variable = columns[4].trimmed();
      QString type = columns[5];
      SourceBC::VariableType variableType;
      int soluteIndex = -1;

      if (!QString::compare(variable, "HEAT", Qt::CaseInsensitive))
      {
        variableType = SourceBC::HeatSource;
        soluteIndex = 0;
      }
      else if (!QString::compare(variable, "FlOW", Qt::CaseInsensitive))
      {
        variableType = SourceBC::FlowSource;
        soluteIndex = 0;
      }
      else
      {
        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            soluteIndex = i;
            variableType = SourceBC::SoluteSource;
            break;
          }
        }
      }

      if (soluteIndex > -1)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[6].toDouble(&valueOk);

          if (valueOk)
          {
            SourceBC *sourceBC = new SourceBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);
            sourceBC->setSoluteIndex(soluteIndex);
            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);

            m_timeSeries[ts->id().toStdString()] = ts;

            sourceBC->setTimeSeries(ts);

            m_boundaryConditions.push_back(sourceBC);
          }
          else
          {
            errorMessage = "Source is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
        {

          std::string tsId = columns[6].toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            SourceBC *sourceBC = new SourceBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);
            sourceBC->setSoluteIndex(soluteIndex);
            sourceBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(sourceBC);
          }
          else
          {
            errorMessage = "Source is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Source is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Source is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool HTSModel::readInputFileHydraulicsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString variableType = columns[2];
      QString valueType = columns[3];
      QString varValue = columns[4].trimmed();

      auto it = m_hydraulicVariableFlags.find(variableType.toStdString());

      if (it != m_hydraulicVariableFlags.end())
      {
        int variableIndex = it->second;

        if (!QString::compare(valueType, "VALUE", Qt::CaseInsensitive))
        {

          bool valueOk;
          double value =  varValue.toDouble(&valueOk);

          if (valueOk)
          {
            HydraulicsBC *hydraulicsBC = new HydraulicsBC(fromElement, toElement, variableIndex, this);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);

            m_timeSeries[ts->id().toStdString()] = ts;

            hydraulicsBC->setTimeSeries(ts);

            m_boundaryConditions.push_back(hydraulicsBC);
          }
          else
          {
            errorMessage = "hydraulics value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = varValue.toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            HydraulicsBC *hydraulicsBC = new HydraulicsBC(fromElement, toElement, variableIndex, this);
            hydraulicsBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(hydraulicsBC);
          }
          else
          {
            errorMessage = "Specified hydraulics timeseries does not exist";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Variable specified for uniform hydraulic is incorrect";
        return false;
      }
    }
    else
    {
      errorMessage = "Uniform hydraulic condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool HTSModel::readInputFileRadiativeFluxesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString type = columns[2];
      QString varValue = columns[3].trimmed();

      if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
      {

        bool valueOk;
        double value =  varValue.toDouble(&valueOk);

        if (valueOk)
        {
          RadiativeFluxBC *radiationFluxBC = new RadiativeFluxBC(fromElement, toElement, this);

          QUuid uid = QUuid::createUuid();
          QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
          ts->addRow(m_startDateTime, value);
          ts->addRow(m_endDateTime, value);
          m_timeSeries[ts->id().toStdString()] = ts;

          radiationFluxBC->setTimeSeries(ts);
          m_boundaryConditions.push_back(radiationFluxBC);
        }
        else
        {
          errorMessage = "Radiation BC value is invalid";
          return false;
        }
      }
      else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
      {
        std::string tsId = varValue.toStdString();
        auto tsIt = m_timeSeries.find(tsId);

        if(tsIt != m_timeSeries.end())
        {
          RadiativeFluxBC *radiationFluxBC = new RadiativeFluxBC(fromElement, toElement, this);
          radiationFluxBC->setTimeSeries(tsIt->second);
          m_boundaryConditions.push_back(radiationFluxBC);
        }
        else
        {
          errorMessage = "Specified timeseries does not exist";
          return false;
        }
      }
    }
    else
    {
      errorMessage = "Boundary condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool HTSModel::readInputFileElementBCTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 6)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      bool sourceOk = false;
      ElementBC::BCSource sourceIndex;

      QString source = columns[3];
      if(!source.compare("MC", Qt::CaseInsensitive))
      {
        sourceIndex = ElementBC::CHANNEL;
        sourceOk = true;
      }
      else if(!source.compare("GR", Qt::CaseInsensitive))
      {
        sourceIndex = ElementBC::GROUND;
        sourceOk = true;
      }

      bool variableOk = false;
      int variableIndex;
      QString variableType = columns[2];

      if(!QString::compare(variableType,"Temperature", Qt::CaseInsensitive))
      {
        variableIndex = -1;
        variableOk = true;
      }
      else
      {
        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variableType.toStdString()))
          {
            variableIndex = i;
            variableOk = true;
          }
        }
      }

      QString valueType = columns[4];
      QString varValue = columns[5].trimmed();

      if(variableOk && sourceOk)
      {
        if (!QString::compare(valueType, "VALUE", Qt::CaseInsensitive))
        {

          bool valueOk;
          double value =  varValue.toDouble(&valueOk);

          if (valueOk)
          {
            ElementBC *elementBC = new ElementBC(fromElement, toElement, sourceIndex,
                                                 variableIndex, this);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);
            m_timeSeries[ts->id().toStdString()] = ts;

            elementBC->setTimeSeries(ts);
            m_boundaryConditions.push_back(elementBC);
          }
          else
          {
            errorMessage = "Uniform bc value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = varValue.toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if(tsIt != m_timeSeries.end())
          {
            ElementBC *elementBC = new ElementBC(fromElement, toElement, sourceIndex,
                                                 variableIndex, this);
            elementBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(elementBC);
          }
          else
          {
            errorMessage = "Specified uniform meteorology filepath does not exist";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Variable specified for uniform meteorology is incorrect";
        return false;
      }
    }
    else
    {
      errorMessage = "Uniform meteorology boundary condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool HTSModel::readInputFileTimeSeriesTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2)
  {
    QFileInfo fileInfo(options[1].trimmed());

    if (fileInfo.isRelative())
      fileInfo = relativePathToAbsolute(fileInfo);

    if(QFile::exists(fileInfo.absoluteFilePath()))
    {
      QSharedPointer<TimeSeries> timeSeries(TimeSeries::createTimeSeries(options[0], fileInfo, this));

      if(!timeSeries.isNull())
      {
        m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
      }
      else
      {
        errorMessage = "Timeseries specified is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified filepath does not exist";
      return false;
    }
  }
  else
  {
    errorMessage = "TimeSeries must have two columns";
    return false;
  }

  return true;
}

bool HTSModel::readOutputVariableOnOff(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2)
  {
    m_outNetCDFVariablesOnOff[options[0].toStdString()] = !QString::compare(options[1].toLower(),"no") ||
        !QString::compare(options[1].toLower(),"false") ? false : true;
  }
  else {
    errorMessage = "Output variable must have two columns";
    return false;
  }

  return true;
}

void HTSModel::writeOutput()
{
  m_currentflushToDiskCount++;

  if (m_currentflushToDiskCount >= m_flushToDiskFrequency)
  {
    m_flushToDisk = true;
    m_currentflushToDiskCount = 0;
  }
  else
  {
    m_flushToDisk = false;
  }

  writeCSVOutput();
  writeNetCDFOutput();
}

void HTSModel::writeCSVOutput()
{
  if(m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    for (size_t i = 0; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->index
                        << ", " << element->x << ", " << element->y << ", " << element->z
                        << ", " << element->depth
                        << ", " << element->width
                        << ", " << element->xSectionArea
                        << ", " << element->temperature.value
                        << ", " << element->totalExternalHeatFluxesBalance
                        << ", " << element->totalRadiationFluxesHeatBalance
                        << ", " << element->totalHeatBalance;

      for (size_t j = 0; j < m_solutes.size(); j++)
      {
        m_outputCSVStream << ", " << element->soluteConcs[j].value
                          << ", " << element->totalExternalSoluteFluxesMassBalance[j]
                             << ", " << element->totalSoluteMassBalance[j];
      }

      m_outputCSVStream << endl;
    }

    if (m_flushToDisk)
    {
      m_outputCSVStream.flush();
    }
  }
}

void HTSModel::writeNetCDFOutput()
{
#ifdef USE_NETCDF

  if(m_outputNetCDF)
  {
    size_t currentTime = m_outNetCDFVariables["time"].getDim(0).getSize();

    //Set current dateTime
    m_outNetCDFVariables["time"].putVar(std::vector<size_t>({currentTime}), m_currentDateTime);

    int nVars = static_cast<int>(m_optionalOutputVariables.size());

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nVars; i++)
    {
      std::string varName = m_optionalOutputVariables[static_cast<size_t>(i)];
      (m_outNetCDFVariablesIOFunctions[varName])(currentTime, m_outNetCDFVariables[varName], m_elements);
    }


    if(m_outNetCDFVariablesOnOff["total_heat_balance"])
      m_outNetCDFVariables["total_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalHeatBalance);

    if(m_outNetCDFVariablesOnOff["total_radiation_flux_heat_balance"])
      m_outNetCDFVariables["total_radiation_flux_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalRadiationHeatBalance);

    if(m_outNetCDFVariablesOnOff["total_external_heat_flux_balance"])
      m_outNetCDFVariables["total_external_heat_flux_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalExternalHeatFluxBalance);


    if(numSolutes())
    {

      if(m_outNetCDFVariablesOnOff["total_solute_mass_balance"])
        m_outNetCDFVariables["total_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalSoluteMassBalance.data());

      if(m_outNetCDFVariablesOnOff["total_external_solute_flux_mass_balance"])
        m_outNetCDFVariables["total_external_solute_flux_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalExternalSoluteFluxMassBalance.data());

    }

    if (m_flushToDisk)
    {
      m_outputNetCDF->sync();
    }
  }

#endif
}

void HTSModel::closeOutputFiles()
{
  closeCSVOutputFile();
  closeOutputNetCDFFile();
}

void HTSModel::closeCSVOutputFile()
{
  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }
}

void HTSModel::closeOutputNetCDFFile()
{
#ifdef USE_NETCDF

  if (m_outputNetCDF && !m_outputNetCDF->ncFile()->isNull())
  {
    m_outputNetCDF->sync();
    delete m_outputNetCDF;
    m_outputNetCDF = nullptr;
  }

#endif

}

QFileInfo HTSModel::relativePathToAbsolute(const QFileInfo &fileInfo)
{
  if (fileInfo.isRelative())
  {
    if (!m_inputFile.filePath().isEmpty() &&
        !m_inputFile.filePath().isNull() &&
        QFile::exists(m_inputFile.absoluteFilePath()))
    {
      QFileInfo absoluteFilePath = m_inputFile.absoluteDir().absoluteFilePath(fileInfo.filePath());

      if (absoluteFilePath.absoluteDir().exists())
      {
        return absoluteFilePath;
      }
    }
  }

  return fileInfo;
}

const unordered_map<string, int> HTSModel::m_inputFileFlags({
                                                              {"[OPTIONS]", 1},
                                                              {"[OUTPUTS]", 2},
                                                              {"[SOLUTES]", 3},
                                                              {"[ELEMENTJUNCTIONS]", 4},
                                                              {"[ELEMENTS]", 5},
                                                              {"[SOURCES]", 6},
                                                              {"[HYDRAULICS]", 7},
                                                              {"[RADIATIVE_FLUXES]", 8},
                                                              {"[BOUNDARY_CONDITIONS]", 9},
                                                              {"[TIMESERIES]", 10},
                                                              {"[OUTPUTVARIABLES]", 11},
                                                            });

const unordered_map<string, int> HTSModel::m_optionsFlags({
                                                            {"START_DATETIME", 1},
                                                            {"END_DATETIME", 2},
                                                            {"REPORT_INTERVAL", 3},
                                                            {"MAX_TIME_STEP", 4},
                                                            {"MIN_TIME_STEP", 5},
                                                            {"NUM_INITIAL_FIXED_STEPS", 6},
                                                            {"USE_ADAPTIVE_TIME_STEP", 7},
                                                            {"TIME_STEP_RELAXATION_FACTOR", 8},
                                                            {"TEMP_SOLVER", 9},
                                                            {"TEMP_SOLVER_ABS_TOL", 10},
                                                            {"TEMP_SOLVER_REL_TOL", 11},
                                                            {"WATER_DENSITY", 12},
                                                            {"WATER_SPECIFIC_HEAT_CAPACITY", 13},
                                                            {"NUM_SOLUTES", 14},
                                                            {"VERBOSE", 15},
                                                            {"FLUSH_TO_DISK_FREQ", 16},
                                                            {"PRINT_FREQ", 17},
                                                            {"SEDIMENT_DENSITY", 18},
                                                            {"SEDIMENT_SPECIFIC_HEAT_CAPACITY", 19},
                                                          });

const unordered_map<string, int> HTSModel::m_solverTypeFlags({{"RK4", 1},
                                                              {"RKQS", 2},
                                                              {"ADAMS", 3},
                                                              {"BDF", 4}});

const unordered_map<string, int> HTSModel::m_hydraulicVariableFlags({{"DEPTH", 1},
                                                                     {"WIDTH", 2},
                                                                     {"MC_ADVECTION_COEFFICIENT", 3},
                                                                     {"GR_DEPTH",4}});



const QRegExp HTSModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");
