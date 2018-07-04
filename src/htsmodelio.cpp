#include "htsmodel.h"
#include "element.h"
#include "elementjunction.h"
#include "pointsrctimeseriesbc.h"
#include "nonpointsrctimeseriesbc.h"
#include "hydraulicstimeseriesbc.h"
#include "radiativefluxtimeseriesbc.h"
#include "boundarycondition.h"
#include "temporal/timedata.h"
#include "threadsafencfile.h"
#include "timeseries.h"

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
                readSuccess = readInputFilePointSourcesTag(line, error);
                break;
              case 7:
                readSuccess = readInputFileNonPointSourcesTag(line, error);
                break;
              case 8:
                readSuccess = readInputFileUniformHydraulicsTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileNonUniformHydraulicsTag(line, error);
                break;
              case 10:
                readSuccess = readInputFileUniformRadiativeFluxesTag(line, error);
                break;
              case 11:
                readSuccess = readInputFileNonUniformRadiativeFluxesTag(line, error);
                break;
              case 12:
                readSuccess = readInputFileUniformBoundaryConditionTag(line, error);
                break;
              case 13:
                readSuccess = readInputFileNonUniformBoundaryConditionTag(line, error);
                break;
            }
          }

          if (!readSuccess)
          {
            errors.push_back("Line " + std::to_string(lineCount) + " : " + error.toStdString());
            file.close();
            return false;
            break;
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
    }

    m_outputCSVStream.flush();

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

#ifdef USE_OPENMP
#pragma omp critical
#endif
  {

    closeOutputNetCDFFile();

    try
    {

      m_outputNetCDF = new ThreadSafeNcFile(m_outputNetCDFFileInfo.absoluteFilePath().toStdString(), NcFile::replace);

      //time variable
      NcDim timeDim =  m_outputNetCDF->ncFile()->addDim("time");
      NcVar timeVar =  m_outputNetCDF->ncFile()->addVar("time", NcType::nc_DOUBLE, timeDim);
      timeVar.putAtt("time:long_name", "time");
      timeVar.putAtt("time:units", "days since 1858-11-17 0:0:0");
      timeVar.putAtt("time:calendar", "modified_julian");

      //Add Solutes
      NcDim solutesDim =  m_outputNetCDF->ncFile()->addDim("solutes", m_solutes.size());
      NcVar solutes =  m_outputNetCDF->ncFile()->addVar("solute_names", NcType::nc_STRING, solutesDim);
      solutes.putAtt("solute_names::long_name", "solute names");

      if (m_solutes.size())
      {
        char **soluteNames = new char *[m_solutes.size()];

        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          string soluteName = m_solutes[i];
          soluteNames[i] = new char[soluteName.size() + 1];
          std::strcpy(soluteNames[i], soluteName.c_str());
        }

        solutes.putVar(soluteNames);

        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          delete[] soluteNames[i];
        }

        delete[] soluteNames;
      }

      //Add element junctions
      NcDim junctionDim =  m_outputNetCDF->ncFile()->addDim("element_junctions", m_elementJunctions.size());

      NcVar junctionIdentifiers =  m_outputNetCDF->ncFile()->addVar("element_junction_id", NcType::nc_STRING, junctionDim);
      junctionIdentifiers.putAtt("element_junction_id::long_name", "element junction identifier");

      NcVar junctionX =  m_outputNetCDF->ncFile()->addVar("x", NcType::nc_DOUBLE, junctionDim);
      junctionX.putAtt("x:long_name", "junction x-coordinate");
      junctionX.putAtt("x:units", "m");

      NcVar junctionY =  m_outputNetCDF->ncFile()->addVar("y", NcType::nc_DOUBLE, junctionDim);
      junctionY.putAtt("y:long_name", "junction y-coordinate");
      junctionY.putAtt("y:units", "m");

      NcVar junctionZ =  m_outputNetCDF->ncFile()->addVar("z", NcType::nc_DOUBLE, junctionDim);
      junctionZ.putAtt("z:long_name", "junction z-coordinate");
      junctionZ.putAtt("z:units", "m");

      double *vertx = new double[m_elementJunctions.size()];
      double *verty = new double[m_elementJunctions.size()];
      double *vertz = new double[m_elementJunctions.size()];
      char **junctionIds = new char *[m_elementJunctions.size()];

      //write other relevant junction attributes here.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (size_t i = 0; i < m_elementJunctions.size(); i++)
      {
        ElementJunction *junction = m_elementJunctions[i];

        junctionIds[i] = new char[junction->id.size() + 1];
        std::strcpy(junctionIds[i], junction->id.c_str());

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
      NcDim elementsDim =  m_outputNetCDF->ncFile()->addDim("elements", m_elements.size());

      NcVar elementIdentifiers =  m_outputNetCDF->ncFile()->addVar("element_id", NcType::nc_STRING, elementsDim);
      elementIdentifiers.putAtt("element_id::long_name", "element identifier");

      NcVar elementFromJunction =  m_outputNetCDF->ncFile()->addVar("from_junction", NcType::nc_INT64, elementsDim);
      elementFromJunction.putAtt("from_junction:long_name", "upstream junction");

      NcVar elementToJunction =  m_outputNetCDF->ncFile()->addVar("to_junction", NcType::nc_INT64, elementsDim);
      elementToJunction.putAtt("to_junction:long_name", "downstream junction");

      int *fromJunctions = new int[m_elements.size()];
      int *toJunctions = new int[m_elements.size()];
      char **elementIds = new char *[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (size_t i = 0; i < m_elements.size(); i++)
      {
        Element *element = m_elements[i];

        elementIds[i] = new char[element->id.size() + 1];
        std::strcpy(elementIds[i], element->id.c_str());

        fromJunctions[i] = element->upstreamJunction->index;
        toJunctions[i] = element->downstreamJunction->index;
      }

      elementIdentifiers.putVar(elementIds);
      elementFromJunction.putVar(fromJunctions);
      elementToJunction.putVar(toJunctions);

      delete[] fromJunctions;
      delete[] toJunctions;

      for (size_t i = 0; i < m_elements.size(); i++)
      {
        delete[] elementIds[i];
      }

      delete[] elementIds;


      NcVar depthVar =  m_outputNetCDF->ncFile()->addVar("depth", "double",
                                                         std::vector<std::string>({"time", "elements"}));
      depthVar.putAtt("depth:long_name", "channel flow depth");
      depthVar.putAtt("depth:units", "m");

      NcVar widthVar =  m_outputNetCDF->ncFile()->addVar("width", "double",
                                                         std::vector<std::string>({"time", "elements"}));
      widthVar.putAtt("width:long_name", "channel flow top width");
      widthVar.putAtt("width:units", "m");

      NcVar xsectAreaVar =  m_outputNetCDF->ncFile()->addVar("xsection_area", "double",
                                                             std::vector<std::string>({"time", "elements"}));
      xsectAreaVar.putAtt("xsection_area:long_name", "flow cross-section area");
      xsectAreaVar.putAtt("xsection_area:units", "m^2");

      NcVar temperatureVar =  m_outputNetCDF->ncFile()->addVar("temperature", "double",
                                                               std::vector<std::string>({"time", "elements"}));
      temperatureVar.putAtt("temperature:long_name", "temperature");
      temperatureVar.putAtt("temperature:units", "°C");


      NcVar channelConductionFluxVar =  m_outputNetCDF->ncFile()->addVar("channel_conduction_heat_flux", "double",
                                                                         std::vector<std::string>({"time", "elements"}));
      channelConductionFluxVar.putAtt("channel_conduction_heat_flux:long_name", "channel conduction heat flux");
      channelConductionFluxVar.putAtt("channel_conduction_heat_flux:units", "W/m2");


      NcVar groundConductionFluxVar =  m_outputNetCDF->ncFile()->addVar("ground_conduction_heat_flux", "double",
                                                                        std::vector<std::string>({"time", "elements"}));
      groundConductionFluxVar.putAtt("ground_conduction_heat_flux:long_name", "ground conduction heat flux");
      groundConductionFluxVar.putAtt("ground_conduction_heat_flux:units", "W/m2");


      NcVar channelAdvectionHeatFluxVar =  m_outputNetCDF->ncFile()->addVar("channel_advection_heat_flux", "double",
                                                                            std::vector<std::string>({"time", "elements"}));
      channelAdvectionHeatFluxVar.putAtt("channel_advection_heat_flux:long_name", "channel advection heat flux");
      channelAdvectionHeatFluxVar.putAtt("channel_advection_heat_flux:units", "W/m2");


      NcVar totalHeatBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_heat_balance", "double",
                                                                    std::vector<std::string>({"time"}));
      totalHeatBalanceVar.putAtt("total_heat_balance:long_name", "total heat balance");
      totalHeatBalanceVar.putAtt("total_heat_balance:units", "KJ");


      NcVar totalRadiationFluxHeatBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_radiation_flux_heat_balance", "double",
                                                                                 std::vector<std::string>({"time"}));
      totalRadiationFluxHeatBalanceVar.putAtt("total_radiation_flux_heat_balance:long_name", "total radiation flux heat balance");
      totalRadiationFluxHeatBalanceVar.putAtt("total_radiation_flux_heat_balance:units", "KJ");

      NcVar totalExternalHeatFluxBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_external_heat_flux_balance", "double",
                                                                                std::vector<std::string>({"time"}));
      totalExternalHeatFluxBalanceVar.putAtt("total_external_heat_flux_balance:long_name", "total external heat flux balance");
      totalExternalHeatFluxBalanceVar.putAtt("total_external_heat_flux_balance:units", "KJ");

      NcVar totalElementHeatBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_element_heat_balance", "double",
                                                                           std::vector<std::string>({"time", "elements"}));
      totalElementHeatBalanceVar.putAtt("total_element_heat_balance:long_name", "total heat balance");
      totalElementHeatBalanceVar.putAtt("total_element_heat_balance:units", "KJ");


      NcVar totalElementRadiationFluxHeatBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_element_radiation_flux_heat_balance", "double",
                                                                                        std::vector<std::string>({"time", "elements"}));
      totalElementRadiationFluxHeatBalanceVar.putAtt("total_element_radiation_flux_heat_balance:long_name", "total element radiation flux heat balance");
      totalElementRadiationFluxHeatBalanceVar.putAtt("total_element_radiation_flux_heat_balance:units", "KJ");

      NcVar totalElementExternalHeatFluxBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_element_external_heat_flux_balance", "double",
                                                                                       std::vector<std::string>({"time", "elements"}));
      totalElementExternalHeatFluxBalanceVar.putAtt("total_element_external_heat_flux_balance:long_name", "total element external heat flux balance");
      totalElementExternalHeatFluxBalanceVar.putAtt("total_element_external_heat_flux_balance:units", "KJ");

      NcVar solutesVar =  m_outputNetCDF->ncFile()->addVar("solute_concentration", "double",
                                                           std::vector<std::string>({"time", "solutes", "elements"}));
      solutesVar.putAtt("solute_concentration:long_name", "solute concentrations");
      solutesVar.putAtt("solute_concentration:units", "kg/m^3");

      NcVar totalSoluteMassBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_solute_mass_balance", "double",
                                                                          std::vector<std::string>({"time", "solutes"}));
      totalSoluteMassBalanceVar.putAtt("total_solute_mass_balance:long_name", "total solute mass balance");
      totalSoluteMassBalanceVar.putAtt("total_solute_mass_balance:units", "kg");

      NcVar totalExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_external_solute_flux_mass_balance", "double",
                                                                                      std::vector<std::string>({"time", "solutes"}));
      totalExternalSoluteFluxMassBalanceVar.putAtt("total_external_solute_flux_mass_balance:long_name", "total external solute flux mass balance");
      totalExternalSoluteFluxMassBalanceVar.putAtt("total_external_solute_flux_mass_balance:units", "kg");

      NcVar totalElementSoluteMassBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_element_solute_mass_balance", "double",
                                                                                 std::vector<std::string>({"time", "solutes", "elements"}));
      totalElementSoluteMassBalanceVar.putAtt("total_element_solute_mass_balance:long_name", "total element solute mass balance");
      totalElementSoluteMassBalanceVar.putAtt("total_element_solute_mass_balance:units", "kg");

      NcVar totalElementExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->ncFile()->addVar("total_element_external_solute_flux_mass_balance", "double",
                                                                                             std::vector<std::string>({"time", "solutes"}));
      totalElementExternalSoluteFluxMassBalanceVar.putAtt("total_element_external_solute_flux_mass_balance:long_name", "total external solute flux mass balance");
      totalElementExternalSoluteFluxMassBalanceVar.putAtt("total_element_external_solute_flux_mass_balance:units", "kg");

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
            m_useAdaptiveTimeStep = QString::compare(options[1], "No", Qt::CaseInsensitive);
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
            m_verbose = QString::compare(options[1], "No", Qt::CaseInsensitive);
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

  if (columns.size() == 4)
  {
    bool foundError = false;

    if (m_addedSoluteCount < (int)m_solutes.size())
    {
      m_solutes[m_addedSoluteCount] = columns[0].toStdString();

      std::string solverType = columns[1].toStdString();
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
        double abs_tol = columns[2].toDouble(&parsed);

        if (parsed)
        {
          m_soluteSolvers[m_addedSoluteCount]->setAbsoluteTolerance(abs_tol);
        }
        else
        {
          errorMessage = "Solute absolute tolerance error";
          return false;
        }

        double rel_tol = columns[3].toDouble(&parsed);

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

        if (columns.size() > 12)
        {
          for (int i = 13; i < columns.size(); i = i+4)
          {
            bool soluteOk ;
            double solute = columns[i].toDouble(&soluteOk);

            bool mainChannelSoluteOk;
            double mainChannelSolute = columns[i+1].toDouble(&mainChannelSoluteOk);

            bool groundSoluteOk;
            double groundSolute = columns[i+2].toDouble(&groundSoluteOk);

            bool soluteDiffCoeffOk;
            double soluteDiffCoeff = columns[i+3].toDouble(&soluteDiffCoeffOk);


            if (soluteOk && mainChannelSoluteOk &&
                groundSoluteOk &&  soluteDiffCoeffOk && i - 13 < (int)m_solutes.size())
            {
              element->soluteConcs[i - 13].value = solute;
              element->groundSoluteConcs[i-13] = groundSolute;
              element->mainChannelSoluteConcs[i-13] = mainChannelSolute;
              element->sedSoluteDiffCoefficients[i-13] = soluteDiffCoeff;
            }
            else
            {
              errorMessage = "Wrong initial solute parameters";
              return false;
            }
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

bool HTSModel::readInputFilePointSourcesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];
    auto it = m_elementsById.find(id.toStdString());

    if (it != m_elementsById.end())
    {
      Element *element = it->second;
      QString variable = columns[1].trimmed();
      QString type = columns[2];

      int variableIndex = -3;

      if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
      {
        variableIndex = -1;
      }
      else
      {
        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            variableIndex = i;
            break;
          }
        }
      }

      if (variableIndex > -2)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[3].toDouble(&valueOk);

          if (valueOk)
          {
            PointSrcTimeSeriesBC *pointSrcTSBC = new PointSrcTimeSeriesBC(element, variableIndex, this);
            pointSrcTSBC->addValue(m_startDateTime, value);
            pointSrcTSBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(pointSrcTSBC);
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
        {

          QString filePath = columns[3];

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                PointSrcTimeSeriesBC *pointSrcTSBC = new PointSrcTimeSeriesBC(element, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  pointSrcTSBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(pointSrcTSBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Point source is invalid";
                return false;
              }
            }
            else
            {
              errorMessage = "Point source is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Point source is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Point source is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool HTSModel::readInputFileNonPointSourcesTag(const QString &line, QString &errorMessage)
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
      int variableIndex = -2;

      if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
      {
        variableIndex = -1;
      }
      else
      {
        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            variableIndex = i;
            break;
          }
        }
      }

      if (variableIndex > -2)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[6].toDouble(&valueOk);

          if (valueOk)
          {
            NonPointSrcTimeSeriesBC *nonPointSrcTSBC = new NonPointSrcTimeSeriesBC(elementFrom, startFactor, elementTo, endFactor, variableIndex, this);
            nonPointSrcTSBC->addValue(m_startDateTime, value);
            nonPointSrcTSBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(nonPointSrcTSBC);
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
        {

          QString filePath = columns[6];

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                NonPointSrcTimeSeriesBC *nonPointSrcTSBC = new NonPointSrcTimeSeriesBC(elementFrom, startFactor, elementTo, endFactor, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  nonPointSrcTSBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(nonPointSrcTSBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Point source is invalid";
                return false;
              }
            }
            else
            {
              errorMessage = "Point source is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Point source is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Point source is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool HTSModel::readInputFileUniformHydraulicsTag(const QString &line, QString &errorMessage)
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
            UniformHydraulicsTimeSeriesBC *uniformHydraulicsBC = new UniformHydraulicsTimeSeriesBC(fromElement, toElement,
                                                                                                   variableIndex, this);
            uniformHydraulicsBC->addValue(m_startDateTime, value);
            uniformHydraulicsBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(uniformHydraulicsBC);
          }
          else
          {
            errorMessage = "Uniform hydraulics value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "FILE", Qt::CaseInsensitive))
        {
          QString filePath = varValue;

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                UniformHydraulicsTimeSeriesBC *uniformHydraulicsBC = new UniformHydraulicsTimeSeriesBC(fromElement, toElement,
                                                                                                       variableIndex, this);
                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  uniformHydraulicsBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(uniformHydraulicsBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified uniform hydraulics filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified uniform hydraulics filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified uniform hydraulics filepath does not exist";
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

bool HTSModel::readInputFileNonUniformHydraulicsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 2)
  {
    auto it = m_hydraulicVariableFlags.find(columns[0].toStdString());

    if (it != m_hydraulicVariableFlags.end())
    {
      int variableIndex = it->second;

      QString filePath = columns[1];

      if (!filePath.isEmpty() && !filePath.isNull())
      {
        QFileInfo fileInfo(filePath);

        if (fileInfo.isRelative())
          fileInfo = relativePathToAbsolute(fileInfo);

        if (QFile::exists(fileInfo.absoluteFilePath()))
        {
          std::map<double, std::vector<double>> timeSeries;
          std::vector<std::string> headers;

          if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
          {
            for (size_t i = 0; i < headers.size(); i++)
            {
              auto eit = m_elementsById.find(headers[i]);

              if (eit != m_elementsById.end())
              {
                Element *element = eit->second;
                HydraulicsTimeSeriesBC *hydraulicsTimeSeries = new HydraulicsTimeSeriesBC(element, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[i];
                  hydraulicsTimeSeries->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(hydraulicsTimeSeries);
              }
              else
              {
                errorMessage = "Specified time varying hydraulic file is invalid";
                return false;
              }
            }
          }
          else
          {
            errorMessage = "Specified time varying hydraulic file is invalid";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified time varying hydraulic file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying hydraulic file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying hydraulic file is invalid";
      return false;
    }
  }
  else
  {
    errorMessage = "Specified time varying hydraulic file is invalid";
    return false;
  }

  return true;
}

bool HTSModel::readInputFileUniformRadiativeFluxesTag(const QString &line, QString &errorMessage)
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
          UniformRadiativeFluxTimeSeriesBC *uniformRadiationFluxBC = new UniformRadiativeFluxTimeSeriesBC(fromElement, toElement, this);
          uniformRadiationFluxBC->addValue(m_startDateTime, value);
          uniformRadiationFluxBC->addValue(m_endDateTime, value);
          m_boundaryConditions.push_back(uniformRadiationFluxBC);
        }
        else
        {
          errorMessage = "Radiation BC value is invalid";
          return false;
        }
      }
      else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
      {
        QString filePath = varValue;

        if (!filePath.isEmpty() && !filePath.isNull())
        {
          QFileInfo fileInfo(filePath);

          if (fileInfo.isRelative())
            fileInfo = relativePathToAbsolute(fileInfo);

          if (QFile::exists(fileInfo.absoluteFilePath()))
          {
            std::map<double, std::vector<double>> timeSeries;
            std::vector<std::string> headers;

            if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
            {
              UniformRadiativeFluxTimeSeriesBC *uniformRadiationFluxBC = new UniformRadiativeFluxTimeSeriesBC(fromElement, toElement, this);

              for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
              {
                double dateTime = it->first;
                double value = it->second[0];
                uniformRadiationFluxBC->addValue(dateTime, value);
              }

              m_boundaryConditions.push_back(uniformRadiationFluxBC);

              timeSeries.clear();
              headers.clear();
            }
            else
            {
              errorMessage = "Specified BC filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified BC filepath does not exist";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified BC filepath does not exist";
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

bool HTSModel::readInputFileNonUniformRadiativeFluxesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 1)
  {
    QString filePath = columns[0];

    if (!filePath.isEmpty() && !filePath.isNull())
    {
      QFileInfo fileInfo(filePath);

      if (fileInfo.isRelative())
        fileInfo = relativePathToAbsolute(fileInfo);

      if (QFile::exists(fileInfo.absoluteFilePath()))
      {
        std::map<double, std::vector<double>> timeSeries;
        std::vector<std::string> headers;

        if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
        {
          for (size_t i = 0; i < headers.size(); i++)
          {
            auto eit = m_elementsById.find(headers[i]);

            if (eit != m_elementsById.end())
            {
              Element *element = eit->second;
              RadiativeFluxTimeSeriesBC *radiativeFluxTimeSeries = new RadiativeFluxTimeSeriesBC(element, this);

              for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
              {
                double dateTime = it->first;
                double value = it->second[i];
                radiativeFluxTimeSeries->addValue(dateTime, value);
              }

              m_boundaryConditions.push_back(radiativeFluxTimeSeries);
            }
            else
            {
              errorMessage = "Specified time varying radiative flux file is invalid";
              return false;
            }
          }
        }
        else
        {
          errorMessage = "Specified time varying radiative flux file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying radiative flux file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying radiative flux file is invalid";
      return false;
    }

  }
  else
  {
    errorMessage = "Specified time varying radiative flux file is invalid";
    return false;
  }

  return true;
}

bool HTSModel::readInputFileUniformBoundaryConditionTag(const QString &line, QString &errorMessage)
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
      int sourceIndex = -1;

      QString source = columns[3];
      if(!source.compare("MC", Qt::CaseInsensitive))
      {
        sourceIndex = 1;
        sourceOk = true;
      }
      else if(!source.compare("GR", Qt::CaseInsensitive))
      {
        sourceIndex = 2;
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
            UniformBoundaryCondition *uniformBC = new UniformBoundaryCondition(fromElement, toElement, sourceIndex,
                                                                               variableIndex, this);
            uniformBC->addValue(m_startDateTime, value);
            uniformBC->addValue(m_endDateTime, value);

            m_boundaryConditions.push_back(uniformBC);
          }
          else
          {
            errorMessage = "Uniform bc value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "FILE", Qt::CaseInsensitive))
        {
          QString filePath = varValue;

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                UniformBoundaryCondition *uniformBC = new UniformBoundaryCondition(fromElement, toElement, sourceIndex,
                                                                                   variableIndex, this);
                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  uniformBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(uniformBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified uniform meteorology filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified uniform meteorology filepath does not exist";
              return false;
            }
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

bool HTSModel::readInputFileNonUniformBoundaryConditionTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 3)
  {

    bool sourceOk = false;
    int sourceIndex = -1;

    QString source = columns[1];

    if(!source.compare("MC", Qt::CaseInsensitive))
    {
      sourceIndex = 1;
      sourceOk = true;
    }
    else if(!source.compare("GR", Qt::CaseInsensitive))
    {
      sourceIndex = 2;
      sourceOk = true;
    }

    bool variableOk = false;
    int variableIndex;
    QString variable = columns[0];

    if(!QString::compare(variable,"Temperature", Qt::CaseInsensitive))
    {
      variableIndex = -1;
      variableOk = true;
    }
    else
    {
      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        std::string solute = m_solutes[i];

        if (!solute.compare(variable.toStdString()))
        {
          variableIndex = i;
          variableOk = true;
        }
      }
    }

    if(variableOk && sourceOk)
    {
      QString filePath = columns[2];

      if(!filePath.isEmpty() && !filePath.isNull())
      {
        QFileInfo fileInfo(filePath);

        if (fileInfo.isRelative())
          fileInfo = relativePathToAbsolute(fileInfo);

        if (QFile::exists(fileInfo.absoluteFilePath()))
        {
          std::map<double, std::vector<double>> timeSeries;
          std::vector<std::string> headers;

          if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
          {
            for (size_t i = 0; i < headers.size(); i++)
            {
              auto eit = m_elementsById.find(headers[i]);

              if (eit != m_elementsById.end())
              {
                Element *element = eit->second;
                BoundaryCondition *boundaryCondition = new BoundaryCondition(element, sourceIndex, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[i];
                  boundaryCondition->addValue(dateTime, value);
                }


                m_boundaryConditions.push_back(boundaryCondition);
              }
              else
              {
                errorMessage = "Specified main channel boundary condition file is invalid";
                return false;
              }
            }
          }
          else
          {
            errorMessage = "Specified main channel boundary condition file is invalid";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified main channel boundary condition file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified main channel boundary condition file is invalid";
        return false;
      }
    }
  }
  else
  {
    errorMessage = "Specified main channel boundary condition file is invalid";
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
  if (m_outputCSVStream.device()->isOpen())
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

#ifdef USE_OPENMP
#pragma omp critical
#endif
  {
    if(m_outputNetCDF)
    {
      size_t currentTime = m_outputNetCDF->getDimSize("time");

      //Set current dateTime
      //            NcVar timeVar =  m_outputNetCDF->ncFile()->getVar("time");
      //            timeVar.putVar(std::vector<size_t>({currentTime}), m_currentDateTime);
      m_outputNetCDF->putVar("time", std::vector<size_t>({currentTime}), m_currentDateTime);

      double *depth = new double[m_elements.size()];
      double *width = new double[m_elements.size()];
      double *xsectArea = new double[m_elements.size()];
      double *temperature = new double[m_elements.size()];
      double *totalElementHeatBalance = new double[m_elements.size()];
      double *totalElementRadiationFluxHeatBalance = new double[m_elements.size()];
      double *totalElementExternalHeatFluxBalance = new double[m_elements.size()];
      double *channelConductionFlux = new double[m_elements.size()];
      double *groundConductionFlux = new double[m_elements.size()];
      double *channelAdvectionHeatFlux = new double[m_elements.size()];
      double *solutes = new double[m_elements.size() * m_solutes.size()];
      double *totalElementSoluteMassBalance = new double[m_elements.size() * m_solutes.size()];
      double *totalElementExternalSoluteFluxMassBalance = new double[m_elements.size() * m_solutes.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (size_t i = 0; i < m_elements.size(); i++)
      {
        Element *element = m_elements[i];
        depth[i] = element->depth;
        width[i] = element->width;
        xsectArea[i] = element->xSectionArea;
        temperature[i] = element->temperature.value;

        channelConductionFlux[i] = element->mainChannelConductionHeat / (element->width * element->length);
        groundConductionFlux[i] = element->groundConductionHeat / (element->width * element->length);
        channelAdvectionHeatFlux[i] = element->advectionHeat / (element->width * element->length);

        totalElementHeatBalance[i] = element->totalHeatBalance;
        totalElementRadiationFluxHeatBalance[i] = element->totalRadiationFluxesHeatBalance;
        totalElementExternalHeatFluxBalance[i] = element->totalExternalHeatFluxesBalance;

        for (size_t j = 0; j < m_solutes.size(); j++)
        {
          solutes[j * m_elements.size() + i] = element->soluteConcs[j].value;
          totalElementSoluteMassBalance[j * m_elements.size() + i] = element->totalSoluteMassBalance[j];
          totalElementExternalSoluteFluxMassBalance[j * m_elements.size() + i] = element->totalExternalSoluteFluxesMassBalance[j];
        }
      }

      //      NcVar depthVar =  m_outputNetCDF->ncFile()->getVar("depth");
      //      NcVar widthVar =  m_outputNetCDF->ncFile()->getVar("width");
      //      NcVar xsectAreaVar =  m_outputNetCDF->ncFile()->getVar("xsection_area");
      //      NcVar temperatureVar =  m_outputNetCDF->ncFile()->getVar("temperature");
      //      NcVar totalHeatBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_heat_balance");
      //      NcVar totalRadiationFluxHeatBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_radiation_flux_heat_balance");
      //      NcVar totalExternalHeatFluxBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_external_heat_flux_balance");
      //      NcVar totalElementHeatBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_element_heat_balance");
      //      NcVar totalElementRadiationFluxHeatBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_element_radiation_flux_heat_balance");
      //      NcVar totalElementExternalHeatFluxBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_element_external_heat_flux_balance");
      //      NcVar solutesVar =  m_outputNetCDF->ncFile()->getVar("solute_concentration");
      //      NcVar totalSoluteMassBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_solute_mass_balance");
      //      NcVar totalExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_external_solute_flux_mass_balance");
      //      NcVar totalElementSoluteMassBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_element_solute_mass_balance");
      //      NcVar totalElementExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->ncFile()->getVar("total_element_external_solute_flux_mass_balance");
      //      NcVar channelConductionFluxVar =  m_outputNetCDF->ncFile()->getVar("channel_conduction_heat_flux");
      //      NcVar groundConductionFluxVar =  m_outputNetCDF->ncFile()->getVar("ground_conduction_heat_flux");
      //      NcVar channelAdvectionHeatFluxVar =  m_outputNetCDF->ncFile()->getVar("channel_advection_heat_flux");

      //      depthVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), depth);
      m_outputNetCDF->putVar("depth", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), depth);

      //      widthVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), width);
      m_outputNetCDF->putVar("width", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), width);

      //      xsectAreaVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), xsectArea);
      m_outputNetCDF->putVar("xsection_area", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), xsectArea);

      //      temperatureVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), temperature);
      m_outputNetCDF->putVar("temperature", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), temperature);

      //      channelConductionFluxVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), channelConductionFlux);
      m_outputNetCDF->putVar("channel_conduction_heat_flux", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), channelConductionFlux);

      //      groundConductionFluxVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), groundConductionFlux);
      m_outputNetCDF->putVar("ground_conduction_heat_flux", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), groundConductionFlux);

      //      channelAdvectionHeatFluxVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), channelAdvectionHeatFlux);
      m_outputNetCDF->putVar("channel_advection_heat_flux", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), channelAdvectionHeatFlux);

      //      totalHeatBalanceVar.putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalHeatBalance);
      m_outputNetCDF->putVar("total_heat_balance", std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalHeatBalance);

      //      totalRadiationFluxHeatBalanceVar.putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalRadiationHeatBalance);
      m_outputNetCDF->putVar("total_radiation_flux_heat_balance", std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalRadiationHeatBalance);

      //      totalExternalHeatFluxBalanceVar.putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalExternalHeatFluxBalance);
      m_outputNetCDF->putVar("total_external_heat_flux_balance", std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalExternalHeatFluxBalance);

      //      totalElementHeatBalanceVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementHeatBalance);
      m_outputNetCDF->putVar("total_element_heat_balance", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementHeatBalance);

      //      totalElementRadiationFluxHeatBalanceVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementRadiationFluxHeatBalance);
      m_outputNetCDF->putVar("total_element_radiation_flux_heat_balance", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementRadiationFluxHeatBalance);

      //      totalElementExternalHeatFluxBalanceVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementExternalHeatFluxBalance);
      m_outputNetCDF->putVar("total_element_external_heat_flux_balance", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementExternalHeatFluxBalance);

      //      solutesVar.putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), solutes);
      m_outputNetCDF->putVar("solute_concentration", std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), solutes);

      //      totalSoluteMassBalanceVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalSoluteMassBalance.data());
      m_outputNetCDF->putVar("total_solute_mass_balance", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalSoluteMassBalance.data());

      //      totalExternalSoluteFluxMassBalanceVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalExternalSoluteFluxMassBalance.data());
      m_outputNetCDF->putVar("total_external_solute_flux_mass_balance", std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalExternalSoluteFluxMassBalance.data());

      //      totalElementSoluteMassBalanceVar.putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), totalElementSoluteMassBalance);
      m_outputNetCDF->putVar("total_element_solute_mass_balance", std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), totalElementSoluteMassBalance);

      //      totalElementExternalSoluteFluxMassBalanceVar.putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), totalElementExternalSoluteFluxMassBalance);
      m_outputNetCDF->putVar("total_element_external_solute_flux_mass_balance", std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), totalElementExternalSoluteFluxMassBalance);

      delete[] depth;
      delete[] width;
      delete[] xsectArea;
      delete[] temperature;
      delete[] totalElementHeatBalance;
      delete[] totalElementRadiationFluxHeatBalance;
      delete[] totalElementExternalHeatFluxBalance;
      delete[] solutes;
      delete[] totalElementSoluteMassBalance;
      delete[] totalElementExternalSoluteFluxMassBalance;
      delete[] channelConductionFlux;
      delete[] groundConductionFlux;
      delete[] channelAdvectionHeatFlux;

      if (m_flushToDisk)
      {
        m_outputNetCDF->sync();
      }
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
                                                              {"[POINT_SOURCES]", 6},
                                                              {"[NON_POINT_SOURCES]", 7},
                                                              {"[UNIFORM_HYDRAULICS]", 8},
                                                              {"[NON_UNIFORM_HYDRAULICS]", 9},
                                                              {"[UNIFORM_RADIATIVE_FLUXES]", 10},
                                                              {"[NON_UNIFORM_RADIATIVE_FLUXES]", 11},
                                                              {"[UNIFORM_BOUNDARY_CONDITION]", 12},
                                                              {"[NON_UNIFORM_BOUNDARY_CONDITION]", 13}
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
