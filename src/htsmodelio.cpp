#include "htsmodel.h"
#include "element.h"
#include "elementjunction.h"

#include <QDir>

using namespace std;

QFileInfo HTSModel::discretizationFile() const
{
  return m_discretizationFile;
}

void HTSModel::setDiscretizationFile(const QFileInfo &discretizationFile)
{
  m_discretizationFile = discretizationFile;
}

QFileInfo HTSModel::hydrodynamicFile() const
{
  return m_hydrodynamicFile;
}

void HTSModel::setHydrodynamicFile(const QFileInfo &hydrodynamicFile)
{
  m_hydrodynamicFile = hydrodynamicFile;
}

QFileInfo HTSModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void HTSModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

bool HTSModel::initializeInputFiles(list<string> &errors)
{
  return true;
}

bool HTSModel::initializeOutputFiles(list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool HTSModel::initializeCSVOutputFile(list<string> &errors)
{
  QString file = m_outputCSVFileInfo.absoluteFilePath();

  if(!file.isEmpty() && !file.isNull() && !m_outputCSVFileInfo.dir().exists())
  {
    errors.push_back("Output shapefile directory does not exist: " + file.toStdString());
    return false;
  }

  if(!file.isEmpty() && !file.isNull())
  {
    if(m_outputCSVStream.device() == nullptr)
    {
      QFile *device = new QFile(file, this);
      m_outputCSVStream.setDevice( device);
    }

    if(m_outputCSVStream.device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_outputCSVStream.setRealNumberPrecision(10);
      m_outputCSVStream.setRealNumberNotation(QTextStream::SmartNotation);
      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Temperature";

      for(size_t i = 0; i < m_solutes.size(); i++)
      {
        m_outputCSVStream << ", " << QString::fromStdString(m_solutes[i]);
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
  return true;
}

void HTSModel::writeOutput()
{
  writeCSVOutput();
  writeNetCDFOutput();
}

void HTSModel::writeCSVOutput()
{
  if(m_outputCSVStream.device()->isOpen())
  {
    for(size_t i = 0 ; i < m_elements.size() ; i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->index
                        << ", " << element->x << ", " << element->y << ", " << element->z
                        <<  ", " << element->temperature.value;

      for(size_t j = 0 ; j < m_solutes.size(); j++)
      {
        m_outputCSVStream << ", " << element->soluteConcs[j].value;
      }

      m_outputCSVStream << endl;
    }
  }
}

void HTSModel::writeNetCDFOutput()
{

}

void HTSModel::closeOutputFiles()
{
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      closeCSVOutputFile();
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      closeOutputNetCDFFile();
    }
  }
}

void HTSModel::closeCSVOutputFile()
{
  if(m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }
}

void HTSModel::closeOutputNetCDFFile()
{

}
