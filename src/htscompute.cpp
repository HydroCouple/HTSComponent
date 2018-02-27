#include "htsmodel.h"
#include "element.h"
#include "elementjunction.h"

void HTSModel::update()
{
  if(m_currentDateTime < m_endDateTime)
  {
    prepareForNextTimeStep();

    m_timeStep = computeTimeStep();

    applyBoundaryConditions(m_currentDateTime + m_timeStep / 86400.0);

    //Solve the transport for each element
    {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {
#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          //Solve heat transport first
          solveHeatTransport(m_timeStep);
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(size_t i = 0 ; i < m_solutes.size(); i++)
          {
            solveSoluteTransport(i, m_timeStep);
          }
        }
      }
    }


    m_currentDateTime +=  m_timeStep / 86400.0;

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime += m_outputInterval / 86400.0;
    }
  }
}

void HTSModel::prepareForNextTimeStep()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    {
      element->prevTemperature.copy(element->temperature);
    }

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      {
        element->prevSoluteConcs[j].copy(element->soluteConcs[j]);
      }
    }
  }

}

void HTSModel::applyInitialConditions()
{


  //Write initial output
  writeOutput();

  //Set next output time
  m_nextOutputTime += m_outputInterval / 86400.0;
}

void HTSModel::applyBoundaryConditions(double dateTime)
{

}

double HTSModel::computeTimeStep()
{
  double timeStep = m_maxTimeStep;

  double maxCourantFactor = 0.0;// Δx / v (s^-1)

  if(m_useAdaptiveTimeStep)
  {


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0 ; i < m_elements.size()  ; i++)
    {
      Element *element = m_elements[i];
      double dispersionFactor = element->computeDispersionFactor();


      if(dispersionFactor > maxCourantFactor)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        maxCourantFactor = dispersionFactor;
      }
    }

    timeStep = maxCourantFactor ? m_timeStepRelaxationFactor / maxCourantFactor : m_maxTimeStep;\

    if(m_numCurrentInitFixedTimeSteps < m_numInitFixedTimeSteps)
    {
      timeStep = std::min(timeStep, m_minTimeStep);
      m_numCurrentInitFixedTimeSteps++;
    }
  }

  double nextTime = m_currentDateTime + timeStep / 86400.0;

  if(nextTime > m_nextOutputTime)
  {
    timeStep = std::max(m_minTimeStep,  (m_nextOutputTime - m_currentDateTime) *  86400.0);
  }

  timeStep = std::min(std::max(timeStep, m_minTimeStep), m_maxTimeStep);

  return timeStep;
}


void HTSModel::solveHeatTransport(double timeStep)
{
  //Allocate memory to store inputs and outputs
  double *currentTemperatures = new double[m_elements.size()];
  double *outputTemperatures = new double[m_elements.size()];

  //Set initial input and output values to current values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentTemperatures[element->index] = element->temperature.value;
    outputTemperatures[element->index] = element->temperature.value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;
  m_solver->solve(currentTemperatures, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                  outputTemperatures, &HTSModel::computeDTDt, &solverUserData);

  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    double outputTemperature = outputTemperatures[element->index];
    element->temperature.value = outputTemperature;
  }

  //Delete allocated memory
  delete[] currentTemperatures;
  delete[] outputTemperatures;
}

void HTSModel::solveSoluteTransport(int soluteIndex, double timeStep)
{
  double *currentSoluteConcs = new double[m_elements.size()];
  double *outputSoluteConcs = new double[m_elements.size()];

  //Set initial values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
    outputSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = soluteIndex;
  m_solver->solve(outputSoluteConcs, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                  outputSoluteConcs, &HTSModel::computeDSoluteDt, &solverUserData);


  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->soluteConcs[soluteIndex].value = outputSoluteConcs[element->index];
  }


  //Delete allocated memory
  delete[] currentSoluteConcs;
  delete[] outputSoluteConcs;

}

void HTSModel::computeDTDt(double t, double y[], double dydt[], void* userData)
{

  SolverUserData *solverUserData = (SolverUserData*) userData;
  HTSModel *modelInstance = solverUserData->model;
  double dt = t - (modelInstance->m_currentDateTime *  86400.0);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    double DTDt = element->computeDTDt(dt,y);
    dydt[element->index] = DTDt;
  }
}

void HTSModel::computeDSoluteDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  HTSModel *modelInstance = solverUserData->model;
  double dt = t - modelInstance->m_currentDateTime *  86400.0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    double DSoluteDt = element->computeDSoluteDt(dt,y,solverUserData->variableIndex);
    dydt[element->index] = DSoluteDt;
  }
}
