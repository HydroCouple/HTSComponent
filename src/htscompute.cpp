#include "htsmodel.h"
#include "element.h"
#include "elementjunction.h"
#include "iboundarycondition.h"
#include "htscomponent.h"

using namespace std;

void HTSModel::update()
{
  if(m_currentDateTime < m_endDateTime)
  {

    applyBoundaryConditions(m_currentDateTime);

    //Retrieve external data from other coupled models
    if(m_retrieveCouplingDataFunction)
    {
      (*m_retrieveCouplingDataFunction)(this, m_currentDateTime);
    }

    if(m_component)
      m_component->applyInputValues();

    m_timeStep = computeTimeStep();

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
          for(int i = 0 ; i < (int)m_solutes.size(); i++)
          {
            solveSoluteTransport(i, m_timeStep);
          }
        }
      }
    }

    m_prevDateTime = m_currentDateTime;
    m_currentDateTime = m_currentDateTime + m_timeStep / 86400.0;

    prepareForNextTimeStep();

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime = std::min(m_nextOutputTime + m_outputInterval / 86400.0 , m_endDateTime);
    }

    if(m_verbose)
    {
      printStatus();
    }
  }
}

void HTSModel::prepareForNextTimeStep()
{

  m_minTemp = std::numeric_limits<double>::max();
  m_maxTemp = std::numeric_limits<double>::lowest();

  std::fill(m_maxSolute.begin(), m_maxSolute.end(), m_maxTemp);
  std::fill(m_minSolute.begin(), m_minSolute.end(), m_minTemp);

  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computeHeatBalance(m_timeStep);
    m_totalHeatBalance += element->totalHeatBalance;
    m_totalRadiationHeatBalance += element->totalRadiationFluxesHeatBalance;
    m_totalExternalHeatFluxBalance += element->totalExternalHeatFluxesBalance;

    element->prevTemperature.copy(element->temperature);

    m_minTemp = min(m_minTemp , element->temperature.value);
    m_maxTemp = max(m_maxTemp , element->temperature.value);

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->computeSoluteBalance(m_timeStep, j);
      m_totalSoluteMassBalance[j] += element->totalSoluteMassBalance[j];
      m_totalExternalSoluteFluxMassBalance[j] += element->totalExternalSoluteFluxesMassBalance[j];

      element->prevSoluteConcs[j].copy(element->soluteConcs[j]);

      m_minSolute[j] = min(m_minSolute[j] , element->soluteConcs[j].value);
      m_maxSolute[j] = max(m_maxSolute[j] , element->soluteConcs[j].value);
    }
  }

}

void HTSModel::applyInitialConditions()
{

  //Initialize heat and solute balance trackers
  m_totalHeatBalance = 0.0;
  m_totalRadiationHeatBalance = 0.0;
  m_totalExternalHeatFluxBalance = 0.0;

  std::fill(m_totalSoluteMassBalance.begin(), m_totalSoluteMassBalance.end(), 0.0);
  std::fill(m_totalExternalSoluteFluxMassBalance.begin(), m_totalExternalSoluteFluxMassBalance.end(), 0.0);

  applyBoundaryConditions(m_currentDateTime);

  //Write initial output
  writeOutput();

  //Set next output time
  m_nextOutputTime += m_outputInterval / 86400.0;
}

void HTSModel::applyBoundaryConditions(double dateTime)
{
  //reset external fluxes
#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->externalHeatFluxes = 0.0;
    element->radiationFluxes = 0.0;

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->externalSoluteFluxes[j] = 0.0;
    }
  }

#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->applyBoundaryConditions(dateTime);
  }
}

double HTSModel::computeTimeStep()
{
  double timeStep = m_maxTimeStep;

  double maxCourantFactor = 0.0;// Δx / v (s^-1)

  if(m_numCurrentInitFixedTimeSteps < m_numInitFixedTimeSteps)
  {
    timeStep = m_minTimeStep;
    m_numCurrentInitFixedTimeSteps++;
  }
  else if(m_useAdaptiveTimeStep)
  {

    for(size_t i = 0 ; i < m_elements.size()  ; i++)
    {
      Element *element = m_elements[i];
      double dispersionFactor = element->computeDispersionFactor();


      if(dispersionFactor > maxCourantFactor)
      {
        maxCourantFactor = dispersionFactor;
      }
    }

    timeStep = maxCourantFactor ? m_timeStepRelaxationFactor / maxCourantFactor : m_maxTimeStep;

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
  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentTemperatures[element->index] = element->temperature.value;
    outputTemperatures[element->index] = element->temperature.value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;

  if(m_heatSolver->solve(currentTemperatures, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                  outputTemperatures, &HTSModel::computeDTDt, &solverUserData))
  {
    printf("HTS Temperature Solver failed \n");
  }

  //Apply computed values;
  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
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
  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
    outputSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = soluteIndex;

  if(m_soluteSolvers[soluteIndex]->solve(outputSoluteConcs, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                  outputSoluteConcs, &HTSModel::computeDSoluteDt, &solverUserData))
  {
    printf("HTS Solute Solver failed \n");
  }

  //Apply computed values;
  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
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
  for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
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
  for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    double DSoluteDt = element->computeDSoluteDt(dt,y,solverUserData->variableIndex);
    dydt[element->index] = DSoluteDt;
  }
}
