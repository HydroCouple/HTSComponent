#include "stdafx.h"
#include "test/htscomponenttest.h"
#include "htsmodel.h"


void HTSComponentTest::kuparuk_river1()
{
  QBENCHMARK_ONCE
  {

    //Error messages
    std::list<std::string> errors;

    //Create model instance
    HTSModel *model = new HTSModel(nullptr);

    //Read input filel
    model->setInputFile(QFileInfo("../../examples/kuparuk_river1/kuparuk_river1.inp"));

    //initialize model
    if(model->initialize(errors))
    {

      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //finalize model
    model->finalize(errors);

    delete model;
  }
}
