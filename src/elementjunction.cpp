#include "stdafx.h"
#include "elementjunction.h"
#include "element.h"
#include "htsmodel.h"

#include <math.h>

ElementJunction::ElementJunction(const std::string &id, double x, double y, double z, HTSModel *model)
  :id(id), x(x), y(y), z(z),
    model(model)
{

}

ElementJunction::~ElementJunction()
{

  while (outgoingElements.size())
  {
    Element *element = *outgoingElements.begin();
    delete element;
  }

  while (incomingElements.size())
  {
    Element *element = *incomingElements.begin();
    delete element;
  }
}
