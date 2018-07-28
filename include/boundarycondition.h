/*!
*  \file    boundarycondition.h
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
*  \todo Test transport on branching networks
*  \warning
*/

#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "abstracttimeseriesbc.h"

class HTSCOMPONENT_EXPORT BoundaryCondition : public AbstractTimeSeriesBC
{
  public:

    BoundaryCondition(Element *element, int sourceIndex, int variableIndex, HTSModel *model);

    virtual ~BoundaryCondition();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *element() const;

    void setElement(Element *element);

  private:

    Element *m_element;
    int m_sourceIndex;
    int m_variableIndex;
};

class HTSCOMPONENT_EXPORT UniformBoundaryCondition : public AbstractTimeSeriesBC
{
  public:

    UniformBoundaryCondition(Element *startElement, Element *endElement, int sourceIndex, int variableIndex, HTSModel *model);

    virtual ~UniformBoundaryCondition();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

  private:
    std::list<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    int m_sourceIndex;
    int m_variableIndex;
};

#endif // BOUNDARYCONDITION_H
