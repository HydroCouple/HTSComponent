/*!
*  \file    radiativefluxtimeseriesbc.h
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

#ifndef RADIATIVEFLUXBC_H
#define RADIATIVEFLUXBC_H

#include "iboundarycondition.h"
#include "htscomponent_global.h"

#include <QObject>
#include <QSharedPointer>
#include <unordered_map>

struct Element;
class TimeSeries;
class HTSModel;
class DataCursor;

class HTSCOMPONENT_EXPORT RadiativeFluxBC: public QObject,
    public virtual IBoundaryCondition
{
    Q_OBJECT

  public:

    RadiativeFluxBC(Element *startElement,
                    Element *endElement,
                    HTSModel *model);

    virtual ~RadiativeFluxBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

    QSharedPointer<TimeSeries> timeSeries() const;

    void setTimeSeries(const QSharedPointer<TimeSeries> &timeseries);

  private:

    std::vector<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    DataCursor *m_dataCursor;
    QSharedPointer<TimeSeries> m_timeSeries;
    HTSModel *m_model;
};



#endif // RADIATIVEFLUXBC_H
