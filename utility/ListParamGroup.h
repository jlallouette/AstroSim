/*----------------------------------------------------------------------------
  AstroSim: Simulation of astrocyte networks Ca2+ dynamics
  Copyright (c) 2016-2017 Jules Lallouette
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#ifndef LISTPARAMGROUP_H_
#define LISTPARAMGROUP_H_

#include "ParamGroup.h"

#include <QPushButton>
#include "ListValueParam.h"

class ListParamGroup: public ParamGroup {
	Q_OBJECT

public:
	ListParamGroup();
	ListParamGroup(QString n, QVector<SingleParam*> params);
	virtual ~ListParamGroup();

	virtual bool ParseParamValues(QString str);
	virtual bool AddToLayout(QLayout *layout);
	virtual QString SaveToString() const;

	virtual bool CheckModifiedState();

protected:
	QPushButton *addbtn;
	QPushButton *rmvbtn;

	QVector<ListValueParam*> listparams;

private slots :
	void addItems();
	void rmvItems();
	void currRowChanged(int row);
};

#endif /* LISTPARAMGROUP_H_ */
