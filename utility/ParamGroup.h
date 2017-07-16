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

#ifndef PARAMGROUP_H_
#define PARAMGROUP_H_

#include <QLabel>
#include <QVector>
#include <QHBoxLayout>
#include "SingleParam.h"

class ParamGroup : public QObject {
public:
	ParamGroup();
	ParamGroup(QString n, QVector<SingleParam*> params);
	ParamGroup(QString n);
	virtual ~ParamGroup();

	virtual bool LoadFromString(QString str);
	virtual bool ParseParamValues(QString str);
	virtual bool AddToLayout(QLayout *layout);
	virtual QString SaveToString() const;

	// Checks wether the parameter group has some values different
	// from default values
	virtual bool CheckModifiedState();

	QString GetName() const;
	QLayout* GetLayout() const
	{ return HLay; }

	static ParamGroup* CreateParamGroupFromString(QString str);

protected:
	QString name;
	QVector<SingleParam*> parameters;

	QHBoxLayout *HLay;
	QLabel *lbl;

	void setLblBold(bool b);
};

#endif /* PARAMGROUP_H_ */
