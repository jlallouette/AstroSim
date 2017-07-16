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

#ifndef PARAMDESCRIPTION_H_
#define PARAMDESCRIPTION_H_

#include <QString>
#include <QMap>
#include <QLayout>
#include <QTimer>

#include "ParamGroup.h"

class ParamDescription : public QObject {
	Q_OBJECT

public:
	ParamDescription();
	virtual ~ParamDescription();

	bool LoadFromFile(QString path);
	bool LoadParametersFromFile(QString path);
	bool SaveParametersToFile(QString path);
	bool SetFullLayout(QLayout *layout);

	void CleanParamGroups();

protected:
	QMap<QString, ParamGroup*> parameterGroups;
	QTimer *checkModTimer;
	QLayout *totLayout;

	bool parseParamLine(QString line);

private slots :
	void checkModState();
};

#endif /* PARAMDESCRIPTION_H_ */
