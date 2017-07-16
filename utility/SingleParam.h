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

#ifndef SINGLEPARAM_H_
#define SINGLEPARAM_H_

#include <QString>
#include <QVector>
#include <QWidget>

class SingleParam : public QObject {
public:
	SingleParam();
	SingleParam(QString t, QString defVal, QVector<QString> allowed);
	virtual ~SingleParam();

	virtual bool LoadFromString(QString str);
	virtual bool ParseValue(QString paramVal) = 0;
	virtual QWidget* GetLayoutItem() = 0;
	virtual QString SaveToString() const = 0;

	virtual bool IsParamModified() const = 0;

	bool IsBoolean() const;
	bool IsList() const;

	static SingleParam* CreateSingleParamFromString(QString str);

protected:
	QString name;
	QString type;
	QString defaultVal;
	QVector<QString> allowedVals;
};

#endif /* SINGLEPARAM_H_ */
