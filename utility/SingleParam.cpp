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

#include "SingleParam.h"
#include <QStringList>

#include <QListWidget>
#include <QComboBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>

#include "SingleValueParam.h"
#include "ListValueParam.h"

SingleParam::SingleParam() {
}

SingleParam::SingleParam(QString t, QString defVal, QVector<QString> allowed) :
	type(t), defaultVal(defVal), allowedVals(allowed) {

}

SingleParam::~SingleParam() {
}

SingleParam* SingleParam::CreateSingleParamFromString(QString str) {
	str = str.trimmed();
	QString type = str.left(str.indexOf(' ')).trimmed();
	SingleParam *result = 0;
	if (type.contains("vector")) {
		result = new ListValueParam();
	} else {
		result = new SingleValueParam();
	}
	if (not result->LoadFromString(str)) {
		delete result;
		result = 0;
	}
	return result;
}

bool SingleParam::LoadFromString(QString str) {
	int posTypeEnd = str.indexOf(' ');
	type = str.left(posTypeEnd).trimmed();
	str = str.mid(posTypeEnd+1);
	int posAllowStart = str.indexOf('[');
	if (posAllowStart != -1) {
		str = str.mid(posAllowStart+1);
		int posAllowEnd = str.indexOf(']');
		if (posAllowEnd == -1)
			return false;
		QStringList allow = str.left(posAllowEnd).split(' ');
		for (QStringList::iterator it = allow.begin() ; it != allow.end() ; ++it)
			allowedVals.push_back(it->trimmed());
		str = str.mid(posAllowEnd+1);
	}
	int posDefaultStart = str.indexOf('(');
	if (posDefaultStart != -1) {
		str = str.mid(posDefaultStart+1);
		int posDefaultEnd = str.indexOf(')');
		if (posDefaultEnd == -1)
					return false;
		defaultVal = str.left(posDefaultEnd).trimmed();
	}
	return true;
}

bool SingleParam::IsBoolean() const {
	return type == "b";
}

bool SingleParam::IsList() const {
	return type.contains("vector");
}
