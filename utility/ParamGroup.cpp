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

#include "ParamGroup.h"
#include "ListParamGroup.h"

#include <QStringList>
#include <QMessageBox>

ParamGroup::ParamGroup() : HLay(0), lbl(0) {
}

ParamGroup::ParamGroup(QString n) : name(n), HLay(0), lbl(0) {
}

ParamGroup::ParamGroup(QString n, QVector<SingleParam*> params)
	: name(n), parameters(params), HLay(0), lbl(0) {

}

ParamGroup::~ParamGroup() {
	if (lbl)
		delete lbl;
	for (QVector<SingleParam*>::iterator it = parameters.begin() ; it != parameters.end() ; ++it)
		if (*it)
			delete *it;
}

bool ParamGroup::LoadFromString(QString str) {
	bool ok = true;
	QStringList params = str.split(", ");
	for (QStringList::iterator it = params.begin() ; it != params.end() ; ++it) {
		parameters.push_back(SingleParam::CreateSingleParamFromString(*it));
		ok &= (parameters.back() != 0);
	}
	return ok;
}

bool ParamGroup::ParseParamValues(QString str) {
	bool ok = true;
	QStringList splitRes = str.split(" ", QString::SkipEmptyParts);
	if (splitRes.size() > 0) {
		if (splitRes.size() <= parameters.size()) {
			int paramNum = 0;
			if (parameters[0]->IsBoolean()) {
				ok &= parameters[0]->ParseValue("");
				paramNum = 1;
			}
			for (QStringList::iterator it = splitRes.begin() ; it != splitRes.end() ; ++it) {
				if (not it->trimmed().isEmpty())
					ok &= parameters[paramNum++]->ParseValue(it->trimmed());
			}
		} else {
			QMessageBox::about(0, "Too much params", "Too much parameter values for parameter " + name);
			ok = false;
		}
	} else {
		if (parameters.size() > 0)
			ok &= parameters[0]->ParseValue("");
		else
			ok = false;
	}

	return ok;
}

bool ParamGroup::AddToLayout(QLayout *layout) {
	if (not HLay)
		HLay = new QHBoxLayout();
	if (not lbl)
		lbl = new QLabel(name);
	HLay->addWidget(lbl);
	for (QVector<SingleParam*>::iterator it = parameters.begin() ; it != parameters.end() ; ++it)
		if (*it)
			HLay->addWidget((*it)->GetLayoutItem());
	layout->addItem(HLay);
	return true;
}

ParamGroup* ParamGroup::CreateParamGroupFromString(QString str) {
	QString paramName;

	QStringList splitRes = str.split(" : ");
	if (splitRes.size() != 2)
		QMessageBox::about(0, "Wrong format", "Wrong format for line : " + str);
	paramName = splitRes.at(0);

	QVector<SingleParam*> parameters;
	QStringList params = splitRes.at(1).split(", ");
	bool containsVects = false;
	for (QStringList::iterator it = params.begin() ; it != params.end() ; ++it) {
		parameters.push_back(SingleParam::CreateSingleParamFromString(*it));
		containsVects |= parameters.back()->IsList();
	}

	if (containsVects)
		return new ListParamGroup(paramName, parameters);
	else
		return new ParamGroup(paramName, parameters);
}

QString ParamGroup::GetName() const {
	return name;
}

bool ParamGroup::CheckModifiedState() {
	bool tempMod = false;
	for (QVector<SingleParam*>::iterator it = parameters.begin() ; it != parameters.end() ; ++it)
		if (*it)
			tempMod |= (*it)->IsParamModified();

	setLblBold(tempMod);

	return tempMod;
}

void ParamGroup::setLblBold(bool b) {
	if (lbl) {
		QFont currFont = lbl->font();
		currFont.setBold(b);
		currFont.setUnderline(b);
		lbl->setFont(currFont);
	}
}

QString ParamGroup::SaveToString() const {
	QString res;
	res += name + " ";
	QString uncertain;
	for (QVector<SingleParam*>::const_iterator it = parameters.begin() ; it != parameters.end() ; ++it)
		if (*it) {
			if ((*it)->IsParamModified()) {
				res += uncertain + (*it)->SaveToString() + " ";
				uncertain = "";
			} else
				uncertain += (*it)->SaveToString() + " ";
		}

	return res;
}

