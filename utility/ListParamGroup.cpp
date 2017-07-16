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

#include "ListParamGroup.h"
#include <QMessageBox>

ListParamGroup::ListParamGroup() {

}

ListParamGroup::ListParamGroup(QString n, QVector<SingleParam*> params) :
		ParamGroup::ParamGroup(n, params), addbtn(0), rmvbtn(0) {
	for (QVector<SingleParam*>::iterator it = parameters.begin() ; it != parameters.end() ; ++it)
		if ((*it)->IsList())
			listparams.push_back(dynamic_cast<ListValueParam*>(*it));
}

ListParamGroup::~ListParamGroup() {
	if (addbtn)
		delete addbtn;
	if (rmvbtn)
		 delete rmvbtn;
}

bool ListParamGroup::ParseParamValues(QString str) {
	bool ok = ParamGroup::ParseParamValues(str);
	// Check that list values has no doublons
	int minNbRow = 999999;
	for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) {
		minNbRow = (minNbRow > (*it)->GetNbRows()) ? (*it)->GetNbRows() : minNbRow;
	}
	QVector<QString> lines;
	QString strTemp = "";
	for (int row = 0 ; row < minNbRow ; ++row) {
		strTemp = "";
		for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) {
			strTemp = strTemp + " " + (*it)->GetValue(row);
		}
		if (lines.indexOf(strTemp) != -1) {
			for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) {
				(*it)->ChangeRow(row);
				(*it)->rmvItem();
			}
			row--;
			minNbRow--;
		} else {
			lines.push_back(strTemp);
		}
	}
	return ok;
}

bool ListParamGroup::AddToLayout(QLayout *layout) {
	if (not HLay)
		HLay = new QHBoxLayout();
	if (not lbl)
		lbl = new QLabel(name);
	HLay->addWidget(lbl);

	ListValueParam *tmp = 0;
	for (QVector<SingleParam*>::iterator it = parameters.begin() ; it != parameters.end() ; ++it) {
		HLay->addWidget((*it)->GetLayoutItem());
		tmp = dynamic_cast<ListValueParam*>(*it);
		if (tmp)
			QObject::connect(*it, SIGNAL(rowChanged(int)), this, SLOT(currRowChanged(int)));
	}

	addbtn = new QPushButton("+");
	rmvbtn = new QPushButton("-");

	HLay->addWidget(addbtn);
	HLay->addWidget(rmvbtn);

	QObject::connect(addbtn, SIGNAL(clicked()), this, SLOT(addItems()));
	QObject::connect(rmvbtn, SIGNAL(clicked()), this, SLOT(rmvItems()));


	layout->addItem(HLay);
	return true;
}

void ListParamGroup::addItems() {
	bool allFieldsFull = true;
	for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) {
		allFieldsFull &= not (*it)->IsFieldEmpty();
	}
	if (allFieldsFull) {
		for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) {
			(*it)->addItem();
		}
	} else {
		QMessageBox::about(0, "Empty Field", "One of the fields is empty or contains a space.");
	}
}

void ListParamGroup::rmvItems() {
	for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) {
		(*it)->rmvItem();
	}
}


void ListParamGroup::currRowChanged(int row) {
	for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) {
		(*it)->ChangeRow(row);
	}
}

bool ListParamGroup::CheckModifiedState() {
	bool tempMod = ParamGroup::CheckModifiedState();
	if (not tempMod) {
		for (QVector<ListValueParam*>::iterator it = listparams.begin() ; it != listparams.end() ; ++it) 
			tempMod |= (*it)->IsParamModified();
	}
	setLblBold(tempMod);
	return tempMod;
}

QString ListParamGroup::SaveToString() const {
	QString res;
	QString baseRes = ParamGroup::SaveToString();
	if (not listparams.empty()) {
		unsigned int nbLines = listparams[0]->GetNbRows();
		for (unsigned int i = 0 ; i < nbLines ; ++i) {
			res += baseRes;
			for (QVector<ListValueParam*>::const_iterator it = listparams.begin() ; it != listparams.end() ; ++it) 
				if (*it)
					res += (*it)->GetValue(i) + " ";
		}
	}
	return res;
}

