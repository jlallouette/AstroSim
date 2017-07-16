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

#include "SingleValueParam.h"

#include <QMessageBox>

SingleValueParam::SingleValueParam() :
	SingleParam::SingleParam(), chkbox(0), lineedt(0), combBox(0) {
}

SingleValueParam::~SingleValueParam() {
	if (chkbox)
		delete chkbox;
	if (lineedt)
		delete lineedt;
	if (combBox)
		delete combBox;
}

bool SingleValueParam::LoadFromString(QString str) {
	return SingleParam::LoadFromString(str);
}

bool SingleValueParam::ParseValue(QString paramVal) {
	if (IsBoolean()) {
		if (chkbox)
			chkbox->setChecked(true);
		return true;
	} else {
		if (not allowedVals.isEmpty()) {
			if (allowedVals.indexOf(paramVal) != -1) {
				if (combBox)
					combBox->setCurrentIndex(combBox->findText(paramVal));
			} else {
				QMessageBox::about(0, "Allowed value error", paramVal + 
					" is not an allowed value.");
				return false;
			}
		} else {
			if (lineedt)
				lineedt->setText(paramVal);
		}
		return true;
	}
}

QWidget* SingleValueParam::GetLayoutItem() {
	// If value is boolean, make a checkbox
	if (IsBoolean()) {
		if (not chkbox)
			chkbox = new QCheckBox();
		chkbox->setChecked(defaultVal == "1");
		return chkbox;
	} else {
		// If there are restrained allowed values make a combobox
		if (not allowedVals.isEmpty()) {
			if (not combBox)
				combBox = new QComboBox();
			for (QVector<QString>::iterator it = allowedVals.begin() ; it != allowedVals.end() ; ++it)
				combBox->addItem(*it);
			if (not defaultVal.isEmpty())
				combBox->setCurrentIndex(combBox->findText(defaultVal));
			return combBox;
		} else { // else, if value is unconstrained, make a line edit
			if (not lineedt)
				lineedt = new QLineEdit();
			lineedt->setText(defaultVal);
			return lineedt;
		}
	}
}

bool SingleValueParam::IsParamModified() const {
	if (IsBoolean())
		return chkbox ? not ((chkbox->isChecked() and (defaultVal == "1")) or (not chkbox->isChecked() and (defaultVal != "1"))) : false;
	else if (not allowedVals.isEmpty())
		return combBox ? (combBox->currentText() != defaultVal) : false;
	else
		return lineedt ? (lineedt->text() != defaultVal) : false;
}

QString SingleValueParam::SaveToString() const {
	if (IsBoolean())
		return " ";
	else if (not allowedVals.isEmpty())
		return combBox ? (combBox->currentText() + " ") : " ";
	else
		return lineedt ? (lineedt->text() + " ") : " ";
}
