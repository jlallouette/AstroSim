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

#include "ListValueParam.h"

#include <QMessageBox>

ListValueParam::ListValueParam() :
	lineedt(0), combBox(0),	lst(0), Vlay(0), fullWidg(0) {
}

ListValueParam::~ListValueParam() {
	if (lineedt)
		delete lineedt;
	if (combBox)
		delete lineedt;
	if (lst)
		delete lst;
	if (Vlay)
		delete Vlay;
	if (fullWidg)
		delete fullWidg;
}

bool ListValueParam::LoadFromString(QString str) {
	return SingleParam::LoadFromString(str);
}

bool ListValueParam::ParseValue(QString paramVal) {
	if (not allowedVals.isEmpty()) {
		if (allowedVals.indexOf(paramVal) != -1) {
			if (lst)
				lst->addItem(paramVal);
		} else {
			QMessageBox::about(0, "Allowed value error", paramVal + " is not an allowed value");
			return false;
		}
	} else {
		if (lst)
			lst->addItem(paramVal);
	}
	return true;
}

QWidget* ListValueParam::GetLayoutItem() {
	fullWidg = new QWidget();
	lst = new QListWidget();
	lst->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	lst->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

	Vlay = new QVBoxLayout();

	if (not allowedVals.isEmpty()) {
		combBox = new QComboBox();
		for (QVector<QString>::iterator it = allowedVals.begin() ; it != allowedVals.end() ; ++it)
			combBox->addItem(*it);
		Vlay->addWidget(combBox);
	} else {
		lineedt = new QLineEdit();
		Vlay->addWidget(lineedt);
	}
	Vlay->addWidget(lst);

	QObject::connect(lst, SIGNAL(itemClicked(QListWidgetItem *)), this, SLOT(currRowChanged(QListWidgetItem *)));

	fullWidg->setLayout(Vlay);

	return fullWidg;
}

void ListValueParam::addItem() {
	if (lineedt)
		lst->addItem(lineedt->text());
	else if (combBox)
		lst->addItem(combBox->currentText());
}

void ListValueParam::rmvItem() {
	QListWidgetItem *tmp = lst->takeItem(lst->currentRow());
	if (tmp)
		delete tmp;
}

void ListValueParam::currRowChanged(QListWidgetItem * ) {
	emit rowChanged(lst->currentRow());
}

void ListValueParam::ChangeRow(int row) {
	if (lst)
		lst->setCurrentRow(row);
}

bool ListValueParam::IsFieldEmpty() const {
	return lineedt and (lineedt->text().trimmed().isEmpty() or lineedt->text().contains(" "));
}

int ListValueParam::GetNbRows() const {
	return lst ? lst->count() : 0;
}

QString ListValueParam::GetValue(int row) const {
	return lst ? lst->item(row)->text() : "";
}

bool ListValueParam::IsParamModified() const {
	return lst ? (lst->count() > 0) : false;
}

QString ListValueParam::SaveToString() const {
	return "";
}

