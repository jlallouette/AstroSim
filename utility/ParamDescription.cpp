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

#include "ParamDescription.h"
#include <QFile>
#include <QTextStream>
#include <QStringList>

#include <QMessageBox>

ParamDescription::ParamDescription() : checkModTimer(0) {

}

ParamDescription::~ParamDescription() {
	CleanParamGroups();
	if (checkModTimer) {
		checkModTimer->stop();
		delete checkModTimer;
	}
}

void ParamDescription::CleanParamGroups() {
	for (QMap<QString, ParamGroup*>::iterator it = parameterGroups.begin() ; it != parameterGroups.end() ; ++it)
		if (it.value())
			delete it.value();

	parameterGroups.clear();
}

bool ParamDescription::LoadFromFile(QString path) {
	CleanParamGroups();

	QFile file(path);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&file);
	ParamGroup *pg = 0;
	// First line is discarded
	in.readLine();
	bool ok = true;
	while (!in.atEnd()) {
		QString line = in.readLine();
		pg = ParamGroup::CreateParamGroupFromString(line);
		if (pg)
			parameterGroups[pg->GetName()] = pg;
		else
			ok = false;
	}
	return ok;
}

bool ParamDescription::SaveParametersToFile(QString path) {
	QFile file(path);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		return false;

	QString paramsTot;
	for (QMap<QString, ParamGroup*>::iterator it = parameterGroups.begin() ; it != parameterGroups.end() ; ++it)
		if (*it and (*it)->CheckModifiedState())
			paramsTot += (*it)->SaveToString();

	QTextStream out(&file);
	out << paramsTot;

	file.close();
	return true;
}

bool ParamDescription::LoadParametersFromFile(QString path) {
	QFile file(path);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&file);
	bool ok = true;
	while (!in.atEnd()) {
		QString line = in.readLine();
		ok &= parseParamLine(line);
	}
	checkModState();
	return ok;
}

bool ParamDescription::SetFullLayout(QLayout *layout) {
	bool ok = true;
	totLayout = layout;
	for (QMap<QString, ParamGroup*>::iterator it = parameterGroups.begin() ; it != parameterGroups.end() ; ++it)
		ok &= it.value()->AddToLayout(layout);

	if (not checkModTimer) {
		checkModTimer = new QTimer();
		checkModTimer->setInterval(500);
		QObject::connect(checkModTimer, SIGNAL(timeout()), this, SLOT(checkModState()));
		checkModTimer->start();
	}
	return ok;
}

bool ParamDescription::parseParamLine(QString line) {
	line = line.trimmed();
	bool ok = true;

	QStringList splitRes = line.split(" ", QString::SkipEmptyParts);
	QString paramName = splitRes.at(0).trimmed();
	QString tempValues = "";
	bool firstTime = true;
	for (QStringList::iterator it = splitRes.begin() ; it != splitRes.end() ; ++it) {
		if (parameterGroups.find(it->trimmed()) != parameterGroups.end()) {
			if (firstTime)
				firstTime = false;
			else
				ok &= parameterGroups[paramName]->ParseParamValues(tempValues.trimmed());
			paramName = it->trimmed();
			tempValues = "";
		} else
			tempValues += " " + *it;
	}
	ok &= parameterGroups[paramName]->ParseParamValues(tempValues.trimmed());

	return ok;
}

void ParamDescription::checkModState() {
	for (QMap<QString, ParamGroup*>::iterator it = parameterGroups.begin() ; it != parameterGroups.end() ; ++it)
	{
		if (not it.value()->CheckModifiedState())
		{
			QLayout* l = it.value()->GetLayout();
			totLayout->removeItem(l);
			totLayout->addItem(l);
		}
	}
}

