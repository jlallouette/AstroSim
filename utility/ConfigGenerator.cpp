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

#include "ConfigGenerator.h"

#include <QFileDialog>
#include <QGridLayout>

ConfigGenerator::ConfigGenerator(QWidget *parent)
    : QMainWindow(parent), paramDescriptionFilePath("./Params.ini")
{
	ui.setupUi(this);
}

ConfigGenerator::~ConfigGenerator()
{

}

void ConfigGenerator::menuBarAction(QAction *action)
{
	if (action == ui.actionOpen_Parameters_Description_File)
	{
		paramDescriptionFilePath = QFileDialog::getOpenFileName(this, tr("Open Parameters Description File"));
		descript.LoadFromFile(paramDescriptionFilePath);
		QGridLayout *layout = new QGridLayout();
		descript.SetFullLayout(layout);
		if (ui.scrollAreaWidgetContents->layout())
			delete ui.scrollAreaWidgetContents->layout();
		ui.scrollAreaWidgetContents->setLayout(layout);
	} else if (action == ui.actionLoad_Parameters) {
		QString paramFilePath = QFileDialog::getOpenFileName(this, tr("Open Parameters Description File"));
		descript.LoadParametersFromFile(paramFilePath);
	} else if (action == ui.actionSave_Parameters) {
		QString paramSaveFilePath = QFileDialog::getSaveFileName(this, tr("Save Parameters To File"));
		descript.SaveParametersToFile(paramSaveFilePath);
	}

}
