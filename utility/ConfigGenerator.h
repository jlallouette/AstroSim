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

#ifndef CONFIGGENERATOR_H
#define CONFIGGENERATOR_H

#include <QtWidgets/QMainWindow>
#include "ui_ConfigGenerator.h"
#include <QAction>

#include "ParamDescription.h"

class ConfigGenerator : public QMainWindow
{
    Q_OBJECT

public:
    ConfigGenerator(QWidget *parent = 0);
    ~ConfigGenerator();

private:
    Ui_ConfigGeneratorClass ui;
    QString paramDescriptionFilePath;

    ParamDescription descript;

private slots :
	void menuBarAction(QAction *action);
};

#endif // SCRIPTGENERATOR_H
