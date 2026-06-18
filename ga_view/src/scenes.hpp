#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// The ga_view scene catalog: construction of all demo models and their rendering into a
// QGraphicsScene. Split out of w_mainwindow.cpp so the main window stays pure window
// chrome. (Implementation in scenes.cpp.)

#include "coordsys_model.hpp"

#include <vector>

class Coordsys;
class w_Coordsys;
class QGraphicsScene;

// Build the full, ordered list of demo models/scenes shown in ga_view. Called once by
// w_MainWindow at construction.
std::vector<Coordsys_model> build_models();

// Create the initial Coordsys (axes, ranges, ticks) the views are drawn on.
Coordsys* get_initial_cs();

// Populate the QGraphicsScene with the renderable items of one model, and wire each
// active scene's key/slot connections to the view.
void populate_scene(Coordsys* cs, w_Coordsys* wcs, Coordsys_model* cm,
                    QGraphicsScene* scene);
