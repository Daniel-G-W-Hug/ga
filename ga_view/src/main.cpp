// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "w_mainwindow.hpp"

#include <QApplication>
#include <QFont>

#include <iostream>


int main(int argc, char* argv[])
{

    try {
        QApplication a(argc, argv);

        // Set an explicit default font (an installed family). Without this Qt falls back
        // to the alias "Sans Serif", which on macOS triggers a slow font-alias lookup and
        // a "missing font family" warning at startup. Arial matches the family the axis
        // labels / legends already use, and is present on both macOS and Windows with
        // near-identical metrics (Helvetica is macOS-only -> Qt substitutes it on
        // Windows, shifting glyph widths). See view_fonts.hpp.
        QApplication::setFont(QFont("Arial", 10));

        w_MainWindow w;
        w.move(400, 100);
        w.show();
        return a.exec();
    }
    catch (std::exception const& e) {
        std::cout << e.what();
    }
}
