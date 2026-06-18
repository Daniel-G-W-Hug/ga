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
        // a "missing font family" warning at startup. Helvetica matches the family the
        // axis labels / legends already use.
        QApplication::setFont(QFont("Helvetica", 10));

        w_MainWindow w;
        w.move(400, 100);
        w.show();
        return a.exec();
    }
    catch (std::exception const& e) {
        std::cout << e.what();
    }
}
