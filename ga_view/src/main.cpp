// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "w_mainwindow.hpp"

#include <QApplication>

#include <iostream>


int main(int argc, char* argv[])
{

    try {
        QApplication a(argc, argv);

        w_MainWindow w;
        w.move(400, 100);
        w.show();
        return a.exec();
    }
    catch (std::exception const& e) {
        std::cout << e.what();
    }
}
