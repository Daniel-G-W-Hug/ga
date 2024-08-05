//
// author: Daniel Hug, 2024
//

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
    catch (const std::exception& e) {
        std::cout << e.what();
    }
}
