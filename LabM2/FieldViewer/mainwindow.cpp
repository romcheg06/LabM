#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "lab2types.h"
#include "extendedslice.h"
#include <QDebug>
#include <QStyle>
#include <QDesktopWidget>
#include <QTableWidget>
#include <QFileSystemWatcher>

#include <vector>
#include <iostream>
#include <fstream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setGeometry(QStyle::alignedRect(Qt::LeftToRight,
                                    Qt::AlignCenter,
                                    size(),
                                    qApp->desktop()->availableGeometry()
                                    )
                );

    updateTable();

    m_fileSystemWatcher = new QFileSystemWatcher(this);
    m_fileSystemWatcher->addPath("./field.dat");
    connect(m_fileSystemWatcher, &QFileSystemWatcher::fileChanged, this, &MainWindow::updateTable);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateTable()
{
    try
    {
        std::ifstream is("field.dat", std::ios::binary);
        boost::archive::binary_iarchive iar(is);
        Field field;
        iar >> field;

        ui->tableWidget->setRowCount(0);
        ui->tableWidget->setRowCount(field.m_height);
        ui->tableWidget->setColumnCount(field.m_width);

        for(Slice& slice : field.m_slices)
        {
            int x = slice.m_globalX;
            int y = slice.m_globalY;
            for(values_t& value : slice.m_values)
            {
                if(value)
                {
                    QTableWidgetItem *newItem = new QTableWidgetItem("*");
                    newItem->setTextAlignment(Qt::AlignCenter);
                    ui->tableWidget->setItem(y, x, newItem);
                }
                ++x;
                if(x == static_cast<int>(slice.m_globalX + slice.m_stride))
                {
                    x = slice.m_globalX;
                    ++y;
                }
            }
        }

        ui->tableWidget->resizeColumnsToContents();
        ui->tableWidget->resizeRowsToContents();
    }
    catch(...)
    {
    }
}
