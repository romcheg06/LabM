#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class QFileSystemWatcher;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void updateTable();

private:
    Ui::MainWindow *ui;
    QFileSystemWatcher* m_fileSystemWatcher;
};

#endif // MAINWINDOW_H
