#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

/*!
 * \brief Абстрактный базовый класс процесса
 */
class Process
{
public:
    /*!
     * \brief Конструктор
     * \param size общее число запущенных процессов
     */
    explicit Process(const int size): m_processesCount(size)
    {}

    Process(const Process&) = delete;
    Process& operator=(const Process&) = delete;

    /*!
     * \brief Запуск процесса
     */
    virtual void execute() = 0;

    /*!
     * \brief Виртуальный деструктор
     */
    virtual ~Process()
    {}
protected:
    const unsigned m_processesCount;//!< Количество исполняемых процессов
}; // end of Process

/*!
 * \brief Рабочий процесс (rank > 0).
 */
class WorkerProcess: public Process
{
public:

    /*!
     * \brief Конструктор
     * \param rank номер процесса
     * \param size общее число запущенных процессов
     */
    WorkerProcess(const int rank, const int size):
        Process(size),
        m_rank(rank)
    {}

protected:
    const int m_rank; //!< rank текущего процесса
}; // end of WorkerProcess

/*!
 * \brief Класс главного процесса.
 */
class MainProcess: public WorkerProcess
{
public:
    /*!
     * \brief Конструктор
     * \param size Количество исполняемых процессов
     */
    MainProcess(const int size): WorkerProcess(0, size)
    {
    }
}; // end of MainProcess

#endif // UTILS_H
