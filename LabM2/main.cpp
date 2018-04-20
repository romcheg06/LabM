#include "utils.h"
#include "extendedslice.h"

#include <memory>
#include <vector>
#include <fstream>

//Если определено, приложение вставляет задержки после каждой итерации
#define DELAYS

//Если определено, производится периодическая запись поля в файл
#define FILE_SAVE

//Если определено, поле заполняется тестовым примером иначе - случайно
#define EXAMPLE

#ifdef DELAYS
#   include <chrono>
#   include <thread>
    const int64_t DELAY = 400;//!< Задержка между итерациями в мс
#endif

const int PERIODIC_FIELD = 1;//!< Периодическое поле (1 - да / 0 - нет)
const size_t FIELD_X_SIZE = 20u;//!< Размер поля по ширине
const size_t FIELD_Y_SIZE = 20u;//!< Размер поля по высоте
const size_t STEPS_PER_ITERATION = 1u;//!<  Количество шагов игры между сбросами состояния поля на диск
const size_t ITERATIONS = 500u;//!< Количество итераций (одна итерация = STEPS_PER_ITERATION шагов)

const int SENT_UP_BOUND_TAG = 1;
const int SENT_DOWN_BOUND_TAG = 2;
const int SENT_LEFT_BOUND_TAG = 3;
const int SENT_RIGHT_BOUND_TAG = 4;

/*!
 * \brief Размерности разбиения процессов в сетке
 */
typedef std::pair<int/*кол-во по горизонтали*/, int/*кол-во по вертикали*/> ProcessesDims;

/*!
 * \brief олучить размерности разбиения процессов в сетке. Приложение рассчитано на 1, 2 или 4 процесса
 * \param processesCount Количество процессов
 * \return размерности процессов.
 */
ProcessesDims getProcessesDims(int processesCount)
{
    int dimX = processesCount > 1 ? 2 : 1;
    int dimY = processesCount == 4 ? 2 : 1;
    return std::make_pair(dimX, dimY);
}

/*!
 * \brief Размерности разбиения поля на куски.
 */
typedef std::pair<size_t/*горизонталь*/, size_t/*вертикаль*/> FieldStrides;

/*!
 * \brief Получить размерности разбиения поля на куски по процессам. Приложение рассчитано на 1, 2 или 4 процесса
 * \param processesCount Количество процессов
 * \return размерности одного куска.
 */
FieldStrides getFieldStrides(int processesCount)
{
    const ProcessesDims processesDims = getProcessesDims(processesCount);

    return std::make_pair(FIELD_X_SIZE / static_cast<size_t>(processesDims.first),
                          FIELD_Y_SIZE / static_cast<size_t>(processesDims.second));
}

/*!
 * \brief Создает коммуникатор декартовой топологии
 * \param comm
 * \param processesCount
 */
void createNetCommunicator(MPI_Comm* comm, int processesCount)
{
    const ProcessesDims processesDims = getProcessesDims(processesCount);

    int dims[2] = {processesDims.first, processesDims.second};
    int periods[2] {PERIODIC_FIELD, PERIODIC_FIELD};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, comm);
}

/*!
 * \brief Проводит несколько шагов игры на куске поля
 * \param extendedSlice кусок поля
 * \param netComm коммуникатор декартовой топологии
 */
void doLifeSteps(ExtendedSlice& extendedSlice, MPI_Comm netComm)
{
    static bool firstRun = true;
    static int upperRank, lowerRank, leftRank, rightRank;

    const int bufSize = FIELD_X_SIZE + FIELD_Y_SIZE + MPI_BSEND_OVERHEAD;
    static values_t buf[bufSize];

    if(firstRun)
    {
        MPI_Cart_shift(netComm, 0, 1, &leftRank, &rightRank);
        MPI_Cart_shift(netComm, 1, 1, &upperRank, &lowerRank);
        MPI_Buffer_attach(buf, bufSize);
        firstRun = false;
    }

    size_t steps = STEPS_PER_ITERATION;
    while(steps)
    {
        MPI_Status status;

        if(upperRank != MPI_PROC_NULL)
        {
            std::vector<values_t> firstRow = extendedSlice.getUpperRow();
            MPI_Bsend(firstRow.data(), firstRow.size(),
                      MPI_VALUES_TYPE, upperRank, SENT_UP_BOUND_TAG, netComm);
        }

        if(lowerRank != MPI_PROC_NULL)
        {
            MPI_Recv(extendedSlice.m_lowerBound.data(), extendedSlice.m_lowerBound.size(),
                     MPI_VALUES_TYPE, lowerRank, SENT_UP_BOUND_TAG, netComm, &status);

            std::vector<values_t> lastRow = extendedSlice.getLowerRow();
            MPI_Bsend(lastRow.data(), lastRow.size(),
                      MPI_VALUES_TYPE, lowerRank, SENT_DOWN_BOUND_TAG, netComm);
        }
        if(upperRank != MPI_PROC_NULL)
        {
            MPI_Recv(extendedSlice.m_upperBound.data(), extendedSlice.m_upperBound.size(),
                     MPI_VALUES_TYPE, upperRank, SENT_DOWN_BOUND_TAG, netComm, &status);
        }

        if(leftRank != MPI_PROC_NULL)
        {
            std::vector<values_t> firstColumn = extendedSlice.getFirstExtendedColumn();
            MPI_Bsend(firstColumn.data(), firstColumn.size(),
                      MPI_VALUES_TYPE, leftRank, SENT_LEFT_BOUND_TAG, netComm);
        }

        if(rightRank != MPI_PROC_NULL)
        {
            MPI_Recv(extendedSlice.m_rightExtendedBound.data(), extendedSlice.m_rightExtendedBound.size(),
                     MPI_VALUES_TYPE, rightRank, SENT_LEFT_BOUND_TAG, netComm, &status);

            std::vector<values_t> lastColumn = extendedSlice.getLastExtendedColumn();
            MPI_Bsend(lastColumn.data(), lastColumn.size(),
                      MPI_VALUES_TYPE, rightRank, SENT_RIGHT_BOUND_TAG, netComm);
        }

        if(leftRank != MPI_PROC_NULL)
        {
            MPI_Recv(extendedSlice.m_leftExtendedBound.data(), extendedSlice.m_leftExtendedBound.size(),
                     MPI_VALUES_TYPE, leftRank, SENT_RIGHT_BOUND_TAG, netComm, &status);
        }

        extendedSlice.lifeStep();

        --steps;
    }
}

/*!
 * \brief Рабочий процесс (rank > 0).
 */
class LabWorkerProcess: public WorkerProcess
{
public:

    /*!
     * \brief Конструктор
     * \param rank номер процесса
     * \param size общее число запущенных процессов
     */
    LabWorkerProcess(int rank, int size): WorkerProcess(rank, size)
    {}

    /*!
     * \brief Основной метод запускаемый в процессе
     */
    virtual void execute()
    {
        double calculationTime = MPI_Wtime();

        MPI_Comm netComm;
        createNetCommunicator(&netComm, m_processesCount);

        const FieldStrides fieldStrides = getFieldStrides(m_processesCount);
        const size_t strideX = fieldStrides.first;
        const size_t strideY = fieldStrides.second;

        Slice slice(0, 0, strideX);
        slice.m_values.resize(strideX * strideY, 0);

        MPI_Status status;

        MPI_Recv(slice.m_values.data(), slice.m_values.size(),
                 MPI_VALUES_TYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, netComm, &status);

        MPI_Barrier(netComm);

        const int mainProcessNetRank = status.MPI_SOURCE;

        ExtendedSlice extendedSlice(slice);

        size_t iterations = ITERATIONS;

        while(iterations)
        {
            doLifeSteps(extendedSlice, netComm);

            MPI_Send(slice.m_values.data(), slice.m_values.size(),
                     MPI_VALUES_TYPE, mainProcessNetRank, 0, netComm);

            --iterations;

#ifdef DELAYS
            MPI_Barrier(netComm);
#endif
        }

        calculationTime = MPI_Wtime() - calculationTime;
        printf("Process %d/%d execution time: %.6f\n", m_rank, m_processesCount, calculationTime);
    }
}; // end of LabWorkerProcess

/*!
 * \brief Класс главного процесса.
 */
class LabMainProcess: public MainProcess
{
public:
    /*!
     * \brief Конструктор
     * \param size Количество исполняемых процессов
     */
    LabMainProcess(int size):
        MainProcess(size),
        m_field(FIELD_X_SIZE, FIELD_Y_SIZE)
    {
        const FieldStrides fieldStrides = getFieldStrides(m_processesCount);
        const size_t strideX = fieldStrides.first;
        const size_t strideY = fieldStrides.second;

        for(size_t globalX = 0; globalX < FIELD_X_SIZE; globalX += strideX)
            for(size_t globalY = 0; globalY < FIELD_Y_SIZE; globalY += strideY)
            {
#ifndef EXAMPLE
                m_field.m_slices.push_back(
                    Slice::makeRandomSlice(globalX, globalY, strideX, strideY));
#else
                Slice slice = Slice::makeZeroSlice(globalX, globalY, strideX, strideY);
                if(globalX == 0 && globalY == 0)
                {
                    //Запустим, к примеру, планер
                    slice.m_values[1] = 1;
                    slice.m_values[strideX + 2] = 1;
                    slice.m_values[strideX * 2] = 1;
                    slice.m_values[strideX * 2 + 1] = 1;
                    slice.m_values[strideX * 2 + 2] = 1;
                }
                m_field.m_slices.push_back(slice);
#endif
            }
    }

    virtual void execute()
    {
        double mainTime = MPI_Wtime();

        MPI_Status status;

        MPI_Comm netComm;
        createNetCommunicator(&netComm, m_processesCount);

        int netRank;
        MPI_Comm_rank (MPI_COMM_WORLD, &netRank);

        const size_t mySliceNumber = netRank;

        //Рассылаем куски поля
        for(int process = 0; process < m_processesCount; ++process)
        {
            if(process == netRank)
                continue;
            MPI_Send(m_field.m_slices[process].m_values.data(), m_field.m_slices[process].m_values.size(),
                     MPI_VALUES_TYPE, process, 0, netComm);
        }

        MPI_Barrier(netComm);

        ExtendedSlice extendedSlice(m_field.m_slices[mySliceNumber]);

        size_t iterations = ITERATIONS;
        while(iterations)
        {
            --iterations;

            doLifeSteps(extendedSlice, netComm);

            //Собираем куски поля
            for(int process = 0; process < m_processesCount - 1; ++process)
            {

                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, netComm, &status);
                const int sender = status.MPI_SOURCE;
                MPI_Recv(m_field.m_slices[sender].m_values.data(), m_field.m_slices[sender].m_values.size(),
                         MPI_VALUES_TYPE, sender, MPI_ANY_TAG, netComm, &status);
            }

#ifdef FILE_SAVE
            {
                std::ofstream os("field.dat", std::ios::binary);
                boost::archive::binary_oarchive oar(os);
                oar << m_field;
            }
#endif

#ifdef DELAYS
            std::this_thread::sleep_for(std::chrono::milliseconds(DELAY));
            MPI_Barrier(netComm);
#endif
        }

        mainTime = MPI_Wtime() - mainTime;
        printf("Main process %d/%d execution time: %.6f\n", m_rank, m_processesCount, mainTime);
    }

private:
    Field m_field;

}; // end of MainProcess

std::unique_ptr<Process> makeProcess(int rank, int size)
{
    if(rank == 0)
        return std::unique_ptr<Process>(new LabMainProcess(size));
    else
        return std::unique_ptr<Process>(new LabWorkerProcess(rank, size));
}

int main(int argc, char *argv[])
{
    int rank, size;

    MPI_Init (&argc, &argv); /* starts MPI */

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */

    if(size > 4 || size == 3)
    {
        std::cout << "Only 1, 2 or 4 processes supported!" << std::endl;
        MPI_Finalize(); /* ends MPI */
        return 0;
    }

    std::unique_ptr<Process> process = makeProcess(rank, size);
    process->execute();

    MPI_Finalize(); /* ends MPI */
    return 0;
}
