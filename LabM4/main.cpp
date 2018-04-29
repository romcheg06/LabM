#include "lab4types.h"
#include "extendedslice.h"
#include "utils.h"

#include <memory>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

//Вывод результирующей матрицы в консоль
//#define PRINT_MATRIX

const values_t X_SIDE = 1.;
const values_t Y_SIDE = 2.;

//Граничные условия
const values_t UPPER_BOUNDARY_VALUE = 100.;
const values_t LOWER_BOUNDARY_VALUE = 0.;
const values_t LEFT_BOUNDARY_VALUE = 0.;
const values_t RIGHT_BOUNDARY_VALUE = 0.;

const size_t X_POINTS_COUNT = 240;//!< Размер поля по ширине
const size_t Y_POINTS_COUNT = 480;//!< Размер поля по высоте

const values_t H = X_SIDE / X_POINTS_COUNT;//!< Натуральный шаг сетки
const values_t EPSILON = .0001;//!< Точность сходимости невязки

const size_t STEPS_PER_ITERATION = 10u;//!<  Количество шагов в одной итерации между которыми не вычисляется невязка
const size_t ITERATIONS = 5000u;//!< Количество итераций (одна итерация = STEPS_PER_ITERATION шагов)

//Тэги сообщений
const int SENT_UP_BOUND_TAG = 1;
const int SENT_DOWN_BOUND_TAG = 2;
const int SENT_LEFT_BOUND_TAG = 3;
const int SENT_RIGHT_BOUND_TAG = 4;
const int SENT_SLICE_TAG = 5;

/*!
 * \brief Размерности разбиения процессов в сетке
 */
typedef std::pair<int/*кол-во по горизонтали*/, int/*кол-во по вертикали*/> ProcessesDims;

/*!
 * \brief олучить размерности разбиения процессов в сетке.
 * \param processesCount Количество процессов
 * \return размерности процессов.
 */
ProcessesDims getProcessesDims(const int processesCount)
{
    int dimX = static_cast<int>(sqrt(processesCount));
    int dimY = processesCount / dimX;

    while(processesCount != dimX * dimY )
    {
        ++dimY;
        dimX = processesCount / dimY;
    }

    return std::make_pair(dimX, dimY);
}

/*!
 * \brief Проверяет соответствие размерности множества точек количеству процессов.
 * множество точек должно делиться на целое количество одинаковых кусков в соответствии с разбиением процессов
 * \param processesCount Количество процессов
 * \return
 */
bool assertFieldSize(const int processesCount)
{
    const ProcessesDims processesDims = getProcessesDims(processesCount);

    const size_t xDiv = X_POINTS_COUNT / processesDims.first;
    if(X_POINTS_COUNT != xDiv * processesDims.first)
        return false;

    const size_t yDiv = Y_POINTS_COUNT / processesDims.second;
    if(Y_POINTS_COUNT != yDiv * processesDims.second)
        return false;

    return true;
}

/*!
 * \brief Размерности разбиения поля на куски.
 */
typedef std::pair<size_t/*горизонталь*/, size_t/*вертикаль*/> FieldStrides;

/*!
 * \brief Получить размерности разбиения поля на куски по процессам.
 * \param processesCount Количество процессов
 * \return размерности одного куска.
 */
FieldStrides getFieldStrides(const int processesCount)
{
    const ProcessesDims processesDims = getProcessesDims(processesCount);

    return std::make_pair(X_POINTS_COUNT / static_cast<size_t>(processesDims.first),
                          Y_POINTS_COUNT / static_cast<size_t>(processesDims.second));
}

/*!
 * \brief Создает коммуникатор декартовой топологии
 * \param comm
 * \param processesCount
 */
void createNetCommunicator(MPI_Comm* const comm, const int processesCount)
{
    const ProcessesDims processesDims = getProcessesDims(processesCount);

    int dims[2] = {processesDims.first, processesDims.second};
    int periods[2] {0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, comm);
}


typedef std::pair<values_t/*невязка*/, size_t/*количество пройденных итераций*/> IterationsResult;

/*!
 * \brief Производит итерации методом Зейделя с красно/черным разбиением точек
 * общая функция для всех процессов
 * \param extendedSlice кусок поля
 * \param netComm коммуникатор декартовой топологии
 */
IterationsResult doZeidelIterations(ExtendedSlice& extendedSlice, const MPI_Comm netComm)
{
    int upperRank, lowerRank, leftRank, rightRank;

    const int bufSize = (X_POINTS_COUNT + Y_POINTS_COUNT) * sizeof(values_t) + MPI_BSEND_OVERHEAD;
    static values_t buf[bufSize];

    size_t globalIndex = 0;

    MPI_Cart_shift(netComm, 0, 1, &leftRank, &rightRank);
    MPI_Cart_shift(netComm, 1, 1, &upperRank, &lowerRank);

    int rank = 0;
    int coordinates[2] = {0, 0};
    MPI_Comm_rank(netComm, &rank);
    MPI_Cart_coords(netComm, rank, 2, coordinates);

    globalIndex = coordinates[0] * extendedSlice.m_strideX +
                  coordinates[1] * extendedSlice.m_strideY * X_POINTS_COUNT;

    MPI_Buffer_attach(buf, bufSize);

    if(upperRank == MPI_PROC_NULL)
      std::fill(extendedSlice.m_upperBound.begin(), extendedSlice.m_upperBound.end(), UPPER_BOUNDARY_VALUE);
    if(lowerRank == MPI_PROC_NULL)
      std::fill(extendedSlice.m_lowerBound.begin(), extendedSlice.m_lowerBound.end(), LOWER_BOUNDARY_VALUE);
    if(leftRank == MPI_PROC_NULL)
      std::fill(extendedSlice.m_leftExtendedBound.begin(), extendedSlice.m_leftExtendedBound.end(), LEFT_BOUNDARY_VALUE);
    if(rightRank == MPI_PROC_NULL)
      std::fill(extendedSlice.m_rightExtendedBound.begin(), extendedSlice.m_rightExtendedBound.end(), RIGHT_BOUNDARY_VALUE);

    values_t oldResidual = 1000;

    size_t iterations = ITERATIONS;
    while(iterations)
    {
        size_t steps = STEPS_PER_ITERATION;
        while(steps)
        {
          for(size_t colorStep = 0; colorStep < 2; ++colorStep)
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
                  std::vector<values_t> firstColumn = extendedSlice.getFirstColumn();
                  MPI_Bsend(firstColumn.data(), firstColumn.size(),
                            MPI_VALUES_TYPE, leftRank, SENT_LEFT_BOUND_TAG, netComm);
              }

              if(rightRank != MPI_PROC_NULL)
              {
                  MPI_Recv(&(extendedSlice.m_rightExtendedBound.data()[1]), extendedSlice.m_rightExtendedBound.size() - 2,
                           MPI_VALUES_TYPE, rightRank, SENT_LEFT_BOUND_TAG, netComm, &status);

                  std::vector<values_t> lastColumn = extendedSlice.getLastColumn();
                  MPI_Bsend(lastColumn.data(), lastColumn.size(),
                            MPI_VALUES_TYPE, rightRank, SENT_RIGHT_BOUND_TAG, netComm);
              }

              if(leftRank != MPI_PROC_NULL)
              {
                  MPI_Recv(&(extendedSlice.m_leftExtendedBound.data()[1]), extendedSlice.m_leftExtendedBound.size() - 2,
                           MPI_VALUES_TYPE, leftRank, SENT_RIGHT_BOUND_TAG, netComm, &status);
              }

              if(colorStep == 0 && steps == STEPS_PER_ITERATION)
              {
                const values_t residualPart = extendedSlice.maxResidual(H);
                values_t residual = 0;
                MPI_Allgather(&residualPart, 1, MPI_VALUES_TYPE, &residual, 1,
                              MPI_VALUES_TYPE, netComm);

                if(std::abs(oldResidual - residual) < EPSILON)
                  return std::make_pair(residual, ITERATIONS - iterations);

                oldResidual = residual;
              }

              const ExtendedSlice::ZeidelStepColor zeidelColor =
                  static_cast<ExtendedSlice::ZeidelStepColor>((globalIndex + colorStep) % ExtendedSlice::INVALID_COLOR);
              extendedSlice.zeidelStep(zeidelColor);

              MPI_Barrier(netComm);
            }
            --steps;
        }
        --iterations;
    }

    return std::make_pair(oldResidual, ITERATIONS - iterations);;
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
    LabWorkerProcess(const int rank, const int size): WorkerProcess(rank, size)
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

        doZeidelIterations(extendedSlice, netComm);

        MPI_Send(slice.m_values.data(), slice.m_values.size(),
                 MPI_VALUES_TYPE, mainProcessNetRank, SENT_SLICE_TAG, netComm);


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
    LabMainProcess(const int size):
        MainProcess(size),
        m_field(X_POINTS_COUNT, Y_POINTS_COUNT)
    {
        const FieldStrides fieldStrides = getFieldStrides(m_processesCount);
        const size_t strideX = fieldStrides.first;
        const size_t strideY = fieldStrides.second;

        for(size_t globalX = 0; globalX < X_POINTS_COUNT; globalX += strideX)
            for(size_t globalY = 0; globalY < Y_POINTS_COUNT; globalY += strideY)
            {
                m_field.m_slices.push_back(
                    Slice::makeZeroSlice(globalX, globalY, strideX, strideY));
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

        //Рассылаем куски
        for(int process = 0; process < m_processesCount; ++process)
        {
            if(process == netRank)
                continue;
            MPI_Send(m_field.m_slices[process].m_values.data(), m_field.m_slices[process].m_values.size(),
                     MPI_VALUES_TYPE, process, 0, netComm);
        }

        MPI_Barrier(netComm);

        ExtendedSlice extendedSlice(m_field.m_slices[mySliceNumber]);

        const IterationsResult iterationsResult = doZeidelIterations(extendedSlice, netComm);

        //Собираем куски
        for(int process = 0; process < m_processesCount - 1; ++process)
        {

            MPI_Probe(MPI_ANY_SOURCE, SENT_SLICE_TAG, netComm, &status);
            const int sender = status.MPI_SOURCE;
            MPI_Recv(m_field.m_slices[sender].m_values.data(), m_field.m_slices[sender].m_values.size(),
                     MPI_VALUES_TYPE, sender, SENT_SLICE_TAG, netComm, &status);
        }

        mainTime = MPI_Wtime() - mainTime;
        printf("Main process %d/%d execution time: %.6f\n", m_rank, m_processesCount, mainTime);

#ifdef PRINT_MATRIX
        std::vector<values_t> wholeField(X_POINTS_COUNT * Y_POINTS_COUNT, 0.);

        for(const Slice& slice : m_field.m_slices)
          for(size_t sliceIndex = 0; sliceIndex < slice.m_values.size(); ++sliceIndex)
          {
            const size_t globalX = slice.m_globalX + sliceIndex % slice.m_stride;
            const size_t globalY = slice.m_globalY + sliceIndex / slice.m_stride;
            const size_t globalIndex = globalX + globalY * X_POINTS_COUNT;
            wholeField[globalIndex] = slice.m_values[sliceIndex];
          }

        std::cout << "Result: " << std::endl;
        for(size_t yIndex = 0; yIndex < Y_POINTS_COUNT; ++yIndex)
        {
          std::ostringstream rowStream;
          for(size_t xIndex = 0; xIndex < X_POINTS_COUNT; ++ xIndex)
          {
            rowStream << std::setw(7) << std::fixed << std::showpoint << std::setprecision(4) << wholeField[xIndex + yIndex * Y_POINTS_COUNT]<< "; ";
          }
          rowStream << std::endl;

          std::cout << rowStream.str();
        }
#endif
        std::ostringstream itResultStream;
        itResultStream << " Residual: " << iterationsResult.first << " Iterations: " << iterationsResult.second << std::endl;
        std::cout << itResultStream.str();

        //Возможная запись в файл
//        {
//            std::ofstream os("field.dat", std::ios::binary);
//            boost::archive::binary_oarchive oar(os);
//            oar << m_field;
//        }
    }

private:
    Field m_field;

}; // end of MainProcess

std::unique_ptr<Process> makeProcess(const int rank, const int size)
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

    if(!assertFieldSize(size))
    {
       if(rank == 0)
            std::cout << "Program can not be run with current field size and " << size << " processes!" << std::endl;
       MPI_Finalize(); /* ends MPI */
       return 0;
    }

    std::unique_ptr<Process> process = makeProcess(rank, size);
    process->execute();

    MPI_Finalize(); /* ends MPI */
    return 0;
}

