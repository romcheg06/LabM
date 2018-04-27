/*
 * Решение СЛАУ методом сопряженных градиентов.
 * Алгоритм взят https://www.intuit.ru/studies/courses/4447/983/lecture/14931?page=10
 * Немного оптимизирован и распараллелен
 */


#include "utils.h"

#include <memory>
#include <vector>
#include <iostream>

#include <ctime>
#include <cstdlib>

const size_t N = 4200;//!< Размер системы уравнений
const size_t ITERATIONS = 5000;

//Типы элементов
typedef double values_t;
#define MPI_VALUES_TYPE MPI_DOUBLE

int getStride(int processesCount)
{
  return N / processesCount;
}

/*!
 * \brief Проверяет соответствие размерности матрицы количеству процессов.
 * Сторона матрицы должна нацело делиться на кол-во процессов.
 * \param processesCount Количество процессов
 * \return
 */
bool assertMatrixSize(const int processesCount)
{
    const int stride = getStride(processesCount);

    return (stride * processesCount) == N;
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
    LabWorkerProcess(const int rank, const int size): WorkerProcess(rank, size),
      m_stride(getStride(size)),
      m_aMatrixPart(N * m_stride, 0.),
      m_xVectorPart(m_stride, 0.),
      m_dirVector(N, 0.),
      m_dirVectorPart(m_stride, 0.),
      m_residualVectorPart(m_stride, 0.),
      aMatrixMulDirVectorPart(m_stride, 0.)
    {
    }

    /*!
     * \brief Основной метод запускаемый в процессе
     */
    virtual void execute()
    {
        double calculationTime = MPI_Wtime();

        MPI_Scatter(NULL, m_stride * N, MPI_VALUES_TYPE,
                    m_aMatrixPart.data(), m_stride * N, MPI_VALUES_TYPE,
                    0, MPI_COMM_WORLD);

        MPI_Scatter(NULL, m_stride, MPI_VALUES_TYPE,
                    m_residualVectorPart.data(), m_stride, MPI_VALUES_TYPE,
                    0, MPI_COMM_WORLD);

        MPI_Bcast(m_dirVector.data(), m_dirVector.size(), MPI_VALUES_TYPE, 0, MPI_COMM_WORLD);

        iterate();

        MPI_Gather(m_xVectorPart.data(), m_stride, MPI_VALUES_TYPE, NULL, m_stride,
                   MPI_VALUES_TYPE, 0, MPI_COMM_WORLD);

        calculationTime = MPI_Wtime() - calculationTime;
        printf("Process %d/%d execution time: %.6f\n", m_rank, m_processesCount, calculationTime);
    }

protected:

    /*!
     * \brief Решение СЛАУ, метод общий для всех процессов
     */
    void iterate()
    {
      const values_t residualDotOldPart = matrixVectorMul(m_residualVectorPart, m_residualVectorPart)[0];
      values_t residualDotOld = 0;
      MPI_Allreduce(&residualDotOldPart, &residualDotOld, 1, MPI_VALUES_TYPE, MPI_SUM, MPI_COMM_WORLD);

      size_t iterations = ITERATIONS;
      while(iterations)
      {
        aMatrixMulDirVectorPart = matrixVectorMul(m_aMatrixPart, m_dirVector);

        const values_t alphaDivisorPart = matrixVectorMul(aMatrixMulDirVectorPart, m_dirVector.begin() + m_rank * m_stride, m_stride)[0];
        values_t alphaDivisor = 0;
        MPI_Allreduce(&alphaDivisorPart, &alphaDivisor, 1, MPI_VALUES_TYPE, MPI_SUM, MPI_COMM_WORLD);

        const values_t alpha = residualDotOld / alphaDivisor;

        for(size_t index = 0; index < static_cast<size_t>(m_stride); ++index)
        {
          m_xVectorPart[index] = m_xVectorPart[index] + alpha * m_dirVector[globalVectorIndex(index)];
          m_residualVectorPart[index] = m_residualVectorPart[index] - alpha * aMatrixMulDirVectorPart[index];
        }

        const values_t residualDotNewPart = matrixVectorMul(m_residualVectorPart, m_residualVectorPart)[0];
        values_t residualDotNew = 0;
        MPI_Allreduce(&residualDotNewPart, &residualDotNew, 1, MPI_VALUES_TYPE, MPI_SUM, MPI_COMM_WORLD);

        const values_t betta = residualDotNew / residualDotOld;
        residualDotOld = residualDotNew;

        for(size_t index = 0; index < static_cast<size_t>(m_stride); ++index)
        {
          m_dirVectorPart[index] = m_residualVectorPart[index] + betta * m_dirVector[globalVectorIndex(index)];
        }

        MPI_Allgather(m_dirVectorPart.data(), m_stride, MPI_VALUES_TYPE, m_dirVector.data(), m_stride,
                      MPI_VALUES_TYPE, MPI_COMM_WORLD);
        --iterations;
      }
    }

    /*!
     * \brief Вывод матрицы (или вектора) в консоль
     * \param matrix
     */
    void printMatrix(const std::vector<values_t>& matrix) const
    {
      size_t columnIndex = 0;
      for(auto value : matrix)
      {
        if(columnIndex == N)
        {
          std::cout << std::endl;
          columnIndex = 0;
        }
        std::cout << value << "; ";
        ++columnIndex;
      }
      std::cout << std::endl;
    }

    /*!
     * \brief Умножение матрицы на вектор (можно и 2 вектора перемножать)
     * \param matrix
     * \param vector
     * \return Результирующий вектор
     */
    std::vector<values_t> matrixVectorMul(const std::vector<values_t>& matrix, const std::vector<values_t>& vector)
    {
      std::vector<values_t> result;
      for(size_t row = 0; row < matrix.size() / vector.size(); ++row)
      {
        double value = 0;
        for(size_t column = 0; column < vector.size(); ++column)
          value += matrix[row * vector.size() + column] * vector[column];
        result.push_back(value);
      }
      return result;
    }

    /*!
     * \brief Умножение матрицы на вектор заданный итератором на начало и размером (можно и 2 вектора перемножать)
     * \param matrix
     * \param vectorIterator Итератор начала вектора
     * \param vectorSize Размер вектора
     * \return Результирующий вектор
     */
    template <typename I>
    std::vector<values_t> matrixVectorMul(const std::vector<values_t>& matrix, const I& vectorIterator, const size_t vectorSize)
    {
      std::vector<values_t> result;
      for(size_t row = 0; row < matrix.size() / vectorSize; ++row)
      {
        double value = 0;
        I iterator = vectorIterator;
        for(size_t column = 0; column < vectorSize; ++column)
          value += matrix[row * vectorSize + column] * (*(iterator++));
        result.push_back(value);
      }
      return result;
    }

    /*!
     * \brief Получение индекса элемента в общем векторе по локальному индексу части вектора в процессе
     * \param localVectorIndex Локальный индекс
     * \return Индекс в целом векторе
     */
    const size_t globalVectorIndex(const size_t localVectorIndex) const
    {
      return m_rank * m_stride + localVectorIndex;
    }

    const int m_stride;//!< Количество строк матрицы на 1 процесс
    std::vector<values_t> m_aMatrixPart;//!< Часть матрицы А
    std::vector<values_t> m_xVectorPart;//!< Часть вектора результатов
    std::vector<values_t> m_dirVector;//!< Вектор направления
    std::vector<values_t> m_dirVectorPart;//!< Часть вектора направления
    std::vector<values_t> m_residualVectorPart;//!< Часть вектора невязки
    std::vector<values_t> aMatrixMulDirVectorPart;//!< Результат aMatrixPart x dirVector

}; // end of LabWorkerProcess

/*!
 * \brief Класс главного процесса.
 */
class LabMainProcess: public LabWorkerProcess
{
public:
    /*!
     * \brief Конструктор
     * \param size Количество исполняемых процессов
     */
    LabMainProcess(const int size):
        LabWorkerProcess(0u, size),
        m_aMatrix(N * N, 0.),
        m_bVector(N, 0.),
        m_residualVector(N, 0.),
        m_xVector(N, 0.)
    {
      randomizeAMatrix();
      randomizeBVector();
    }

    virtual void execute()
    {
        double mainTime = MPI_Wtime();

        MPI_Scatter(m_aMatrix.data(), m_stride * N, MPI_VALUES_TYPE,
                    m_aMatrixPart.data(), m_stride * N, MPI_VALUES_TYPE,
                    0, MPI_COMM_WORLD);

        //m_xVector = 0
        m_residualVector = m_bVector;
        m_dirVector = m_residualVector;

        MPI_Scatter(m_residualVector.data(), m_stride, MPI_VALUES_TYPE,
                    m_residualVectorPart.data(), m_stride, MPI_VALUES_TYPE,
                    0, MPI_COMM_WORLD);

        MPI_Bcast(m_dirVector.data(), m_dirVector.size(), MPI_VALUES_TYPE, 0, MPI_COMM_WORLD);

        iterate();

        MPI_Gather(m_xVectorPart.data(), m_stride, MPI_VALUES_TYPE, m_xVector.data(), m_stride,
                   MPI_VALUES_TYPE, 0, MPI_COMM_WORLD);

        mainTime = MPI_Wtime() - mainTime;
        printf("Main process %d/%d execution time: %.6f\n", m_rank, m_processesCount, mainTime);

        double check = 0;
        const std::vector<values_t>& multVector = matrixVectorMul(m_aMatrix, m_xVector);
        for(size_t index = 0; index < N; ++index)
          check += m_bVector[index] - multVector[index];

        std::cout << "RESULT:" << std::endl;
        printMatrix(m_xVector);
        std::cout << "Check delta sum: " << check << std::endl;
    }

private:
    /*!
     * \brief Заполняет матрицу А случайными значениями (так чтобы матрица была симметрична)
     */
    void randomizeAMatrix()
    {
      srand(time(0));

      for(int row = 0; row < N; ++row)
        for(int column = N - 1; column >= 0; --column)
        {
          if(column >= row)
            aValue(column, row) = rand() % 10;
          else
            aValue(column, row) = aValue(row, column);
        }
    }

    /*!
     * \brief Заполняет вектор B случайными значениями
     */
    void randomizeBVector()
    {
      srand(time(0));
      for(auto& value : m_bVector)
        value = rand() % 10;
    }

    /*!
     * \brief Доступ к элементам матрицы A
     * \param column
     * \param row
     * \return
     */
    values_t& aValue(size_t column, size_t row)
    {
      return m_aMatrix[row * N + column];
    }

    /*!
     * \brief Константный доступ к элементам матрицы A
     * \param column
     * \param row
     * \return
     */
    const values_t aValue(size_t column, size_t row) const
    {
      return m_aMatrix[row * N + column];
    }

    std::vector<values_t> m_aMatrix;//!< Матрица A
    std::vector<values_t> m_bVector;//!< Вектор B
    std::vector<values_t> m_residualVector;//!< Вектор невязки
    std::vector<values_t> m_xVector;//!< Вектор результата

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

    if(!assertMatrixSize(size))
    {
       if(rank == 0)
            std::cout << "Program can not be run with current matrix size and " << size << " processes!" << std::endl;
       MPI_Finalize(); /* ends MPI */
       return 0;
    }

    std::unique_ptr<Process> process = makeProcess(rank, size);
    process->execute();

    MPI_Finalize(); /* ends MPI */
    return 0;
}
