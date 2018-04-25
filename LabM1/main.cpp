#include "utils.h"

#include <assert.h>
#include <memory>
#include <vector>
#include <iostream>
#include <cmath>

static const double REAL_RADIUS = 2.;
static const int a = 50000;
static const int b = 20000;
static const size_t MAX_ITERATIONS = 100;

typedef std::pair<double, double> ComplexNumber;

ComplexNumber getComplex(const int i, const int j)
{
    static const double X_STEP = REAL_RADIUS * 2. / static_cast<double>((a - 1));
    static const double Y_STEP = REAL_RADIUS * 2. / static_cast<double>((b - 1));

    const double re = (0. - REAL_RADIUS) + static_cast<double>(i) * X_STEP;
    const double im = (0. - REAL_RADIUS) + static_cast<double>(j) * Y_STEP;

    return std::make_pair(re, im);
}

int getI(const int index)
{
    assert(index >= 0);
    assert(index <= a * b);
    return index % a;
}

int getJ(const int index)
{
    assert(index >= 0);
    assert(index <= a * b);
    return index / a;
}

double getComplexMod(const ComplexNumber& complex)
{
    return sqrt(complex.first * complex.first + complex.second * complex.second);
}

ComplexNumber operator+(const ComplexNumber& lhs, const ComplexNumber& rhs)
{
    return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

ComplexNumber complexPow2(const ComplexNumber& complex)
{
    return std::make_pair(complex.first * complex.first - complex.second * complex.second,
                          complex.second * complex.first + complex.first * complex.second);
}

bool checkComlex(const ComplexNumber& complex)
{
    ComplexNumber currentComplex = std::make_pair(0., 0.);
    bool result = true;
    for(size_t iteration = 0; iteration < MAX_ITERATIONS; ++iteration)
    {
        currentComplex = complexPow2(currentComplex) + complex;
        if(getComplexMod(currentComplex) > 2.)
        {
            result = false;
            break;
        }
    }
    return result;
}

std::vector<int> getResults(const int startingIndex, const int maxIndex, const int step)
{
    std::vector<int> results;

    for(int index = startingIndex; index < maxIndex; index += step)
    {
        if(checkComlex(getComplex(getI(index), getJ(index))))
            results.push_back(index);
    }

    return results;
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
        const std::vector<int>& results = getResults(m_rank, a * b, m_processesCount);
        calculationTime = MPI_Wtime() - calculationTime;
        printf("Process %d/%d execution time: %.6f\n", m_rank, m_processesCount, calculationTime);

        MPI_Send(results.data(), results.size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
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
    LabMainProcess(const int size): MainProcess(size)
    {
    }

    virtual void execute()
    {
        double mainTime = MPI_Wtime();

        std::vector<int> results = getResults(m_rank, a * b, m_processesCount);

        const double calculationTime = MPI_Wtime() - mainTime;
        printf("Process %d/%d execution time: %.6f\n", m_rank, m_processesCount, calculationTime);

        for(unsigned process = 1; process < m_processesCount; ++process)
        {
            const size_t resultsOldSize = results.size();
            MPI_Status status;
            int messageSize = 0;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &messageSize);

            results.resize(resultsOldSize + messageSize);
            MPI_Recv(&(results.data()[resultsOldSize]), messageSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }

        mainTime = MPI_Wtime() - mainTime;

        std::cout << "Found " << results.size() << " points of total " << (a * b) <<" :" << std::endl;
        for(const auto& index : results)
        {
            const int i = getI(index);
            const int j = getJ(index);
            const ComplexNumber complexNumber = getComplex(i, j);
            //std::cout << "a: " << i << " b: " << j << " re: " << complexNumber.first << " im: " << complexNumber.second << std::endl;
        }

        std::cout << "Results end" << std::endl;
        std::cout << "Main process execution time: " << mainTime << std::endl;
        std::cout << "Processes: " << m_processesCount << std::endl;
    }
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

    std::unique_ptr<Process> process = makeProcess(rank, size);
    process->execute();

    MPI_Finalize(); /* ends MPI */
    return 0;
}
