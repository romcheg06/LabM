#ifndef EXTENDEDSLICE_H
#define EXTENDEDSLICE_H

#include "slice.h"

/*!
 * \brief расширенный кусок поля с границами и возможностью сделать один шаг игры
 */
struct ExtendedSlice
{
    explicit ExtendedSlice(Slice& slice):
                m_slice(slice),
                m_strideX{m_slice.m_stride},
                m_strideY{m_slice.m_values.size() / m_strideX}
    {
        m_upperBound.resize(m_strideX, 0);
        m_lowerBound.resize(m_strideX, 0);
        m_leftExtendedBound.resize(m_strideY + 2, 0);
        m_rightExtendedBound.resize(m_strideY + 2, 0);
    }

    values_t value(const size_t x, const size_t y) const
    {
        if(x == 0)
            return m_leftExtendedBound[y];

        if(x == m_strideX + 1)
            return m_rightExtendedBound[y];

        if(y == 0)
            return m_upperBound[x - 1];

        if(y == m_strideY + 1)
            return m_lowerBound[x - 1];

        return m_slice.m_values[(y - 1) * m_strideX + x - 1];
    }

    values_t& value(const size_t x, const size_t y)
    {
        if(x == 0)
            return m_leftExtendedBound[y];

        if(x == m_strideX + 1)
            return m_rightExtendedBound[y];

        if(y == 0)
            return m_upperBound[x - 1];

        if(y == m_strideY + 1)
            return m_lowerBound[x - 1];

        return m_slice.m_values[(y - 1) * m_strideX + x - 1];
    }

    std::vector<values_t> getUpperRow() const
    {
        return std::vector<values_t>{m_slice.m_values.begin(), m_slice.m_values.begin() + m_strideX};
    }

    std::vector<values_t> getLowerRow() const
    {
        return std::vector<values_t>{m_slice.m_values.end() - m_strideX, m_slice.m_values.end()};
    }

    std::vector<values_t> getFirstExtendedColumn() const
    {
        std::vector<values_t> column;
        column.push_back(m_upperBound[0]);
        for(size_t i = 0; i < m_strideY; ++i)
            column.push_back(m_slice.m_values[i * m_strideX]);
        column.push_back(m_lowerBound[0]);

        return column;
    }

    std::vector<values_t> getLastExtendedColumn() const
    {
        std::vector<values_t> column;
        column.push_back(m_upperBound[m_strideX - 1]);
        for(size_t i = 0; i < m_strideY; ++i)
            column.push_back(m_slice.m_values[i * m_strideX + m_strideX - 1]);
        column.push_back(m_lowerBound[m_strideX - 1]);

        return column;
    }

    void lifeStep()
    {
        std::vector<values_t> nextStepSliceValues(m_slice.m_values);

        for(size_t x = 1; x <= m_strideX; ++x)
            for(size_t y = 1; y <= m_strideY; ++y)
            {
                size_t aliveNeighbors = 0;
                for(size_t nx = x - 1; nx < x + 2; ++nx)
                    for(size_t ny = y - 1; ny < y + 2; ++ny)
                    {
                        if(nx == x && ny == y)
                            continue;
                        if(value(nx, ny))
                           ++aliveNeighbors;
                    }
                if(value(x, y))
                {
                    if(aliveNeighbors < 2 || aliveNeighbors > 3)
                        nextStepSliceValues[m_strideX * (y - 1) + x - 1] = 0;
                }
                else
                {
                    if(aliveNeighbors == 3)
                        nextStepSliceValues[m_strideX * (y - 1) + x - 1] = 1;
                }
            }

        m_slice.m_values = nextStepSliceValues;
    }

    Slice& m_slice;
    const size_t m_strideX;
    const size_t m_strideY;
    std::vector<values_t> m_upperBound;
    std::vector<values_t> m_lowerBound;
    std::vector<values_t> m_leftExtendedBound;
    std::vector<values_t> m_rightExtendedBound;
};

#endif // EXTENDEDSLICE_H
