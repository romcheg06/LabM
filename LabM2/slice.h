#ifndef SLICE_H
#define SLICE_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>
#include <ctime>
#include <cstdlib>

typedef char values_t;
#define MPI_VALUES_TYPE MPI_CHAR

/*!
 * \brief Кусок поля
 */
struct Slice
{
    Slice():
        m_globalX{0},
        m_globalY{0},
        m_stride{0}
    {}

    Slice(const size_t globalX, const size_t globalY, const size_t stride):
        m_globalX{globalX},
        m_globalY{globalY},
        m_stride{stride}
    {}

    static Slice makeRandomSlice(const size_t globalX, const size_t globalY, const size_t strideX, const size_t strideY)
    {
        Slice slice(globalX, globalY, strideX);
        slice.m_values.resize(strideX * strideY, 0);

        srand(time(0));
        for(auto& value : slice.m_values)
        {
            int randomval = rand() % 10;
            if(randomval > 7)
                value = 1;
        }

        return slice;
    }

    static Slice makeZeroSlice(const size_t globalX, const size_t globalY, const size_t strideX, const size_t strideY)
    {
        Slice slice(globalX, globalY, strideX);
        slice.m_values.resize(strideX * strideY, 0);
        return slice;
    }

    size_t m_globalX;
    size_t m_globalY;
    size_t m_stride;
    std::vector<values_t> m_values;
};

/*!
 * \brief Поле (разбитое на куски)
 */
struct Field
{
    Field():
        m_width{0},
        m_height{0}
    {}

    Field(const size_t width, const size_t height):
        m_width{width},
        m_height{height}
    {}

    size_t m_width;
    size_t m_height;
    std::vector<Slice> m_slices;
};

bool operator==(const Slice& lhs, const Slice& rhs)
{
    return lhs.m_globalX==rhs.m_globalX &&
           lhs.m_globalY==rhs.m_globalY &&
           lhs.m_stride==rhs.m_stride;
}

bool operator==(const Field& lhs, const Field& rhs)
{
    return lhs.m_width==rhs.m_width &&
           lhs.m_height==rhs.m_height;
}

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        void serialize(Archive& ar, Slice& slice, const unsigned int version)
        {
            ar & slice.m_globalX;
            ar & slice.m_globalY;
            ar & slice.m_stride;
            ar & slice.m_values;
        }

        template<class Archive>
        void serialize(Archive& ar, Field& field, const unsigned int version)
        {
            ar & field.m_width;
            ar & field.m_height;
            ar & field.m_slices;
        }
    }
}

#endif // SLICE_H
