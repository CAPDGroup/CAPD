///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdexcept>

template<typename size_type, size_type max_dimension>
class C2_Jet_Indices_Buffer
{
public:
    static size_type constexpr compute_index(size_type dim, size_type j, size_type c)
    {
        if (dim > max_dimension)
        {
            return compute_index_unoptimized(dim, j, c);
        }

        return g_buffer.impl[dim][j][c];
    }

private:
    struct Buffer
    {
        size_type impl[max_dimension+1][max_dimension][max_dimension] {};
    };

    static size_type constexpr compute_index_unoptimized(size_type dim, size_type j, size_type c)
    {
        if (dim > j)
        {
            if (j <= c)
            {
                return 1+c + ((j+1)*(dim*2-j))/2;
            }
            else
            {
                return 1+j + ((c+1)*(dim*2-c))/2;
            }
        }

        throw std::logic_error("Dimension must be greater than j!");
    }

    static constexpr Buffer initialize()
    {
        Buffer ret {};
        for (size_type dim = 0; dim <= max_dimension; ++dim)
        {
            for (size_type j = 0; j < max_dimension; ++j)
            {
                for (size_type c = 0; c < max_dimension; ++c)
                {
                    if (j < dim)
                    {
                        ret.impl[dim][j][c] = compute_index_unoptimized(dim, j, c);
                    }
                }
            }
        }

        return ret;
    }

    static constexpr Buffer g_buffer
    {
        initialize()
    };
};
