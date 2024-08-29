///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdexcept>

template<typename size_type, size_type max_dimension>
class C3_Jet_Indices_Buffer
{
public:
    static size_type constexpr compute_index(size_type dim, size_type j, size_type c, size_type k)
    {
        if (dim > max_dimension)
        {
            return compute_index_unoptimized(dim, j, c, k);
        }

        return g_buffer.impl[dim][j][c][k];
    }

private:
    struct Buffer
    {
        size_type impl[max_dimension+1][max_dimension][max_dimension][max_dimension] {};
    };
        
    static size_type constexpr compute_index_unoptimized(size_type dim, size_type j, size_type c, size_type k)
    {
        // assume j<=c<=k
        if (j <= c && c <= k)
        {
            return k - c +
                (
                    (1+dim)*(2+dim) +
                    (j*(j-1)*(j-2))/3 +
                    j*dim*(dim-j+2) +
                    (j-c)*(c+j-2*dim-1)
                ) /2;
        }

        throw std::logic_error("Unexpected args!");
        // throw std::logic_error("capd::autodiff::index(size_type,size_type,size_type) - indices are not ordered");

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
                    for (size_type k = 0; k < max_dimension; ++k)
                    {
                        if (j <= c && c <= k)
                        {
                            ret.impl[dim][j][c][k] = compute_index_unoptimized(dim, j, c, k);
                        }
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
