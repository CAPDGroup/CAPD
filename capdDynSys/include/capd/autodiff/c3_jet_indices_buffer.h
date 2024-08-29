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
        if (c <= k)
        {
            size_t ret = k - c;
            if (dim > max_dimension)
            {
                ret +=  compute_internal(dim, j, c);
            }
            else
            {
                ret += g_buffer.impl[dim][j][c];
            }
            return ret;
        }
        
        throw std::logic_error("c must be less or equal k!");
    }

private:
    struct Buffer
    {
        size_type impl[max_dimension+1][max_dimension][max_dimension] {};
    };
        
    static size_type constexpr compute_internal(size_type dim, size_type j, size_type c)
    {
        if (j <= c)
        {
            return 
                (
                    (1+dim)*(2+dim) +
                    (j*(j-1)*(j-2))/3 +
                    j*dim*(dim-j+2) +
                    (j-c)*(c+j-2*dim-1)
                ) /2;
        }

        throw std::logic_error("capd::autodiff::index(size_type,size_type,size_type) - indices are not ordered");

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
                    if (j <= c)
                    {
                        ret.impl[dim][j][c] = compute_internal(dim, j, c);
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
