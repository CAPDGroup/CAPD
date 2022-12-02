from sage.all import *
import sys

def matrix_to_str(matrix):
    rows = [('{' + ', '.join(map(str,r)) + '}') for r in matrix]
    return '{' + ', '.join(rows) + '}'

m = random_matrix(ZZ, 15, 15)
s = m.smith_form()[0]


print '{\nint matrix' + '[{}][{}] = {};'.format(m.nrows(), m.ncols(), matrix_to_str(m)) +'\nmatrices_matrix.push_back(Matrix<int, 0, 0>(matrix));\n}'

print '{\nint smithMatrix' + '[{}][{}] = {};'.format(s.nrows(), s.ncols(), matrix_to_str(s)) +'\nmatrices_smithMatrix.push_back(Matrix<int, 0, 0>(smithMatrix));\n}'
