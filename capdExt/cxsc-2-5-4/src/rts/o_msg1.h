/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: o_msg1.h,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/*------------------------------------------------------*/
/* Notes:                                               */
/*                                                      */
/* -Data initialization is used by "b_glbl.c".          */
/*                                                      */
/* -Data structure e_mtyp is defined in "o_type.h".     */
/*                                                      */
/* -All msgid numbers (1.component) must be different.  */
/*  1 <= msgid < E_TRES                                 */
/*                                                      */
/* -Data initializations may be ordered arbitrarily.    */
/*                                                      */
/* -There must be a final comma for each item.          */
/*                                                      */
/* -"\" is only allowed in conjunction with "r" or "n". */
/*                                                      */
/* - "\r" advances to the next line if '\r' is not the  */
/*   first character and displays e_hlen Blanks if      */
/*   '\r' is not the first or last character.           */
/*   If '\r' is the first character, the current line   */
/*   is continued and no indentation is done.           */
/*                                                      */
/* - "\n" advances to the next line and displays        */
/*   e_head if '\n' is not the last character.          */
/*                                                      */
/*------------------------------------------------------*/

{    1,  "Division by zero\n" },
{    2,  "\r0/0\n" },
{    3,  "Invalid operation : " },
{    4,  "\rinfinity/infinity\n" },
{    5,  "\rSignaling NaN as operand\n" },
{    6,  "Overflow occurred\n" },
{    7,  "Underflow occurred\n" },
{    8,  "Inexact\n" },
{    9,  "\rinfinity-infinity\n" },
{   10,  "\r0*infinity\n" },
{   11,  "Allocation failed : " },
{   12,  "Evaluation error possibly caused by invalid argument.\n" },
{   13,  "Unexpected infinity operand.\n" },
{   14,  "Unexpected quiet NaN operand.\n" },
{   15,  "Range of integer data type exceeded.\n" },
{   16,  "Error in I/O operation : " },
{   17,  "\rNo device assigned.\n" },
{   18,  "\rDevice not opened for reading.\n" },
{   19,  "\rDevice not a TEXT device.\n" },
{   20,  "\rUnexpected End-Of-File.\n" },
{   21,  "\rInvalid syntax of integer value.\n" },
{   22,  "Non-positive modulo value.\n" },
{   23,  "Index out of range.\n" },
{   24,  "Unexpected NULL pointer.\n" },
{   25,  "\r------ Processing aborted ------\n" },
{   26,  "One-dimensional dynamic array expected.\n" },
{   27,  "Equal length of dynamic vectors expected.\n" },
{   28,  "Substring destination array shorter than required.\n" },
{   29,  "Internal buffer too small : " },
{   30,  "\rFilename too long\n" },
{   31,  "\rUnable to open file for reading\n" },
{   32,  "\rUnable to open file for writing\n" },
{   33,  "\rStandard I/O must not be used for binary I/O.\n" },
{   34,  "\rDevice not opened for writing.\n" },
{   35,  "\rError writing data to file.\n" },
{   36,  "\rDevice not a binary device.\n" },
{   37,  "\rMissing command line argument.\n" },
{   38,  "\rNo filename has previously been assigned.\n" },
{   39,  "\rDynamic mantissa too long.\n" },
{   40,  "\rDotprecision variable not allocated.\n" },
{   41,  "\rTrap handler not stacked.\n" },
{   42,  "\rDynamic array not allocated.\n" },
{   43,  "Undefined device specification.\n" },
{   44,  "Invalid width of output field.\n" },
{   45,  "\rInvalid filename.\n" },
{   46,  "\rDynamic object not allocated.\n" },
{   47,  "Mantissa out of range.\n" },
{   48,  "Exponent too large (infinity returned).\n" },
{   49,  "Mantissa bits lost on generating denormalized number.\n" },
{   50,  "Exponent too small (zero returned).\n" },
{   51,  "\rInvalid read/write mode.\n" },
{   52,  "\rInvalid syntax of hexadecimal value.\n" },
{   53,  "\rUnexpected End-Of-Line.\n" },
{   54,  "\rDynamic string not allocated.\n" },
{   55,  "\rReading a dynamic string.\n" },
{   56,  "\rConverting string to real.\n" },
{   57,  "\rConverting real to string.\n" },
{   58,  "\rInvalid syntax of real value.\n" },
{   59,  "\r Inexact conversion of decimal constant.\n" },
{   60,  "\r Inexact conversion of decimal input data.\n" },
{   61,  "\rMissing variable name.\n" },
{   62,  "\rEmpty string.\n" },
{   63,  "\rNo digits in string.\n" },
{   64,  "Exponent range restricted.\n" },
{   65,  "\rDynamic variable.\n" },
{   66,  "\rConverting string to interval.\n" },
{   67,  "Mismatching index ranges.\n" },
{   68,  "Message : " },
{   69,  "\n" },
{   70,  "\rInvalid data in tenbyte operation.\n" },



/*------------------------------------------------------*/
/* These codes are explicitly generated by PASCAL-XSC.  */
/*------------------------------------------------------*/
{  400,  "Summation of vectors with different lengths.\n" },
{  401,  "Scalar product of vectors with different lengths.\n" },
{  402,  "Function call with vectors of different lengths.\n" },
{  403,  "Mismatching inner lengths in a matrix-vector product.\n" },
{  404,  "Summation of matrices with different row lengths.\n" },
{  405,  "Summation of matrices with different column lengths.\n" },
{  406,  "Function call with matrices of different row lengths.\n" },
{  407,  "Function call with matrices of different column lengths.\n" },
{  408,  "Mismatching inner lengths of arguments in a function call.\n" },
{  409,  "Mismatching inner lengths in a matrix-matrix product.\n" },

/*------------------------------------------------------*/
/* These codes must be even multiples of E_TRES=256     */
/* within the range                                     */
/*           (2*256=) 512<=msgid<=32256 (=126*256).     */
/*------------------------------------------------------*/

{  512,  "\r left operand : " },                  /* E_TEXT(1)    */
{ 1024,  "\rright operand : " },                  /* E_TEXT(2)    */
{ 1536,  "\r       result : " },                  /* E_TEXT(3)    */
{ 2048,  "\r        index : " },                  /* E_TEXT(4)    */
{ 2560,  "\r  lower bound : " },                  /* E_TEXT(5)    */
{ 3072,  "\r  upper bound : " },                  /* E_TEXT(6)    */
{ 3584,  "\r     argument : " },                  /* E_TEXT(7)    */
{ 4096,  "\r    file name : " },                  /* E_TEXT(8)    */
{ 4608,  "\rfile variable : " },                  /* E_TEXT(9)    */
{ 5120,  "\r   input data : " },                  /* E_TEXT(10)   */
{ 5632,  "\r    dimension : " },                  /* E_TEXT(11)   */
{ 6144,  "\rvector length : " },                  /* E_TEXT(12)   */
{ 6656,  "\r   bit number : " },                  /* E_TEXT(13)   */
{ 7168,  "\r  set element : " },                  /* E_TEXT(14)   */
{ 7680,  "\rstring length : " },                  /* E_TEXT(15)   */
{ 8192,  "\r   error code : " },                  /* E_TEXT(16)   */
{ 8704,  "\r        basis : " },                  /* E_TEXT(17)   */
{ 9216,  "\r     mantissa : " },                  /* E_TEXT(18)   */
{ 9728,  "\r     exponent : " },                  /* E_TEXT(19)   */





