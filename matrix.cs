using System;
using System.Net;
/*
 *
 *
 *
 *this class involves mathematical methods used in solving matrix related equations.
 * it handles matrix multiplication, matrix inverses and finding determinant of matrices
 * 
 *
 *this class is used in balancing chemical equation.
 *
 *
 *
 * 
 */
namespace mainCodes
{
    public class matrix
    {
        
        /*
         * methods to solve the matrices.
         * includes:
         *         finding determinant
         *         finding inverse
         *         multiplication
         */
        
        public double[][] MultiplyMatrices(double[][] a, double[][] b) { 
            //multiplies two matrices together
            var m1 = a.Length; //stores number of rows of first matrix in m1
            var n1 = a[0].Length; //stores number of columns of first matrix in n1
            var n2 = b[0].Length; //stores number of rows of second matrix in n2
            double[][] c = new double[m1][]; //creates empty matrix of size m1xn2
            for (var i = 0; i<m1; i++)
            {
                c[i] = new double[n2];
            }
            for (var i = 0; i < m1; i++) { //for loop iterates for each row of first matrix
                for (var j = 0; j < n2; j++) { //for loop iterates for each column of second matrix
                    for (var k = 0; k < n1; k++) { //for loop iterates for each column of first matrix
                        c[i][j] += a[i][k] * b[k][j];
                        //finds i, jth entry of product of the matrices by taking dot product of ith row of first
                        //matrix and jth column of second matrix
                    } //end for loop
                } //end for loop
            } //end for loop
            return c; //returns product of matrices
        } //end multiplyMatrices

        public double FindDet(double[][] matrix) {
            var l = matrix.Length;
            if (l == 2) {
                double det = (matrix[0][0])*(matrix[1][1]) - (matrix[0][1])*(matrix[1][0]);
                //calculate determinant of 2x2 matrix with a*d - b*c and store in double det
                return det;
            }
            else { //if matrix is larger than 2x2
                double det = 0; //initialize det = 0 to store accumulated value of determinant
                for (int i = 0; i < matrix.Length; i++) { //for loop iterates for each column of input matrix
                    det += Math.Pow(-1,i)*matrix[0][i]*FindDet(FindSubmatrix(matrix,0,i));
                    //adds to running value of determinant by multiplying (-1)^n with the ith entry of the first
                    //row of the input matrix and with the determinant of the submatrix obtained by omitting the
                    //1st row and ith column
                }
                return det;
            }
        }
        private double[][] FindSubmatrix(double[][] matrix, int m, int n) {
            //returns submatrix obtained by omitting mth row and nth column of input matrix
            double[][] submatrix = new double[matrix.Length - 1][];///jagged array, i have to iterate to initialize each second array
            for (int i = 0; i < matrix.Length - 1; i++)
            {
                submatrix[i] = new double[matrix.Length-1];
            }
            //creates empty matrix of size (m-1)x(n-1)
		
            for (int j = 0; j < m; j++) { //for loop iterates for each row before omitted row
                for (int k = 0; k < n; k++) { //for loop iterates for each column before omitted column
                    submatrix[j][k] = matrix[j][k]; //stores submatrix values
                } //end for loop
                for (int k = submatrix.Length; k > n; k--) {
                    //for loop iterates for each column after omitted column, starting from right side of input matrix
                    submatrix[j][k-1] = matrix[j][k]; //stores submatrix values
                } //end for loop
            } //end for loop
		
            for (int j = submatrix.Length; j > m; j--) {
                //for loop iterates for each row after omitted row, starting from the bottom of input matrix
                for (int k = 0; k < n; k++) { //for loop iterates for each column before omitted column
                    submatrix[j-1][k] = matrix[j][k]; //stores submatrix values
                } //end for loop
                for (int k = submatrix.Length; k > n; k--) {
                    //for loop iterates for each column after omitted column, starting from right side of input matrix
                    submatrix[j-1][k-1] = matrix[j][k]; //stores submatrix values
                } //end for loop
            } //end for loop
            return submatrix; //returns submatrix
        } //end findSubmatrix
        
        public double[][] FindInverse(double[][] matrix) { //calculates determinant of input matrix
            double[][] inverse = new double[matrix.Length][];
            for (int i = 0; i < matrix.Length; i++)
            {
                inverse[i] = new double[matrix.Length];
            }
            //creates empty matrix with same dimensions as input matrix
            double det = FindDet(matrix); //finds determinant of input matrix
		
            if (matrix.Length == 1) {
                inverse[0][0] = 1/matrix[0][0]; //calculates inverse of 1x1 matrix by taking inverse of the entry
            }
            else if (matrix.Length == 2) {
                inverse[0][0] = matrix[1][1]/det + 0.0;
                inverse[0][1] = -matrix[0][1]/det + 0.0;
                inverse[1][0] = -matrix[1][0]/det + 0.0;
                inverse[1][1] = matrix[0][0]/det + 0.0;
                //calculates inverse of 2x2 matrix with 1/det*([d, -b], [-c, a])
            }
            else {
                double[][] cofactor = FindCofactors(matrix); //finds cofactor of input matrix
                for (int i = 0; i < inverse.Length; i++) { //for loop iterates for each row of input matrix
                    for (int j = 0; j < inverse.Length; j++) { //for loop iterates for each column of input matrix
                        inverse[i][j] = cofactor[j][i]/det + 0.0;
                        //calculates i, jth entry by calculating adjoint of input matrix and dividing by determinant
                    }
                }
            }
            return inverse;
        }
        private double[][] FindCofactors(double[][] matrix) {
            double[][] cofactors = new double[matrix.Length][];
            for (int i = 0; i < matrix.Length; i++)
            {
                cofactors[i] = new double[matrix.Length];
            }
            for (int i = 0; i < matrix.Length; i++) {
                for (int j = 0; j < matrix.Length; j++) {
                    cofactors[i][j] = Math.Pow(-1, i + j)
                                      *FindDet(FindSubmatrix(matrix, i, j));
                }
            }
            return cofactors;
        }
    }
}