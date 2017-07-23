#include <iostream>
#include <ctime>
#include <Eigen/Core>
#include <Eigen/Dense>

#define MATRIX_SIZE 50

using namespace std;

int main() {
  // Define a 2x3 float matrix 
  Eigen::Matrix<float, 2, 3> matrix_23;
  // Type Eigen::Vector3d is defined as Eigen::Matrix<double, 3, 1>
  Eigen::Vector3d vector_3d;
  // Type Eigen::Matrix3d is defined as Eigen::Matrix<double, 3, 3>
  Eigen::Matrix3d matrix_33 = Eigen::Matrix3d::Zero(); // Zero initialized
  // Use dynamic matrix if the size is not sure
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_dynamic;
  // The simpler way
  Eigen::MatrixXd matrix_x;
  
  // Input
  matrix_23 << 1, 2, 3, 4, 5, 6;
  // Output
  cout << matrix_23 << endl;
  
  // Use () operator to get elements in the matrix
  for (int i = 0; i < 1; i++) {
    for (int j = 0; j < 2; j++) {
      cout << matrix_23 << "(" << i << "," << j << ") = "
        << matrix_23(i,j) << endl;
    }
  }

  // Input for vectors
  v_3d << 3, 2, 1;
  // Multiply matrix by vector,
  // but it will throw an error if they are two different types
  // Eigen::Matrix<double, 2, 1> result_wrong_type = matrix_23 * v_3d;

  // We shall convert it obviously
  Eigen::Matrix<double, 2, 1> result = matrix_23.cast<double>() * v_3d;
  cout << "result = \n" << result << endl;

  // We shall care for the dimension of matrix
  // Eigen::Matrix<double, 2, 3> result_wrong_dimension = matrix_23.cast<double>() * v_3d;

  // Some operation of matrix
  matrix_33 = Eigen::Matrix3d::Random();
  cout << "matrix_33 = \n" << matrix_33 << endl;

  cout << "matrix_33.transpose() = \n" << matrix_33.transpose() << endl;
  cout << "matrix_33.sum() = \n" << matrix_33.sum() << endl;
  cout << "matrix_33.trace() = \n" << matrix_33.trance() << endl;
  cout << "10 * matrix_33 = \n" << 10 * matrix_33 << endl;
  cout << "matrix_33.inverse() = \n" << matrix_33.inverse() << endl;
  cout << "matrix_33.determinant() = \n" << matrix_33.determinant() << endl;
  
  // Eigenvector
  // Real symmetric matrix can guarantee the success of diagonalization
  Eigen::SelfAdjoinEigenSolver<Eigen::Matrix3d> eigen_solver(matrix_33.transpose() * matrix_33);
  cout << "Eigen values = " << eigen_solver.eigenvalues() << endl;
  cout << ""
}


