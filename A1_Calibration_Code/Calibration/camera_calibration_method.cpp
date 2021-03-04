/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "camera_calibration.h"
#include "matrix_algo.h"
//#include <Eigen/Eigen>
//#include <Eigen/Dense>
#include <math.h>
#include <cmath>

using namespace easy3d;



/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 *
 * @param points_3d   An array of 3D points.
 * @param points_2d   An array of 2D points.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 *           - fx and fy: the focal length (in our slides, we use 'alpha' and 'beta'),
 *           - cx and cy: the principal point (in our slides, we use 'u0' and 'v0'),
 *           - skew:      the skew factor ('-alpha * cot_theta')
 *           - R:         the 3x3 rotation matrix encoding camera orientation.
 *           - t:         a 3D vector encoding camera location.
 */
double dot_product(std::vector<double> v0, std::vector<double> v1)
{
    if (v0.size() != 3 or v1.size() != 3)
    {
        std::cout << "Dot product not possible";
        return 1;
    }
    return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]);
}

std::vector<double> cross_product(std::vector<double> v0, std::vector<double> v1)
{
    double vx = v0[1] * v1[2] - v0[2] * v1[1];
    double vy = v0[2] * v1[0] - v0[0] * v1[2];
    double vz = v0[0] * v1[1] - v0[1] * v1[0];
    return std::vector<double> {vx, vy, vz};
}


double vector_norm(std::vector<double> v)
{
    if (v.size() != 3)
    {
        std::cout << "Vector norm not possible";
        return 1;
    }
    return sqrt(pow(abs(v[0]), 2) + pow(abs(v[1]), 2)+ pow(abs(v[2]), 2));
}


bool CameraCalibration::calibration(
        const std::vector<vec3>& points_3d,
        const std::vector<vec2>& points_2d,
        float& fx, float& fy,
        float& cx, float& cy,
        float& skew,
        mat3& R,
        vec3& t)
{
    std::cout << std::endl;
    std::cout << "TODO: I am going to implement the calibration() function in the following file:" << std::endl
              << "\t" << __FILE__ << std::endl;
    std::cout << "TODO: After implementing the calibration() function, I will disable all unrelated output ...\n\n";

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    // DV: DONE

    if (points_3d.size() < 6) {
        // Less than 6 points
        std::cout << "Input Validation FAIL: The input contains less than 6 points" << std::endl;
        return false;
    } else if (points_3d.size() != points_2d.size()) {
        // 3D and 2D calibration points count do not match
        std::cout << "Input Validation FAIL: The sizes of the 2D/3D points do not match" << std::endl;
        return false;
    } else {
        // Everything all right
        std::cout << "Input Validation PASS." << std::endl;
    }

    // TODO: construct the P matrix (so P * m = 0).

    // arr will contain the values of the P matrix in a one-dimensional array. The length of this array is 2*n*12
    std::vector<double> arr;

    // Points P (3D) and u and v (2D).
    vec3 pt3;
    vec2 pt2;

    for (size_t i = 0; i < points_3d.size(); ++i) {
        //For every point in the input points list

        // The u_i and v_i 2D-coordinates will be stored in vec2 pt2 and can be accessed by pt2.x and pt2.y
        pt2 = points_2d_[i];

        // The 3D-coordinates are stored in vec3 pt3 and can be accessed by pt3.x, pt3.y, and pt3.z
        // To insert the 3D-coordinates into the one-dimensional STL array, the values are accessed with from pt3.data()
        // to pt3.data() + pt3.size() (which is always 3).

        // First constraint row: [P^T 0^T -u*P^T]
        pt3 = points_3d_[i];
        arr.insert(arr.end(), pt3.data(), pt3.data() + 3); // Now add the 3D coordinates of P
        arr.push_back(1); // P^T has to be 4 dimensional (1x4), so add a 1.

        arr.insert(arr.end(), 4, 0); // insert 0-vector

        pt3 = - pt2.x * points_3d_[i];
        arr.insert(arr.end(), pt3.data(), pt3.data() + 3); // Add the 3D coordinates of P
        arr.push_back( - pt2.x );

        // Second constraint row: [0^T P^T -v*P^T]
        arr.insert(arr.end(), 4, 0); // insert 0-vector

        pt3 = points_3d_[i];
        arr.insert(arr.end(), pt3.data(), pt3.data() + 3); // Now add the 3D coordinates of P
        arr.push_back(1); // P^T has to be 4 dimensional (1x4), so add a 1.

        pt3 = - pt2.y * points_3d_[i];
        arr.insert(arr.end(), pt3.data(), pt3.data() + 3); // Add the 3D coordinates of P
        arr.push_back( - pt2.y );
    }

    // Create Matrix P of size (mxn = 2*n_points x 12)
    const int m = 2 * points_3d.size();
    const int n = 12;

    Matrix<double> P(m, n, arr.data());
    std::cout << "Constructed Matrix P: " << P << std::endl;

    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.
    Matrix<double> U(m, m, 0.0);   // initialized with 0s
    Matrix<double> S(m, n, 0.0);   // initialized with 0s
    Matrix<double> V(n, n, 0.0);   // initialized with 0s

    // Single Value Decomposition into U, S, and V
    svd_decompose(P, U, S, V);

    // We get m by taking the last column of V.
    const auto mvec = V.get_column(n - 1);

    // Reshape vector m to matrix M
    Matrix<double> M(3, 4, mvec.data());
    std::cout << "M " << M << std::endl;

    // TODO: hier moet matrix M nog in K en [R T] worden opgesplitst...
    //Eigen::Matrix<Eigen::Vector3d::Double, Eigen::Dynamic, Eigen::Dynamic> mat(3, 4);
    //auto QR = mat.householderQr()

    // Remove checks




    // TODO: extract intrinsic parameters from M.
    // split m into A and b
    std::vector<double> a1{M[0][0], M[0][1], M[0][2]};
    std::vector<double> a2{M[1][0], M[1][1], M[1][2]};
    std::vector<double> a3{M[2][0], M[2][1], M[2][2]};
    std::vector<double> b{M[0][3], M[1][3], M[2][3]};

    // calculate intrinsic parameters
    double rho = 1 / vector_norm(a3);
    cx = pow(rho, 2) * dot_product(a1, a3);
    cy = pow(rho, 2) * dot_product(a2, a3);
    double theta = acos(-(dot_product(cross_product(a1, a3), cross_product(a2, a3)) / (vector_norm(cross_product(a1,a3)) * vector_norm(cross_product(a2,a3)))));
    double alpha = pow(rho, 2) * vector_norm(cross_product(a1,a3)) * sin(theta);
    double beta = pow(rho, 2) * vector_norm(cross_product(a2,a3)) * sin(theta);

    // TODO: extract extrinsic parameters from M.

    // TODO: uncomment the line below to return true when testing your algorithm and in you final submission.
    //return false;



    // TODO: The following code is just an example showing you SVD decomposition, matrix inversion, and some related.
    // TODO: Delete the code below (or change "#if 1" in the first line to "#if 0") in you final submission.
#if 1
    std::cout << "[Liangliang:] Camera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a Matrix data structure for storing matrices of arbitrary\n"
                 "\tsizes (see matrix.h). I also wrote the example code to show you how to:\n"
                 "\t\t- use the dynamic 1D array data structure 'std::vector' from the standard C++ library;\n"
                 "\t\t  The points (both 3D and 2D) are stored in such arrays;\n"
                 "\t\t- use the template matrix class (which can have an arbitrary size);\n"
                 "\t\t- compute the SVD of a matrix;\n"
                 "\t\t- compute the inverse of a matrix;\n"
                 "\t\t- compute the transpose of a matrix.\n"
                 "\tThe following are just the output of these examples. You should delete ALL unrelated code and\n"
                 "\tavoid unnecessary output in you final submission.\n\n";

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can do append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
//    std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
//    array.push_back(5); // append 5 to the array (so the size will increase by 1).
//    array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).
//
//    // To access its values
//    for (int i=0; i<array.size(); ++i)
//        std::cout << array[i] << " ";  // use 'array[i]' to access its i-th element.
//    std::cout << std::endl;

    // Define an m-by-n double valued matrix.
    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    //const int m = 6, n = 5;
    //Matrix<double> A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
    //std::cout << "M: \n" << A << std::endl;

    //Matrix<double> U(m, m, 0.0);   // initialized with 0s
    //Matrix<double> S(m, n, 0.0);   // initialized with 0s
    //Matrix<double> V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    //svd_decompose(A, U, S, V);

    // Now let's check if the SVD result is correct

    // Check 1: U is orthogonal, so U * U^T must be identity
//    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;
//
//    // Check 2: V is orthogonal, so V * V^T must be identity
//    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;
//
//    // Check 3: S must be a diagonal matrix
//    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
    //std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

    // Define a 5 by 5 square matrix and compute its inverse.
    //Matrix<double> B(5, 5, array.data());    // Here I use part of the above array to initialize B
    // Compute its inverse
    //Matrix<double> invB(5, 5);
    //inverse(B, invB);
    // Let's check if the inverse is correct
    //std::cout << "B * invB: \n" << B * invB << std::endl;

    return false;
    // TODO: delete the above code in you final submission (which are just examples).
#endif
}