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
    std::cout << "calibration() function is executing with the following file:" << std::endl
              << "\t" << __FILE__ << std::endl;

    // check if input is valid
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

    // construct the P matrix
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

    // solve for M using SVD decomposition.
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

    // extract intrinsic parameters from M.
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
    fx = pow(rho, 2) * vector_norm(cross_product(a1,a3)) * sin(theta);
    fy = pow(rho, 2) * vector_norm(cross_product(a2,a3)) * sin(theta);

    // extract extrinsic parameters from M.
    // Use the A and b vectors from intrinsic parameter calculation above
    // Calculate extrinsic parameters. First R matrix's constituent parts, r1, r2, and r3:
    std::vector<double> r_1 = cross_product(a2, a3) / vector_norm(cross_product(a2, a3));
    std::vector<double> r_3 = rho * a3;
    std::vector<double> r_2 = cross_product(r_3, r_1);
    // Fill matrix R (defined at start of code) with constituent parts
    R(0, 0) = r_1[0];
    R(0, 1) = r_1[1];
    R(0, 2) = r_1[2];

    R(1, 0) = r_2[0];
    R(1, 1) = r_2[1];
    R(1, 2) = r_2[2];

    R(2, 0) = r_3[0];
    R(2, 1) = r_3[1];
    R(2, 2) = r_3[2];
    // print matrix R
    std::cout << "R: \n" << R << std::endl;

    // Must now calculate extrinsic parameter t, where t = rho * K^-1 * b
    // Intrinsic Matrix K is needed. And for K, skew is needed, where skew = cot(theta):
    skew = -fx * (1 / tan(theta));

    // Define matrix K and make array of what belongs in K:
    std::vector<double> K_array = {fx, skew, cx, 0, fy, cy, 0, 0, 1};
    Matrix<double> K(3, 3, K_array.data());
    std::cout << "K: \n" << K << std::endl;

    // Compute matrix K's inverse, invK:
    Matrix<double> invK(3, 3);
    inverse(K, invK);
    // Check to see if inverse is correct:
    std::cout << "K * invK: \n" << K * invK << std::endl;

    // Calculate t vector with rho, inverse of matrix K, and b
    Matrix<double> T(3, 1);
    Matrix<double> B(3, 1,
        std::vector<double> {M(0, 3), M(1, 3), M(2, 3)}.data());
    T = (rho * invK * B);
    t = {float(T[0][0]), float(T[1][0]), float(T[2][0])};

    // Print extrinsic parameters:
    std::cout << "Extrinsic parameters: \n"
        << "R: \n" << R
        << "T: \n" << t;

    return false;

}