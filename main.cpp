#include <amino.h>
#include <amino/mat.h>
#include <cmath>

void
aafeq(const char *name, double a, double b, double tol)
{
    if (!aa_feq(a, b, tol)) {
        fprintf(stderr, "FAILED: %s\n", name);
        fprintf(stderr, "a: %f, b: %f\n", a, b);
        abort();
    }
}

static inline void
admeq( const char *name, const double *A, const double *B, double tol, size_t size)
{
    aafeq( name, aa_la_ssd(size,A,B), 0, tol );
}

static inline void
mat_print_raw(const double *data, const int rows, const int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << data[rows*j+i] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

static inline void
mat_print_pretty(const double *data, const int rows, const int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double x = data[rows*j+i];
            if (fabs(x) < 1e-6) {
                std::cout << 0 << "  ";
            } else {
                std::cout << x << "  ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

static inline void
array_print(const double *data, const int n)
{
    for (int i = 0; i < n; i++) {
        std::cout << data[i];
        std::cout << "\n";
    }
    std::cout << "\n";
}

/* ik7dof takes the DH parameter link lengths, desired position of tool,
   desired orientation of tool, and returns the joint angles needed */
void ik7dof(const double &d1,
            const double &d3,
            const double &d5,
            const double &d7,
            const double p_T_data[3],
            const struct amino::Quat &qu_T,
            const double &elbow_sign_param,
            const double &elbow_ang_param,
            const double &arm_sign_param,
            const double &wrist_sign_param) 
{
    std::cout << "Got into ik7dof() successfully!\n\n";

    /* Following  "Analytical Inverse Kinematics and Self-Motion 
    Application for 7-DOF Redundant Manipulator" - M. Gong et al.*/

    /* Verifying input parameter of elbow joint sign */
    double abs_val = fabs(elbow_sign_param);
    if (abs_val - 1 != 0) {
        std::cout << "ERROR: Elbow joint sign parameter needs to be either -1 or 1.\n\n";
        assert(abs_val - 1 == 0);
    }

    /* Verifying input parameter of elbow self motion angle*/
    if (elbow_ang_param < -M_PI || elbow_ang_param > M_PI) {
        std::cout << "ERROR: Elbow self-motion angle parameter needs to be from -pi to pi, inclusive.\n\n";
        assert(elbow_ang_param >= -M_PI && elbow_ang_param <= M_PI);
    }

    /* Verifying input parameter of arm sign */
    double abs_val_2 = fabs(arm_sign_param);
    if (abs_val_2 - 1 != 0) {
        std::cout << "ERROR: Arm sign parameter needs to be either -1 or 1.\n\n";
        assert(abs_val_2 - 1 == 0);
    }

    /* Verifying input parameter of wrist sign */
    double abs_val_3 = fabs(wrist_sign_param);
    if (abs_val_3 - 1 != 0) {
        std::cout << "ERROR: Wrist sign parameter needs to be either -1 or 1.\n\n";
        assert(abs_val_3 - 1 == 0);
    }

    /* See Figure 1 in M. Gong et al. for clarification on the DH parameters */
    double alpha1 = -M_PI_2; // radians
    double alpha2 = M_PI_2;
    double alpha3 = -M_PI_2;
    double alpha4 = M_PI_2;
    double alpha5 = -M_PI_2;
    double alpha6 = M_PI_2;
    double alpha7 = 0;
    double a1 = 0; // meter
    double a2 = 0;
    double a3 = 0;
    double a4 = 0;
    double a5 = 0;
    double a6 = 0;
    double a7 = 0;
    double d2 = 0;
    double d4 = 0;
    double d6 = 0;
    double d_BS = d1; // B=base, S=shoulder, E=elbow, W=wrist, T=tool
    double d_SE = d3;
    double d_EW = d5;
    double d_WT = d7;

    /* Position of shoulder */
    struct amino::Vec3 p_S{0, 0, d_BS};
    std::cout << "Position of S:\n";
    array_print(p_S.data, 3);

    /* User-assigned position of tool */
    struct amino::Vec3 p_T{p_T_data};
    std::cout << "Position of T:\n";
    array_print(p_T.data, 3);

    // We have orientation of tool as qu_T

    /* User-assigned pose/transformation of tool */
    struct amino::QuatTran qutr_T{qu_T, p_T};

    /* Equation 2 of M. Gong et al. */
    double temp_data[] = {0, 0, -d_WT};
    double p_W_data[3];
    qutr_T.transform(temp_data, p_W_data);
    struct amino::Vec3 p_W{p_W_data};
    std::cout << "Position of W:\n";
    array_print(p_W.data, 3);

    /* Check limit of d_SW */
    double v_SW_data[3];
    aa_la_vsub(3, p_W.data, p_S.data, v_SW_data);
    struct amino::Vec3 v_SW{v_SW_data};
    std::cout << "Displacement from S to W:\n";
    array_print(v_SW.data, 3);
    double d_SW = aa_la_norm(3, v_SW.data);
    double difference = d_SW - (d_SE + d_EW);
    if (difference >= -0.02 && difference <= 0) {
        std::cout << "WARNING: d_SW approximately equal d_SE + d_EW --> pose is near workspace limit!\n";
    } else if (difference > 0) {
        std::cout << "ERROR: d_SW > d_SE + d_EW\n";
        assert(d_SW - (d_SE + d_EW) <= 0);
    }

    /* Find joint 4 */
    double cos_SEW  = (pow(d_SE,2)+pow(d_EW,2)-pow(d_SW,2)) / (2*d_SE*d_EW);
    double q4 = elbow_sign_param * (M_PI - acos(cos_SEW));
    // while () // joint limits
    std::cout << "Joint q4 = " << q4 << "\n\n";


    /* Columns for coordinate system sigma-D when psi (elbow self-motion angle) = 0 */
    double X_sigma_data[3];
    memcpy(X_sigma_data, v_SW_data, 3*sizeof(double));
    aa_la_normalize(3, X_sigma_data);
    // std::cout << "X_sigma_data:\n";
    // array_print(X_sigma_data, 3);

    double Z_0_data[3] = {0, 0, 1};
    double Y_sigma_data[3];
    aa_la_cross(Z_0_data, X_sigma_data, Y_sigma_data);
    aa_la_normalize(3, Y_sigma_data);
    double temp = aa_la_norm(3, Y_sigma_data);
    if (temp <= AA_EPSILON) { // special case when X_sigma and Z_0 are parallel
        double Y_0_data[3] = {0, 1, 0};
        memcpy(Y_sigma_data, Y_0_data, 3*sizeof(double));
    }
    // std::cout << "Y_sigma_data:\n";
    // array_print(Y_sigma_data, 3);

    double Z_sigma_data[3];
    aa_la_cross(X_sigma_data, Y_sigma_data, Z_sigma_data);
    aa_la_normalize(3, Z_sigma_data);
    // std::cout << "Z_sigma_data:\n";
    // array_print(Z_sigma_data, 3);

    struct amino::RotMat R_sigma_D_psi_0{X_sigma_data[0], Y_sigma_data[0], Z_sigma_data[0],
                                         X_sigma_data[1], Y_sigma_data[1], Z_sigma_data[1],
                                         X_sigma_data[2], Y_sigma_data[2], Z_sigma_data[2]}; // aa_tf_rotmat
    struct amino::Quat qu_sigma_D_psi_0{R_sigma_D_psi_0}; // aa_tf_quat


    /* Find coordinate system sigma-D with given psi (elbow self-motion angle) */
    /* One way */
    // double v_SW_data_normed[3];
    // memcpy(v_SW_data_normed, v_SW_data, 3*sizeof(double));
    // aa_la_normalize(3, v_SW_data_normed);
    // struct amino::AxisAngle axang_SW_psi{v_SW_data_normed, elbow_ang_param}; // aa_tf_axang
    // struct amino::Quat qu_SW_psi{axang_SW_psi}; // aa_tf_quat
    // double qu_sigma_D_data_1[4];
    // aa_tf_qmulnorm(qu_SW_psi.data, qu_sigma_D_psi_0.data, qu_sigma_D_data_1); // with angle-axis, order of matmul not important
    // std::cout << "Checking 1st way result:\n";
    // array_print(qu_sigma_D_data_1, 4);
    /* Another way */
    struct amino::XAngle x_angle_psi{elbow_ang_param}; 
    struct amino::Quat qu_x_psi{x_angle_psi}; // aa_tf_quat
    double qu_sigma_D_data_2[4];
    aa_tf_qmulnorm(qu_sigma_D_psi_0.data, qu_x_psi.data, qu_sigma_D_data_2);
    std::cout << "Checking 2nd way result for sigma_D quaternion:\n";
    array_print(qu_sigma_D_data_2, 4);


    /* Find coordinate system sigma-0 */
    double p_W_data_2[3] = {d_EW*sin(q4), 0, d_BS+d_SE+d_EW*cos(q4)};
    std::cout << "Position of W when calculating sigma_0:\n";
    array_print(p_W_data_2, 3);

    double v_SW_data_2[3];
    aa_la_vsub(3, p_W_data_2, p_S.data, v_SW_data_2);
    std::cout << "Displacement from S to W when calculating sigma_0:\n";
    array_print(v_SW_data_2, 3);

    double X_sigma_data_2[3];
    memcpy(X_sigma_data_2, v_SW_data_2, 3*sizeof(double));
    aa_la_normalize(3, X_sigma_data_2);
    // std::cout << "X_sigma_data_2:\n";
    // array_print(X_sigma_data_2, 3);

    // double Z_0_data[3] = {0, 0, 1};
    double Y_sigma_data_2[3];
    aa_la_cross(Z_0_data, X_sigma_data_2, Y_sigma_data_2);
    aa_la_normalize(3, Y_sigma_data_2);
    double temp_2 = aa_la_norm(3, Y_sigma_data_2);
    if (temp_2 <= AA_EPSILON) { // special case when X_sigma and Z_0 are parallel
        double Y_0_data[3] = {0, 1, 0};
        memcpy(Y_sigma_data_2, Y_0_data, 3*sizeof(double));
    }
    // std::cout << "Y_sigma_data_2:\n";
    // array_print(Y_sigma_data_2, 3);

    double Z_sigma_data_2[3];
    aa_la_cross(X_sigma_data_2, Y_sigma_data_2, Z_sigma_data_2);
    aa_la_normalize(3, Z_sigma_data_2);
    // std::cout << "Z_sigma_data_2:\n";
    // array_print(Z_sigma_data_2, 3);

    struct amino::RotMat R_sigma_0{X_sigma_data_2[0], Y_sigma_data_2[0], Z_sigma_data_2[0],
                                   X_sigma_data_2[1], Y_sigma_data_2[1], Z_sigma_data_2[1],
                                   X_sigma_data_2[2], Y_sigma_data_2[2], Z_sigma_data_2[2]}; // aa_tf_rotmat
    struct amino::Quat qu_sigma_0{R_sigma_0}; // aa_tf_quat
    // array_print(qu_sigma_0.data, 4);


    /* Find rotation of shoulder spherical joint */
    // double qu_sigma_0_inv_data[4];
    // aa_tf_qinv(qu_sigma_0.data, qu_sigma_0_inv_data);
    // double qu_S[4];
    // aa_tf_qmulnorm(qu_sigma_D_data_2, qu_sigma_0_inv_data, qu_S);
    // array_print(qu_S, 4);
    double qu_S_data[4];
    aa_tf_qmulc(qu_sigma_D_data_2, qu_sigma_0.data, qu_S_data);
    // array_print(qu_S_data, 4);
    double R_S_data[9];
    aa_tf_quat2rotmat(qu_S_data, R_S_data);
    std::cout << "Matrix R_S = R_0_3:\n";
    mat_print_raw(R_S_data, 3, 3);

    /* Calculate joint 2 */
    double q2 = arm_sign_param * acos(R_S_data[8]);
    std::cout << "Joint q2 = " << q2 << "\n\n";

    /* Calculate joint 1*/
    double q1 = atan2( arm_sign_param*R_S_data[7] , arm_sign_param*R_S_data[6] );
    std::cout << "Joint q1 = " << q1 << "\n\n";

    /* Calculate joint 3*/
    double q3 = atan2( arm_sign_param*R_S_data[5] , -arm_sign_param*R_S_data[2] );
    std::cout << "Joint q3 = " << q3 << "\n\n";


    /* Find rotation of last three joints (5, 6, 7) combined, qu_W */
    struct amino::YAngle y_angle_q4{q4};
    struct amino::Quat qu_y_q4{y_angle_q4}; // aa_tf_quat
    double qu_0_4_data[4];
    aa_tf_qmulnorm(qu_S_data, qu_y_q4.data, qu_0_4_data);
    /* One way*/
    // double qu_W_data[4];
    // aa_tf_qmulc(qu_T.data, qu_0_4_data, qu_W_data); // qu_T means qu_0_7, rotation from 0 to 7
    // std::cout << "Quaternion of rot from 4 to 7:\n";
    // array_print(qu_W_data, 4);
    /* Another way, matches frame transitions better on paper */
    double qu_W_data[4];
    aa_tf_qcmul(qu_0_4_data, qu_T.data, qu_W_data); // qu_T means qu_0_7, rotation from 0 to 7
    std::cout << "Quaternion of R_4_7, 2nd way:\n";
    array_print(qu_W_data, 4);

    double R_W_data[9];
    aa_tf_quat2rotmat(qu_W_data, R_W_data);
    std::cout << "Matrix R_W = R_4_7:\n";
    mat_print_raw(R_W_data, 3, 3);

    /* Calculate joint 6 */
    double q6 = wrist_sign_param * acos(R_W_data[8]);
    std::cout << "Joint q6 = " << q6 << "\n\n";

    /* Calculate joint 5 */
    double q5 = atan2( wrist_sign_param*R_W_data[7] , wrist_sign_param*R_W_data[6] );
    std::cout << "Joint q5 = " << q5 << "\n\n";

    /* Calculate joint 7 */
    double q7 = atan2( wrist_sign_param*R_W_data[5] , -wrist_sign_param*R_W_data[2] );
    std::cout << "Joint q7 = " << q7 << "\n\n";
    // q7 = 0;


    // /* Check with forward kinematics, using modified (proximal) DH per M. Gong et al. */
    // double T_0_1[12];
    // aa_tf_dhprox2tfmat(0, 0, d1, q1, T_0_1);
    // double T_1_2[12];
    // aa_tf_dhprox2tfmat(alpha1, a1, d2, q2, T_1_2);
    // double T_2_3[12];
    // aa_tf_dhprox2tfmat(alpha2, a2, d3, q3, T_2_3);
    // double T_3_4[12];
    // aa_tf_dhprox2tfmat(alpha3, a3, d4, q4, T_3_4);
    // double T_4_5[12];
    // aa_tf_dhprox2tfmat(alpha4, a4, d5, q5, T_4_5);
    // double T_5_6[12];
    // aa_tf_dhprox2tfmat(alpha5, a5, d6, q6, T_5_6);
    // double T_6_7[12];
    // aa_tf_dhprox2tfmat(alpha6, a6, d7, q7, T_6_7);   
    // double T_0_2[12];  
    // aa_tf_12chain(T_0_1, T_1_2, T_0_2);
    // double T_0_3[12];  
    // aa_tf_12chain(T_0_2, T_2_3, T_0_3);
    // double T_0_4[12];  
    // aa_tf_12chain(T_0_3, T_3_4, T_0_4);
    // double T_0_5[12];  
    // aa_tf_12chain(T_0_4, T_4_5, T_0_5);
    // double T_0_6[12];  
    // aa_tf_12chain(T_0_5, T_5_6, T_0_6);
    // double T_0_7[12];  
    // aa_tf_12chain(T_0_6, T_6_7, T_0_7);

    // std::cout << "Transf matrix from result:\n";
    // mat_print_pretty(T_0_7, 3, 4);

    // double T_0_7_given[12];
    // aa_tf_qv2tfmat(qu_T.data, p_T_data, T_0_7_given);
    // std::cout << "Transf matrix given:\n";
    // mat_print_raw(T_0_7_given, 3, 4);

    // // struct aa_dmat T_result = AA_DMAT_INIT(3, 4, T_0_7, 3);
    // // struct aa_dmat T_given = AA_DMAT_INIT(3, 4, T_0_7_given, 3);
    // admeq( "result T == given T", T_0_7, T_0_7_given, AA_EPSILON, 12 );


    /* Check with forward kinematics, using modified (proximal) DH per M. Gong et al. */
    double qu_0_1[4];
    double v_0_1[3];
    aa_tf_dhprox2qv(0, 0, d1, q1, qu_0_1, v_0_1);
    double qu_1_2[4];
    double v_1_2[3];
    aa_tf_dhprox2qv(alpha1, a1, d2, q2, qu_1_2, v_1_2);
    double qu_2_3[4];
    double v_2_3[3];
    aa_tf_dhprox2qv(alpha2, a2, d3, q3, qu_2_3, v_2_3);
    double qu_3_4[4];
    double v_3_4[3];
    aa_tf_dhprox2qv(alpha3, a3, d4, q4, qu_3_4, v_3_4);
    double qu_4_5[4];
    double v_4_5[3];
    aa_tf_dhprox2qv(alpha4, a4, d5, q5, qu_4_5, v_4_5);
    double qu_5_6[4];
    double v_5_6[3];
    aa_tf_dhprox2qv(alpha5, a5, d6, q6, qu_5_6, v_5_6);
    double qu_6_7[4];
    double v_6_7[3];
    aa_tf_dhprox2qv(alpha6, a6, d7, q7, qu_6_7, v_6_7);

    double qu_0_2[4];
    double v_0_2[3];
    aa_tf_qv_chain(qu_0_1, v_0_1, qu_1_2, v_1_2, qu_0_2, v_0_2);
    double qu_0_3[4];
    double v_0_3[3];
    aa_tf_qv_chain(qu_0_2, v_0_2, qu_2_3, v_2_3, qu_0_3, v_0_3);
    double qu_0_4[4];
    double v_0_4[3];
    aa_tf_qv_chain(qu_0_3, v_0_3, qu_3_4, v_3_4, qu_0_4, v_0_4);
    double qu_0_5[4];
    double v_0_5[3];
    aa_tf_qv_chain(qu_0_4, v_0_4, qu_4_5, v_4_5, qu_0_5, v_0_5);
    double qu_0_6[4];
    double v_0_6[3];
    aa_tf_qv_chain(qu_0_5, v_0_5, qu_5_6, v_5_6, qu_0_6, v_0_6);
    double qu_0_7[4];
    double v_0_7[3];
    aa_tf_qv_chain(qu_0_6, v_0_6, qu_6_7, v_6_7, qu_0_7, v_0_7);

    double T_0_7[12];
    aa_tf_qv2tfmat(qu_0_7, v_0_7, T_0_7);
    std::cout << "Transf matrix from result:\n";
    mat_print_pretty(T_0_7, 3, 4);

    double T_0_7_given[12];
    aa_tf_qv2tfmat(qu_T.data, p_T_data, T_0_7_given);
    std::cout << "Transf matrix given:\n";
    mat_print_raw(T_0_7_given, 3, 4);

    // struct aa_dmat T_result = AA_DMAT_INIT(3, 4, T_0_7, 3);
    // struct aa_dmat T_given = AA_DMAT_INIT(3, 4, T_0_7_given, 3);
    admeq( "result T == given T", T_0_7, T_0_7_given, AA_EPSILON, 12 );


    std::cout << "Got to the end of ik7dof()!\n\n";
    return;
}


int main(int argc, char ** argv) 
{
    std::cout << "Got into main() successfully!\n\n";

    /* Following  "Analytical Inverse Kinematics and Self-Motion 
    Application for 7-DOF Redundant Manipulator" - M. Gong et al.*/

    /* See Figure 1 in M. Gong et al. for clarification on the DH parameters below*/
    double d1 = 0.381; // meter. These dimensions are for the Schunk arm
    double d3 = 0.33;
    double d5 = 0.32;
    double d7 = 0.29;

    /* User-assigned position of tool */
    double p_T_data[] = {d3*cos(M_PI_4) + d5 + d7*cos(M_PI_4), 
                         0, 
                         d1 + d3*cos(M_PI_4) - d7*cos(M_PI_4)};

    /* User-assigned orientation/rotation of tool */
    // double angle = M_PI / 2;
    double angle = M_PI*3/4;
    struct amino::YAngle y_angle{angle};
    struct amino::Quat qu_T{y_angle};

    /* Additional parameters*/
    double elbow_sign_param = 1; // either 1 or -1
    double elbow_ang_param = 0; // from -pi to pi
    double arm_sign_param = 1; // either 1 or -1
    double wrist_sign_param = 1; // either 1 or -1

    /* Call ik function */
    ik7dof(d1,
           d3,
           d5,
           d7, 
           p_T_data, 
           qu_T, 
           elbow_sign_param, 
           elbow_ang_param, 
           arm_sign_param, 
           wrist_sign_param);

    std::cout << "Got to the end of main()! Returning 0 next...\n\n";

    return 0;
}