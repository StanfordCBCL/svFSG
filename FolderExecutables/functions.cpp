#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#include "vessel.h"
#include "functions.h"
#include "matfun.h"

using std::string;
using std::vector;
using std::cout;

template <typename T>
//C++ template to print vessel contents
std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
    for (int i = 0; i < v.size(); ++i) { 
        os << v[i]; 
        if (i != v.size() - 1) 
            os << " "; 
    }
    return os; 
}

void update_time_step(vessel& curr_vessel) {
    //Solves equilibrium equations at the current time point and updates kinetic variables
    //Find current time step
    int n_alpha = curr_vessel.n_alpha;
    int nts = curr_vessel.nts;
    int sn = curr_vessel.sn;
    double s = sn * curr_vessel.dt;
    double rhoR_s0 = 0, rhoR_s1 = 0; //Total mass at current step
    double tol = 1E-14; //Convergence tolerance
    int iter = 0, equil_check = 0, unloaded_check = 0;
    double delta_sigma, delta_tauw, mb_equil;

    //Get initial guess for radius
    curr_vessel.a_mid[sn] = curr_vessel.a_mid[sn - 1];
    curr_vessel.a_act[sn] = curr_vessel.a_act[sn - 1];
    update_kinetics(curr_vessel);

    //Find new equilibrium geometry.
    //std::cout << curr_vessel.a_mid[sn]<< std::endl;
    equil_check = find_iv_geom(&curr_vessel);

    //Find real mass density production at the current time step iteratively
    double mass_check = 0.0;
    do {
        iter++;
        //Value from previous prediction
        rhoR_s0 = curr_vessel.rhoR[sn];

        //Update prediction
        update_kinetics(curr_vessel);
        equil_check = find_iv_geom(&curr_vessel);

        rhoR_s1 = curr_vessel.rhoR[sn];
        mass_check = abs((rhoR_s1 - rhoR_s0) / rhoR_s0);
    } while (mass_check > tol&& iter < 100);

    //--------------------------------------------------------------------
    ////For passive configuration calculation
    //curr_vessel.num_exp_flag = 1;
    
    //curr_vessel.T_act = 0.0;
    //curr_vessel.a_mid[sn] = curr_vessel.a_mid_pas[sn - 1];
    //equil_check = find_iv_geom(&curr_vessel);
    //curr_vessel.h_pas[sn] = curr_vessel.h[sn];
    //curr_vessel.a_pas[sn] = curr_vessel.a[sn];
    //curr_vessel.a_mid_pas[sn] = curr_vessel.a_mid[sn];

    ////Reset to active config
    //curr_vessel.T_act = curr_vessel.T_act_h;
    //curr_vessel.a_mid[sn] = curr_vessel.a_mid[sn - 1];
    //curr_vessel.a_act[sn] = curr_vessel.a_act[sn - 1];
    //equil_check = find_iv_geom(&curr_vessel);

    //curr_vessel.num_exp_flag = 0;
    //--------------------------------------------------------------------

    ////Find current mechanobiological perturbation


    delta_sigma = (curr_vessel.sigma_inv / curr_vessel.sigma_inv_h) - 1;
    delta_tauw = (curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1;

    mb_equil = 1 + curr_vessel.K_sigma_p_alpha_h[2] * delta_sigma -
                   curr_vessel.K_tauw_p_alpha_h[2] * delta_tauw;

    //Print current state
    // printf("%s %f %s %f %s %f %s %f\n", "Time:", s, "a: ", curr_vessel.a[sn], 

    //       "h:", curr_vessel.h[sn], "mb_equil:", mb_equil);
    //printf("%s %f\n", "Delta_sigma:", delta_sigma);
    //printf("%s %f\n", "sigma_inv:", curr_vessel.sigma_inv);
    //printf("%s %f\n", "sigma_inv_h:", curr_vessel.sigma_inv_h);

    //printf("%s %f\n", "Delta_sigma:", delta_sigma);
    //printf("%s %f\n", "delta_tauw:", curr_vessel.bar_tauw);
    //printf("%s %f\n", "delta_tauw_h:", curr_vessel.bar_tauw_h);
    //fflush(stdout);

}

void update_time_step_explicit(vessel& curr_vessel) {
    //Solves equilibrium equations at the current time point and updates kinetic variables
    //Find current time step
    int n_alpha = curr_vessel.n_alpha;
    int nts = curr_vessel.nts;
    int sn = curr_vessel.sn;
    double s = sn * curr_vessel.dt;
    double rhoR_s0 = 0, rhoR_s1 = 0; //Total mass at current step
    double tol = 1E-14; //Convergence tolerance
    int iter = 0, equil_check = 0, unloaded_check = 0;
    double delta_sigma, delta_tauw, mb_equil;



    //Get initial guess for radius
    if (sn > 0){
        curr_vessel.a_mid[sn] = curr_vessel.a_mid[sn - 1];
        curr_vessel.a_act[sn] = curr_vessel.a_act[sn - 1];
    }

    update_kinetics(curr_vessel);
    //Find new equilibrium geometry.
    equil_check = find_iv_geom(&curr_vessel);

    double mass_check = 0.0;
    do {
        iter++;
        //Value from previous prediction
        rhoR_s0 = curr_vessel.rhoR[sn];

        //Update prediction
        update_kinetics(curr_vessel);
        equil_check = find_iv_geom(&curr_vessel);

        rhoR_s1 = curr_vessel.rhoR[sn];
        mass_check = abs((rhoR_s1 - rhoR_s0) / rhoR_s0);
    } while (mass_check > tol && iter < 100);


    delta_sigma = (curr_vessel.sigma_inv / curr_vessel.sigma_inv_h) - 1;
    delta_tauw = (curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1;

    //std::cout << "delta_sigma " << delta_sigma << std::endl;
    //std::cout << "delta_tauw " << delta_tauw << std::endl;

    mb_equil = 1 + curr_vessel.K_sigma_p_alpha_h[2] * delta_sigma -
                   curr_vessel.K_tauw_p_alpha_h[2] * delta_tauw;

    //Print current state
    //printf("%s %f %s %f %s %f %s %f\n", "Time:", s, "a: ", curr_vessel.a[sn], 
    //        "h:", curr_vessel.h[sn], "mb_equil:", mb_equil);
    //printf("%s %f\n", "Delta_sigma:", delta_sigma);
    //printf("%s %f\n", "sigma_inv:", curr_vessel.sigma_inv);
    //printf("%s %f\n", "sigma_inv_h:", curr_vessel.sigma_inv_h);
    //printf("%s %f\n", "Delta_sigma:", delta_sigma);
    //printf("%s %f\n", "delta_tauw:", curr_vessel.bar_tauw);
    //printf("%s %f\n", "delta_tauw_h:", curr_vessel.bar_tauw_h);
    //fflush(stdout);


}


void update_time_step_handshake(vessel& curr_vessel, int iter_arg) {
    //Solves equilibrium equations at the current time point and updates kinetic variables
    //Find current time step
    int n_alpha = curr_vessel.n_alpha;
    int nts = curr_vessel.nts;
    int sn = curr_vessel.sn;
    double s = sn * curr_vessel.dt;
    //double rhoR_s0 = 0, rhoR_s1 = 0; //Total mass at current step
    double tol = 1E-14; //Convergence tolerance
    int iter = 0, equil_check = 0, unloaded_check = 0;
    double delta_sigma, delta_tauw, mb_equil;
    double rhoR_n;
    double rho_n;

    //First up, update F from FSI
    //vector<double> F_s(9,0);
    //for (int i = 0; i < 9; i++) {
    //     F_s[i] = curr_vessel.F_s[sn*9+i];
    //}
    rhoR_n = curr_vessel.rhoR[sn];
    rho_n = curr_vessel.rho[sn];
    
    if (!iter_arg) {
        sn = sn + 1;
        curr_vessel.sn = sn;
        s = sn * curr_vessel.dt;
        curr_vessel.s = s;
        curr_vessel.rho[sn] = curr_vessel.rho[sn-1];
        curr_vessel.rhoR[sn] = curr_vessel.rhoR[sn-1];
    } 

    //F_s = mat_mul(curr_vessel.F_curr,F_s);

    for (int i = 0; i < 9; i++) {
        curr_vessel.F_s[sn*9+i] = curr_vessel.F_curr[i];
    }

    update_kinetics(curr_vessel);
    update_sigma_handshake(&curr_vessel);

    curr_vessel.J_di = (curr_vessel.rhoR[sn] / curr_vessel.rho[sn]) /
    (rhoR_n / rho_n);

    ////Find current mechanobiological perturbation
    delta_sigma = (curr_vessel.sigma_inv / curr_vessel.sigma_inv_h) - 1;
    delta_tauw = (curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1;

    //Reset deformation gradient
    curr_vessel.F_curr = {1,0,0,0,1,0,0,0,1};

    //Print current state
    //printf("%s %f %s %f %s %f\n", "Time:", s, "delta_sigma: ", delta_sigma, 
    //       "delta_tauw:", delta_tauw);
    //fflush(stdout);

}

int ramp_pressure_test(void* curr_vessel, double P_low, double P_high) {
    int sn =((struct vessel*) curr_vessel)->sn;
    int equil_check = 0;
    int num_P = 100;
    double P_incr = (P_high - P_low) / num_P;


    ((struct vessel*) curr_vessel)->P = P_low;
    //std::cout << "Start ramp..." << std::endl;
    for (int i = 0; i < num_P; i++) {
        equil_check = find_iv_geom(curr_vessel);
        ((struct vessel*) curr_vessel)->P += P_incr;
    }
    return equil_check;
}

int run_pd_test(vessel& curr_vessel, double P_low, double P_high, double lambda_z_test) {

    int sn = curr_vessel.sn;

    double lambda_z_store = curr_vessel.lambda_z_curr;
    double P_store = curr_vessel.P;

    int equil_check = 0;
    int num_P = 100;
    double P_incr = (P_high - P_low) / num_P;

    curr_vessel.lambda_z_curr = lambda_z_test;
    curr_vessel.P = P_low;

    for (int i = 0; i < num_P; i++) {

        equil_check = find_iv_geom(&curr_vessel);
        curr_vessel.Exp_out << curr_vessel.P << "\t" << curr_vessel.a_mid[sn] << "\n";
        curr_vessel.P += P_incr;

    }

    curr_vessel.lambda_z_curr = lambda_z_store;
    curr_vessel.P = P_store;

    return 0;
}

int find_equil_geom(void* curr_vessel) {
    //Finds the mechanobiologically equilibrated geometry for a given set of loads inclduing
    //pressure, flow, and axial stretch with a set of G&R parameter values from the original
    //homeostatic state

    //Loading changes
    double gamma = ((struct vessel*) curr_vessel)->P 
        / ((struct vessel*) curr_vessel)->P_h; //Fold change in pressure from homeostatic
    double epsilon = ((struct vessel*) curr_vessel)->Q 
        / ((struct vessel*) curr_vessel)->Q_h; //Fold chance in flow from homeostatic
    double lambda = ((struct vessel*) curr_vessel)->lambda_z_curr 
        / ((struct vessel*) curr_vessel)->lambda_z_h; //Fold change in axial stretch from homeostatic

    //Homeostatic geometry
    double a_h = ((struct vessel*) curr_vessel)->a_h;
    double h_h = ((struct vessel*) curr_vessel)->h_h;

    //Initial guesses based on the loading changes
    double a_e_guess = pow(epsilon, 1.0 / 3.0) * a_h;
    double h_e_guess = gamma * pow(epsilon, 1.0 / 3.0) * h_h;
    double rho_c_e_guess = ((struct vessel*) curr_vessel)->rhoR_alpha_h[2] + ((struct vessel*) curr_vessel)->rhoR_alpha_h[3] + 
                           ((struct vessel*) curr_vessel)->rhoR_alpha_h[4] + ((struct vessel*) curr_vessel)->rhoR_alpha_h[5];
    double f_z_e_guess = ((struct vessel*) curr_vessel)->f_h* (h_e_guess * (2 * a_e_guess + h_e_guess)) / (h_h * (2 * a_h + h_h));

    //For psuedo-time dependent evolutions
    //if (0.98 * ((struct vessel*) curr_vessel)->a_h > ((struct vessel*) curr_vessel)->a_e || ((struct vessel*) curr_vessel)->a_e > 1.02 * ((struct vessel*) curr_vessel)->a_h){
        //a_e_guess = ((struct vessel*) curr_vessel)-> a_e;
        //h_e_guess = ((struct vessel*) curr_vessel)-> h_e;
        //rho_c_e_guess = ((struct vessel*) curr_vessel)-> rho_c_e;
        //f_z_e_guess = ((struct vessel*) curr_vessel)-> f_z_e;
    //}
    
    //printf("%s %f %s %f %s %f\n", "Time:", ((struct vessel*) curr_vessel)->s, "a_e_guess: ", ((struct vessel*) curr_vessel)->a_e, "h_e_guess:", ((struct vessel*) curr_vessel)->h_e);
    //fflush(stdout);

    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;

    int status;
    size_t iter = 0;

    const size_t n = 4;

    gsl_multiroot_function f = { &equil_obj_f, n, curr_vessel };
    double x_init[4] = {a_e_guess, h_e_guess, rho_c_e_guess, f_z_e_guess};
    gsl_vector* x = gsl_vector_alloc(n);

    for (int i = 0; i < n; i++) {
        gsl_vector_set(x, i, x_init[i]);
    }

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, n);

    gsl_multiroot_fsolver_set(s, &f, x);

    //print_state_mr(iter, s);

    double epsabs = 1e-7;
    double epsrel = 1e-3;
    int check_tol = 0;
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        //print_state_mr(iter, s);

        if (status)
            break;

        //status = gsl_multiroot_test_delta(s->dx, s->x, epsabs, epsrel);
        status = gsl_multiroot_test_residual(s->f, epsabs);
        //for (int i = 0; i < n; i++)
            //check_tol += (gsl_vector_get(s->f, i) > epsabs);

        //if (check_tol >= 0)

    } while (status == GSL_CONTINUE && iter < 1000);
    
    //printf("status = %s \n", gsl_strerror(status));

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;

}

int equil_obj_f(const gsl_vector* x, void* curr_vessel, gsl_vector* f) {
    //Mechanobiologically equilibrated objective function
    //Unknown input variables
    const double a_e_guess = gsl_vector_get(x, 0);
    const double h_e_guess = gsl_vector_get(x, 1);
    const double rho_c_e_guess = gsl_vector_get(x, 2);
    const double f_z_e_guess = gsl_vector_get(x, 3);

    //Equations for J1
    //Stress values from equilibrium equations
    double sigma_e_th_lmb = ((struct vessel*) curr_vessel)->P * a_e_guess / h_e_guess;
    double sigma_e_z_lmb = f_z_e_guess / (M_PI * h_e_guess * (2 * a_e_guess + h_e_guess));

    //WSS from Pousielle flow, constant viscosity
    int sn = ((struct vessel*) curr_vessel)->sn;
    double a_store = ((struct vessel*) curr_vessel)->a[sn];
    ((struct vessel*) curr_vessel)->a[sn] = a_e_guess;
	double mu = get_app_visc(curr_vessel, sn);
    ((struct vessel*) curr_vessel)->a[sn] = a_store;
    double bar_tauw_e = 4*mu*((struct vessel*) curr_vessel)->Q/(3.14159265*pow(a_e_guess*100, 3)); //((struct vessel*) curr_vessel)->Q / pow(a_e_guess, 3);

    //Ratio of stress:WSS mediated matrix production
    double eta_K = ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[1] /
        ((struct vessel*) curr_vessel)->K_tauw_p_alpha_h[1];

    //Stress and WSS deviations from HS state
    double delta_sigma = (sigma_e_th_lmb + sigma_e_z_lmb) /
        (((struct vessel*) curr_vessel)->sigma_h[1] + ((struct vessel*) curr_vessel)->sigma_h[2]) - 1;
    double delta_tauw = bar_tauw_e / ((struct vessel*) curr_vessel)->bar_tauw_h - 1;

    //Equations for J2
    //Homeostatic geometry
    double a_h = ((struct vessel*) curr_vessel)->a_h;
    double h_h = ((struct vessel*) curr_vessel)->h_h;

    //Equilibrated stretches
    double lambda_r_e = h_e_guess / h_h;
    double lambda_th_e = (a_e_guess + h_e_guess / 2) / (a_h + h_h / 2);
    double lambda_z_e = ((struct vessel*) curr_vessel)->lambda_z_curr;
    double F_e[3] = { lambda_r_e, lambda_th_e, lambda_z_e };

    //Equilibrated volume change
    double J_e = lambda_r_e * lambda_th_e * lambda_z_e; 

    //Equilibrated mass densities
    //Equilibrated elastin density, account for possible elastin degradation
    double rho_el_e = ((struct vessel*) curr_vessel)->rhoR_alpha[0 * sn + sn] / J_e; 

    //Ratio of degrdation rate for smc:col
    double eta_q = ((struct vessel*) curr_vessel)->k_alpha_h[1] / 
        ((struct vessel*) curr_vessel)->k_alpha_h[2];

    //Ratio of stress med prod for smc:col must be same for WSS
    double eta_ups = ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[1] /
        ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[2];

    //Equilibrated muscle density
    double rhoR_c_curr = 0;
    double rhoR_c_h_total = 0;
    double phi_hat_c_curr = 0;
    double rho_c_k_guess[4] = {0, 0, 0, 0};
    for (int k = 2; k < 6; k++){
        rhoR_c_curr = ((struct vessel*) curr_vessel)->rhoR_alpha_h[k];
        rhoR_c_h_total += rhoR_c_curr; 
    }
    for (int k = 2; k < 6; k++){
        rhoR_c_curr = ((struct vessel*) curr_vessel)->rhoR_alpha_h[k];
        phi_hat_c_curr = rhoR_c_curr / rhoR_c_h_total;
        rho_c_k_guess[k - 2] = phi_hat_c_curr * rho_c_e_guess; 
    }
 
    double rho_m_e = ((struct vessel*) curr_vessel)->rhoR_alpha_h[1] / J_e *
        pow(J_e * rho_c_e_guess / rhoR_c_h_total, eta_q * eta_ups);

    //Array of equilibrated densities
    //Add individual collagen densities
    double rho_alpha[6] = { rho_el_e , rho_m_e , rho_c_k_guess[0], rho_c_k_guess[1], rho_c_k_guess[2], rho_c_k_guess[3] };
    double rho_h = ((struct vessel*) curr_vessel)->rhoR_h;

    //Equations for J3 & J4
    //number of constituents
    int n_alpha = ((struct vessel*) curr_vessel)->n_alpha;

    //cauchy stress hat for each const. in each direction
    vector<double> hat_sigma_alpha_dir(3 * n_alpha, 0);
    vector<double> sigma_e_dir(3, 0);

    //equilibrated stresses
    vector<double> sigma_e(3, 0);

    //equilibrated active stress
    double C = ((struct vessel*) curr_vessel)->CB -
        ((struct vessel*) curr_vessel)->CS * delta_tauw;
    double lambda_act = 1.0;
    double parab_act = 1 - pow((((struct vessel*) curr_vessel)->lambda_m - lambda_act) /
        (((struct vessel*) curr_vessel)->lambda_m - ((struct vessel*) curr_vessel)->lambda_0), 2);
    double hat_sigma_act_e = ((struct vessel*) curr_vessel)->T_act * (1 - exp(-pow(C, 2))) * lambda_act * parab_act;

    //equilibrated constituent strech
    double lambda_alpha_ntau_s = 0.0;

    //equilbirated 2nd PK stress hat
    double hat_S_alpha = 0.0;

    for (int alpha = 0; alpha < n_alpha; alpha++) {
        for (int dir = 0; dir < 3; dir++) {

            //Check if vessel is anisotropic
            if (((struct vessel*) curr_vessel)->eta_alpha_h[alpha] >= 0) {
                //Constituent stretch is equal to deposition stretch
                lambda_alpha_ntau_s = ((struct vessel*) curr_vessel)->g_alpha_h[alpha];

                //2nd PK stress hat at equilibrium
                hat_S_alpha = ((struct vessel*) curr_vessel)->c_alpha_h[2 * alpha] * (pow(lambda_alpha_ntau_s, 2) - 1) *
                    exp(((struct vessel*) curr_vessel)->c_alpha_h[2 * alpha + 1] * pow(pow(lambda_alpha_ntau_s, 2) - 1, 2));

                //Cauchy stress hat at equilibrium
                hat_sigma_alpha_dir[alpha * dir + dir] = ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir] *
                hat_S_alpha * ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir];

            }
            else {
                //2nd PK stress hat at equilibrium
                hat_S_alpha = ((struct vessel*) curr_vessel)->c_alpha_h[2 * alpha];

                //Cauchy stress hat at equilibrium
                hat_sigma_alpha_dir[alpha * dir + dir] = ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir] * hat_S_alpha *
                    ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + dir];

                //Check if the consituent is present from the initial time point
                //Account for volume change and mixture deformation
                if (((struct vessel*) curr_vessel)->k_alpha_h[alpha] == 0) {
                    hat_sigma_alpha_dir[alpha * dir + dir] = F_e[dir] * hat_sigma_alpha_dir[alpha * dir + dir] * F_e[dir];
                }

            }

            sigma_e_dir[dir] += rho_alpha[alpha] / rho_h * hat_sigma_alpha_dir[alpha * dir + dir];

            if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1 && dir == 1) {
                sigma_e_dir[dir] += rho_alpha[alpha] / rho_h * hat_sigma_act_e;
            }
        }
    }

    //Four objective equations
    double J1 = eta_K*delta_sigma - delta_tauw; //mechano-mediated matrix production equation
    double J2 = rho_el_e + rho_m_e + rho_c_e_guess - rho_h; //mixture mass balance equation
    double J3 = sigma_e_dir[1] - sigma_e_dir[0] - sigma_e_th_lmb; //circumferential const - lmb
    double J4 = sigma_e_dir[2] - sigma_e_dir[0] - sigma_e_z_lmb; //axial const - lmb

    gsl_vector_set(f, 0, J1);
    gsl_vector_set(f, 1, J2);
    gsl_vector_set(f, 2, J3);
    gsl_vector_set(f, 3, J4);

    //Store equilibrated results  
    //if (a_e_guess == a_e_guess) {
        ((struct vessel*) curr_vessel)->a_e = a_e_guess;
        ((struct vessel*) curr_vessel)->h_e = h_e_guess;
        ((struct vessel*) curr_vessel)->rho_c_e = rho_c_e_guess * J_e;
        ((struct vessel*) curr_vessel)->rho_m_e = rho_m_e * J_e;
        ((struct vessel*) curr_vessel)->f_z_e = f_z_e_guess;
        ((struct vessel*) curr_vessel)->mb_equil_e = 1 + ((struct vessel*) curr_vessel)->K_sigma_p_alpha_h[2] * delta_sigma -
                                                     ((struct vessel*) curr_vessel)->K_tauw_p_alpha_h[2] * delta_tauw;
    //}     

    return GSL_SUCCESS;
}

int print_state_mr(size_t iter, gsl_multiroot_fsolver* s)
{

    printf("iter = %3u", iter);
    printf("x = ");
    for (int i = 0; i < s->x->size; i++) {
        printf("%.3f ", gsl_vector_get(s->x, i));
    }
    printf("f(x) = ");
    for (int i = 0; i < s->x->size; i++) {
        printf("%.3f ", gsl_vector_get(s->f, i));
    }
    printf("\n");

    return 0;
}

int find_tf_geom(void* curr_vessel) {
    //For geometries defined in the loaded configuration this code solves for the unloaded geometric 
    //variables at the current time point and stores them in the current v
    //vessel structure. Unloaded is also
    //referred to as traction free (but not stress free).

    //Local vars to store current pressure, force, and stretch
    double P_temp = ((struct vessel*)curr_vessel)->P;
    double f_temp = ((struct vessel*)curr_vessel)->f;
    double lambda_z_temp = ((struct vessel*)curr_vessel)->lambda_z_curr;

    //Get initial guesses from the loaded geometry
    double lambda_th_ul = 0.95;
    double lambda_z_ul = 0.95 * ((struct vessel*) curr_vessel)->lambda_z_curr;

    //Update vesseel loads to zero for traction free
    ((struct vessel*) curr_vessel)->P = 0.0;
    ((struct vessel*) curr_vessel)->f = 0.0;
    //((struct vessel*) curr_vessel)->T_act = 0.0;

    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;

    int status;
    size_t iter = 0;

    const size_t n = 2;

    gsl_multiroot_function f = { &tf_obj_f, n, curr_vessel };
    double x_init[2] = { lambda_th_ul, lambda_z_ul };
    gsl_vector* x = gsl_vector_alloc(n);

    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, 2);

    gsl_multiroot_fsolver_set(s, &f, x);

    //print_state_mr(iter, s);

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        //print_state_mr(iter, s);

        if (status)   // check if solver is stuck
            break;

        status = gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 100);

    //printf("status = %s\n", gsl_strerror(status));

    //Store the traction free results
    int sn = ((struct vessel*) curr_vessel)->sn;
    ((struct vessel*) curr_vessel)->A_mid[sn] = gsl_vector_get(s->x, 0) * ((struct vessel*) curr_vessel)->a_mid[sn];
    ((struct vessel*) curr_vessel)->lambda_z_pre[sn] = 1 / gsl_vector_get(s->x, 1);
    ((struct vessel*) curr_vessel)->H[sn] = 1.0 / (gsl_vector_get(s->x, 0) * gsl_vector_get(s->x, 1)) * ((struct vessel*) curr_vessel)->h[sn];
    ((struct vessel*)curr_vessel)->A[sn] = ((struct vessel*)curr_vessel)->A_mid[sn] - ((struct vessel*)curr_vessel)->H[sn] / 2;

    //Return current vars to the loaded conditions
    ((struct vessel*) curr_vessel)->P = P_temp;
    ((struct vessel*) curr_vessel)->f = f_temp;
    //((struct vessel*) curr_vessel)->T_act = ((struct vessel*) curr_vessel)->T_act_h;
    ((struct vessel*) curr_vessel)->lambda_z_curr = lambda_z_temp;

    status = find_iv_geom(curr_vessel);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;

}

int tf_obj_f(const gsl_vector* x, void* curr_vessel, gsl_vector* f) {

    //Seperate out inputs
    const double lambda_th_ul_guess = gsl_vector_get(x, 0);
    const double lambda_z_ul_guess = gsl_vector_get(x, 1);

    //Finds the difference in the theoretical stress from Laplace for deformed mixture
    //from the stress calculated from the mixture equations
    int sn = ((struct vessel*) curr_vessel)->sn;

    //Reference config.
    double a_mid_0 = ((struct vessel*) curr_vessel)->a_mid[0];
    double lambda_z_0 = ((struct vessel*)curr_vessel)->lambda_z_tau[0];

    //Current loaded config
    double a_mid = ((struct vessel*) curr_vessel)->a_mid[sn];

    //Current stretches ref -> loaded
    double lambda_th_ref = a_mid / a_mid_0;
    double lambda_z_ref = lambda_z_0;

    //Update current total stretches
    ((struct vessel*) curr_vessel)->lambda_th_curr = lambda_th_ul_guess * lambda_th_ref;
    ((struct vessel*) curr_vessel)->lambda_z_curr = lambda_z_ul_guess * lambda_z_ref;

    update_sigma(curr_vessel);

    //Should be 0 for the traction-free configuration
    double J1 = ((struct vessel*) curr_vessel)->sigma[1];
    double J2 = ((struct vessel*) curr_vessel)->sigma[2];

    gsl_vector_set(f, 0, J1);
    gsl_vector_set(f, 1, J2);

    return GSL_SUCCESS;

}

int find_iv_geom(void* curr_vessel) {
    //Finds the loaded conifiguration for a given pressure and axial stretch

    int status;
    int iter = 0;
    int max_iter = 10000;

    int sn = ((struct vessel*) curr_vessel)->sn;

    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
    gsl_function f = { &iv_obj_f, curr_vessel };

    //Set search range for new mid radius
    double a_mid_act, a_mid_high, a_mid_low;
    a_mid_low = 0.90 * ((struct vessel*) curr_vessel)->a_mid[sn];
    a_mid_high = 1.1 * ((struct vessel*) curr_vessel)->a_mid[sn];

    gsl_root_fsolver_set(s, &f, a_mid_low, a_mid_high);

    //printf("Using %s method \n", gsl_root_fsolver_name(s));
    //printf("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err (est)");

    a_mid_act = 0.0;
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        a_mid_act = gsl_root_fsolver_root(s);
        a_mid_low = gsl_root_fsolver_x_lower(s);
        a_mid_high = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(a_mid_low, a_mid_high, 0, pow(10, -12));

        if (status == GSL_SUCCESS) {
            //printf("Loaded Config Converged:\n");
            //printf("%5d [%.7f, %.7f] %.7f %.7fs\n", iter, a_mid_low, a_mid_high, a_mid_act, a_mid_high - a_mid_low);
        }

        //printf("%5d [%.7f, %.7f] %.7f %.7fs\n", iter, a_mid_low, a_mid_high, a_mid_act, a_mid_high - a_mid_low);
    } while (status == GSL_CONTINUE && iter < max_iter);

    //if (iter > max_iter-2){
    //    printf("status = NOT CONVERGED\n");
    //}
    //printf("status = %s\n", gsl_strerror(status));

    //int set_iv = iv_obj_f(a_mid_act, curr_vessel);
    //((struct vessel*) curr_vessel)->a_mid[sn] = a_mid_act;

    return status;
}

double iv_obj_f(double a_mid_guess, void* curr_vessel) {
    //Finds the difference in the theoretical stress from Laplace for deformed mixture
    //from the stress calculated from the mixture equations

    int sn = ((struct vessel*) curr_vessel)->sn;
    double a = 0.0, h = 0.0, lambda_t = 0.0, lambda_z = 0.0, J_s = 0.0;
    double mu = 0.0;

    if (sn > 0 || ((struct vessel*) curr_vessel)->num_exp_flag == 1) {
        lambda_t = a_mid_guess / ((struct vessel*) curr_vessel)->a_mid_h;
        lambda_z = ((struct vessel*) curr_vessel)->lambda_z_curr;
        J_s = ((struct vessel*) curr_vessel)->rhoR[sn] / ((struct vessel*) curr_vessel)->rho[sn];

        //Update vessel geometry for calculation of next time step
        h = J_s / (lambda_t * lambda_z) * ((struct vessel*) curr_vessel)->h_h;
        a = a_mid_guess - h / 2;
        ((struct vessel*) curr_vessel)->a_mid[sn] = a_mid_guess;
        ((struct vessel*) curr_vessel)->a[sn] = a;
        ((struct vessel*) curr_vessel)->h[sn] = h;
        ((struct vessel*) curr_vessel)->lambda_th_curr = lambda_t;
        ((struct vessel*) curr_vessel)->lambda_z_curr = lambda_z;
        //Update WSS from Q Flow
        if (((struct vessel*) curr_vessel)->wss_calc_flag > 0) {
            mu = get_app_visc(curr_vessel, sn);
            ((struct vessel*) curr_vessel)->bar_tauw = 4 * mu * ((struct vessel*)curr_vessel)->Q / (3.14159265 * pow(a * 100, 3));
        }

    }
    else {
        //This only checks that the current state is the initial equilibrium state
        lambda_t = a_mid_guess / ((struct vessel*) curr_vessel)->a_mid_h;
        lambda_z = ((struct vessel*) curr_vessel)->lambda_z_curr;
        J_s = 1.0;

        //Update vessel geometry for calculation of next time step
        h = J_s / (lambda_t * lambda_z) * ((struct vessel*) curr_vessel)->h_h;
        a = a_mid_guess - h / 2;
        ((struct vessel*) curr_vessel)->a_mid[sn] = a_mid_guess;
        ((struct vessel*) curr_vessel)->a_act[sn] = a;
        ((struct vessel*) curr_vessel)->a[sn] = a;
        ((struct vessel*) curr_vessel)->h[sn] = h;
        ((struct vessel*) curr_vessel)->lambda_th_curr = 1.0;
        ((struct vessel*) curr_vessel)->lambda_z_curr = lambda_z;
        //Update WSS from Q Flow (ELS)
        //((struct vessel*)curr_vessel)->bar_tauw = 4 * 0.04 * ((struct vessel*)curr_vessel)->Q / (3.14159265 * pow(a * 100, 3));
    }

    //Calculating sigma_t_th from pressure P
    double sigma_t_th = ((struct vessel*) curr_vessel)->P * a / h;

    update_sigma(curr_vessel);
    double sigma_t_calc = ((struct vessel*) curr_vessel)->sigma[1];

    ((struct vessel*) curr_vessel)->f = M_PI * h * (2 * a + h) *
        ((struct vessel*) curr_vessel)->sigma[2];

    double J = sigma_t_calc - sigma_t_th;

    return J;
}

void update_kinetics(vessel& curr_vessel) {
//    std cout << "Begin update kinetics ----------------------------" << std::endl;

    //This function updates the kinetics for G&R.
    int n_alpha = curr_vessel.n_alpha;
    int nts = curr_vessel.nts;
    double dt = curr_vessel.dt;
    double s = curr_vessel.s;
    int sn = curr_vessel.sn;
    int taun_min = 0;

    double tau_max = 100 * (1 / curr_vessel.k_alpha_h[3]); //max time of 10 half-lives

    //Differences in current mechanical state from the reference state
    //Stress invariant
    double delta_sigma = ((1 - curr_vessel.K_delta_sigma) * curr_vessel.sigma_inv / curr_vessel.sigma_inv_h) - 1;
    //Wall shear stress
    double delta_tauw = ((1 - curr_vessel.K_delta_tauw) * curr_vessel.bar_tauw / curr_vessel.bar_tauw_h) - 1;

    //Initialize pars for looping later
    double K_sigma_p = 0, K_tauw_p = 0, K_sigma_d = 0, K_tauw_d = 0;
    double upsilon_p = 0, upsilon_d = 0;

    double k_alpha_s = 0;
    double mR_alpha_s = 0;
    double rhoR_alpha_calc = 0;

    double mq_0 = 0, mq_1 = 0, mq_2;
    double q_0 = 0, q_1 = 0, q_2;
    double k_0 = 0, k_1 = 0, k_2;
    double rhoR_s = 0, rhoR_alpha_s = 0;
    double J_s = 0;

    int n = 0; //number of points in integration interval

    bool deg_check = 0;

    //Check if we've exceeded the initial history
    if (s > tau_max) {
        taun_min = sn - int(tau_max / dt);
    }
    else {
        taun_min = 0;
    }

    n = (sn - taun_min) + 1; //find number of integration pts

    //Loop through each constituent to update its mass density
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Find if constituent degrades
        deg_check = curr_vessel.k_alpha_h[alpha] > 0;

        if ((sn > 0 && deg_check) && curr_vessel.pol_only_flag == 0) {

            if (curr_vessel.alpha_infl[alpha] == 0) {

                //Get the gains for the current constituent
                K_sigma_p = curr_vessel.K_sigma_p_alpha_h[alpha];
                K_tauw_p = curr_vessel.K_tauw_p_alpha_h[alpha];

                K_sigma_d = curr_vessel.K_sigma_d_alpha_h[alpha];
                K_tauw_d = curr_vessel.K_tauw_d_alpha_h[alpha];

                //Update the stimulus functions for each constituent
                upsilon_p = 1 + K_sigma_p * delta_sigma - K_tauw_p * delta_tauw;
                if (upsilon_p < 0.1) {
                    upsilon_p = 0.1;
                }

                upsilon_d = 1 + K_sigma_d * pow(delta_sigma, 2) + K_tauw_d * pow(delta_tauw, 2);

                if (curr_vessel.rhoR_alpha[nts * alpha + sn] <= curr_vessel.rhoR_alpha_h[alpha]) {
                    rhoR_alpha_calc = curr_vessel.rhoR_alpha_h[alpha];
                }
                else {
                    rhoR_alpha_calc = curr_vessel.rhoR_alpha[nts * alpha + sn];
                }
            }
            else {

                upsilon_p = curr_vessel.ups_infl_p[sn];
                upsilon_d = 1 + curr_vessel.ups_infl_d[sn];

                rhoR_alpha_calc = curr_vessel.rhoR_alpha_h[alpha];

            }

            //Make sure productions don't become negative
            upsilon_p = (upsilon_p > 0)* upsilon_p;
            upsilon_d = (upsilon_d > 0)* upsilon_d;

            //Reset these to zero for each constituent
            k_alpha_s = 0;
            mR_alpha_s = 0;
            rhoR_alpha_s = 0;

            //Update kinetic values for the current time
            k_alpha_s = curr_vessel.k_alpha_h[alpha] * upsilon_d;
            
            mR_alpha_s = curr_vessel.k_alpha_h[alpha] * rhoR_alpha_calc * upsilon_p; //curr_vessel.rhoR_alpha_h[alpha];
            curr_vessel.k_alpha[nts * alpha + sn] = k_alpha_s;
            curr_vessel.mR_alpha[nts * alpha + sn] = mR_alpha_s;

            k_2 = curr_vessel.k_alpha[nts * alpha + sn];
            q_2 = 1.0;
            mq_2 = curr_vessel.mR_alpha[nts * alpha + sn] * q_2;

            //loop through and update constituent densities from previous time points
            //starting from the current time point and counting down is more efficient
            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) {

                //Trapazoidal rule     
                k_1 = curr_vessel.k_alpha[nts * alpha + taun];
                q_1 = exp(-(k_2 + k_1) * dt / 2) * q_2;
                mq_1 = curr_vessel.mR_alpha[nts * alpha + taun] * q_1;

                rhoR_alpha_s += (mq_2 + mq_1) * dt / 2;

                k_2 = k_1;
                q_2 = q_1;
                mq_2 = mq_1;

            }

            //Account for the cohort of material present initially
            if (taun_min == 0) {
                rhoR_alpha_s += curr_vessel.rhoR_alpha[nts * alpha + 0] * q_1;
            }

            //Update referential volume fraction
            curr_vessel.epsilonR_alpha[curr_vessel.nts * alpha + sn] = rhoR_alpha_s / curr_vessel.rho_hat_alpha_h[alpha];
            J_s += curr_vessel.epsilonR_alpha[curr_vessel.nts * alpha + sn];
        }
        else {
            //Precalculate maintenance or loss of cosntituents not produced
            rhoR_alpha_s = curr_vessel.rhoR_alpha[nts * alpha + sn];
            J_s += curr_vessel.epsilonR_alpha[curr_vessel.nts * alpha + sn];
        }

        curr_vessel.rhoR_alpha[nts * alpha + sn] = rhoR_alpha_s;
        rhoR_s += rhoR_alpha_s;

    }
    //Update spatial volume fractions
    for (int alpha = 0; alpha < n_alpha; alpha++) {
        curr_vessel.epsilon_alpha[curr_vessel.nts * alpha + sn] = curr_vessel.epsilonR_alpha[curr_vessel.nts * alpha + sn] / J_s;
    }

    curr_vessel.rhoR[sn] = rhoR_s;
    curr_vessel.rho[sn] = rhoR_s / J_s;

}


void update_sigma(void* curr_vessel) {

    //Get current time index
    double s = ((struct vessel*)curr_vessel)->s;
    int sn = ((struct vessel*)curr_vessel)->sn;
    int nts = ((struct vessel*)curr_vessel)->nts;
    double dt = ((struct vessel*)curr_vessel)->dt;
    int taun_min = 0;

    double tau_max = 10000 * (1 / ((struct vessel*)curr_vessel)->k_alpha_h[3]); //max time of 10 half-lives

    //Specify vessel geometry
    double a0 = ((struct vessel*)curr_vessel)->a[0];
    double h0 = ((struct vessel*)curr_vessel)->h[0];

    //Calculate vessel stretches
    double lambda_th_s = ((struct vessel*)curr_vessel)->lambda_th_curr;
    double lambda_z_s = ((struct vessel*)curr_vessel)->lambda_z_curr;

    //Calculate constituent specific stretches for evolving constituents at the current time
    int n_alpha = ((struct vessel*)curr_vessel)->n_alpha;
    vector<double> lambda_alpha_s(n_alpha, 0);
    double eta_alpha = 0;
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Check to see if constituent is isotropic
        eta_alpha = ((struct vessel*)curr_vessel)->eta_alpha_h[alpha];
        if (eta_alpha >= 0) {

            //Stretch is equal to the sqrt of I4
            lambda_alpha_s[alpha] = sqrt(pow(lambda_z_s * cos(eta_alpha), 2)
                + pow(lambda_th_s * sin(eta_alpha), 2));

            //Update stored current stretch if not numerical experiment
            if (((struct vessel*)curr_vessel)->num_exp_flag == 0) {
                ((struct vessel*)curr_vessel)->lambda_alpha_tau[nts * alpha + sn] = lambda_alpha_s[alpha];
            }
        }
    }

    //Find the current deformation gradient
    double J_s = ((struct vessel*)curr_vessel)->rhoR[sn] / ((struct vessel*)curr_vessel)->rho[sn];
    double F_s[3] = { J_s / (lambda_th_s * lambda_z_s), lambda_th_s, lambda_z_s };

    //Find the mechanical contributions of each constituent for each direction
    double a, h;
    double lambda_th_tau = 0;
    double lambda_z_tau = 0; //Assume constant axial stretch
    double J_tau = 1;
    double F_tau[3] = { 1, 1, 1 };
    double lambda_alpha_ntau_s = 0;
    double Q1 = 0, Q2 = 0;
    double F_alpha_ntau_s = 0;
    double hat_S_alpha = 0;
    double sigma[3] = { 0 };
    double lagrange = 0;
    double pol_mod = 0;

    //Local active variables
    double C = 0;
    double lambda_act = 0;
    double parab_act = 0;
    double hat_sigma_act = 0, sigma_act = 0;

    //Stiffness variables
    double hat_dSdC_alpha = 0;
    double hat_dSdC_act = 0;
    double Cbar[3] = { 0 };
    double Cbar_act = 0;
    vector<double> constitutive_return = { 0, 0 };

    //Integration variables
    //For mass
    double mq_0 = 0, mq_1 = 0, mq_2 = 0;
    double q_0 = 1.0, q_1 = 1.0, q_2 = 1.0;
    double k_0 = 0, k_1 = 0, k_2 = 0;

    int n = 0; //number of pts in integration interval

    //For stress
    vector<double> hat_sigma_0 = { 0, 0, 0 }, hat_sigma_1 = { 0, 0, 0 }, hat_sigma_2 = { 0, 0, 0 };
    //For active stress
    double a_act = 0;
    double k_act = ((struct vessel*) curr_vessel)->k_act;
    double q_act_0 = 0, q_act_1 = 0, q_act_2 = 0;
    double a_0 = 0, a_1 = 0, a_2 = 0;

    //For stiffness
    vector<double> hat_Cbar_0 = { 0, 0, 0 }, hat_Cbar_1 = { 0, 0, 0 }, hat_Cbar_2 = { 0, 0, 0 };

    //Boolean for checks
    bool deg_check = 0;

    //Determine if beyond initial time history
    if (s > tau_max) {
        taun_min = sn - int(tau_max / dt);
    }
    else {
        taun_min = 0;
    }

    n = (sn - taun_min) + 1;; //number of integration pts

    //Similar integration to that used for kinematics
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Trapz rule allows for fast heredity integral evaluation
        k_2 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + sn];
        q_2 = 1.0;
        mq_2 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + sn];

        //Find active radius from current cohort
        if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
            a_2 = ((struct vessel*) curr_vessel)->a[sn];
            q_act_2 = 1.0;
        }

        //Kinematics
        F_tau[0] = F_s[0], F_tau[1] = F_s[1], F_tau[2] = F_s[2];
        J_tau = J_s;

        //Find stress from current cohort
        for (int dir = 0; dir < 3; dir++) {

            constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, sn);
            hat_S_alpha = constitutive_return[0];
            hat_dSdC_alpha = constitutive_return[1];
            F_alpha_ntau_s = F_s[dir] / F_tau[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
            hat_sigma_2[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
            hat_Cbar_2[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;
        }

        //Boolean for whether the constituent increases ref mass density
        deg_check = ((struct vessel*)curr_vessel)->mR_alpha_h[alpha] > 0;

        //Check if during G&R or at initial time point
        if (sn > 0 && deg_check) {

            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) {

                //Find the 1st intermediate deformation gradient
                a = ((struct vessel*)curr_vessel)->a[taun];
                h = ((struct vessel*)curr_vessel)->h[taun];
                lambda_th_tau = (a + h / 2) / (a0 + h0 / 2);
                lambda_z_tau = ((struct vessel*)curr_vessel)->lambda_z_tau[taun];
                J_tau = ((struct vessel*)curr_vessel)->rhoR[taun] / ((struct vessel*)curr_vessel)->rho[taun];
                F_tau[0] = J_tau / (lambda_th_tau * lambda_z_tau);
                F_tau[1] = lambda_th_tau;
                F_tau[2] = lambda_z_tau;

                //Find 1st intermediate kinetics
                k_1 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + taun];
                q_1 = exp(-(k_2 + k_1) * dt / 2) * q_2;
                //mRq -> referential
                mq_1 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + taun] * q_1;

                //Find intermediate active state
                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    a_1 = a;
                    q_act_1 = exp(-k_act * dt) * q_act_2;
                }

                //Find 1st intermeidate stress component in each direction
                for (int dir = 0; dir < 3; dir++) {

                    constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, taun);
                    hat_S_alpha = constitutive_return[0];
                    hat_dSdC_alpha = constitutive_return[1];
                    F_alpha_ntau_s = F_s[dir] / F_tau[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
                    hat_sigma_1[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                    hat_Cbar_1[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;
                }

                //Find component in each direction
                for (int dir = 0; dir < 3; dir++) {

                    constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, taun - 1);
                    hat_S_alpha = constitutive_return[0];
                    hat_dSdC_alpha = constitutive_return[1];
                    F_alpha_ntau_s = F_s[dir] / F_tau[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
                    hat_sigma_0[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                    hat_Cbar_0[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

                    //Add to the stress and stiffness contribution in the given direction
                    sigma[dir] += (mq_2 * hat_sigma_2[dir] + mq_1 * hat_sigma_1[dir])
                        / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * dt / 2;
                    Cbar[dir] += (mq_2 * hat_Cbar_2[dir] + mq_1 * hat_Cbar_1[dir])
                        / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * dt / 2;
                }

                //Store active vars for next iteration
                //Find intermediate active state
                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    a_act += k_act * (q_act_2 * a_2 + q_act_1 * a_1) * dt / 2;
                    a_2 = a_1;
                    q_act_2 = q_act_1;
                }

                //Store intermediate kinetics for next iteration
                k_2 = k_1;
                q_2 = q_1;
                mq_2 = mq_1;

                //Store intermediate stress and stiffness for next iteration
                hat_sigma_2 = hat_sigma_1;
                hat_Cbar_2 = hat_Cbar_1;

            }

            //Add in the stress and stiffness contributions of the initial material
            if (taun_min == 0) {
                for (int dir = 0; dir < 3; dir++) {
                    sigma[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + 0]
                        / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * q_1 * hat_sigma_1[dir];
                    Cbar[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + 0]
                        / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * q_1 * hat_Cbar_1[dir];
                }
            }

        }
        //Initial time point and constituents with prescribed degradation profiles
        else {
            //Find stress from initial cohort          
            for (int dir = 0; dir < 3; dir++) {

                constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, 0);
                hat_S_alpha = constitutive_return[0];
                hat_dSdC_alpha = constitutive_return[1];
                F_alpha_ntau_s = F_s[dir] * ((struct vessel*)curr_vessel)->G_alpha_h[3 * alpha + dir];
                hat_sigma_2[dir] = F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s / J_s;
                hat_Cbar_2[dir] = F_alpha_ntau_s * F_alpha_ntau_s * hat_dSdC_alpha * F_alpha_ntau_s * F_alpha_ntau_s / J_s;

                sigma[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + sn] /
                    ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * hat_sigma_2[dir];
                Cbar[dir] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + sn] /
                    ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * hat_Cbar_2[dir];
            }

        }

        if (taun_min == 0 && ((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
            a_act += ((struct vessel*) curr_vessel)->a_act[0] * q_act_1;
        }


    }

    //Find active stress contribtion
    //add in initial active stress radius contribution
    C = ((struct vessel*) curr_vessel)->CB -
        ((struct vessel*) curr_vessel)->CS * (((struct vessel*) curr_vessel)->bar_tauw /
        ((struct vessel*) curr_vessel)->bar_tauw_h - 1);

    lambda_act = ((struct vessel*) curr_vessel)->a[sn] / ((struct vessel*) curr_vessel)->a_act[sn];

    if (sn == 0) {
        a_act = ((struct vessel*) curr_vessel)->a_act[0];
        C = ((struct vessel*) curr_vessel)->CB;
        lambda_act = 1.0;
    }

    parab_act = 1 - pow((((struct vessel*) curr_vessel)->lambda_m - lambda_act) /
        (((struct vessel*) curr_vessel)->lambda_m - ((struct vessel*) curr_vessel)->lambda_0), 2);

    hat_sigma_act = ((struct vessel*) curr_vessel)->T_act * (1 - exp(-pow(C, 2))) * lambda_act * parab_act;

    sigma_act = ((struct vessel*) curr_vessel)->rhoR_alpha[nts * 1 + sn] / J_s / 
                ((struct vessel*) curr_vessel)->rhoR_h * hat_sigma_act;

    hat_dSdC_act = ((struct vessel*) curr_vessel)->T_act * (pow(lambda_act, -2) / 2 * 
                   ((((struct vessel*) curr_vessel)->lambda_m - lambda_act) / 
                   pow(((struct vessel*) curr_vessel)->lambda_m - ((struct vessel*) curr_vessel)->lambda_0, 2)) 
                   - pow(lambda_act, -3) / 4 * (parab_act));

    Cbar_act = ((struct vessel*) curr_vessel)->rhoR_alpha[nts * 1 + sn] / J_s / 
                ((struct vessel*) curr_vessel)->rhoR_h * 
                lambda_act * lambda_act * lambda_act * lambda_act * hat_dSdC_act;

    //The Lagrange multiplier is the radial stress component
    //subtract from each direction
    lagrange = sigma[0] + ((struct vessel*) curr_vessel)->P/2;
    for (int dir = 0; dir < 3; dir++) {
        
        //Accounting for active stress
        if (dir == 1) {
            sigma[dir] += sigma_act;
            Cbar[dir] += Cbar_act;
        }
        //Calculating stiffness using extra stresses
        Cbar[dir] = 2 * sigma[dir] + 2 * Cbar[dir];
        //Calculating full cauchy stress
        sigma[dir] = sigma[dir] - lagrange;
        
        ((struct vessel*)curr_vessel)->sigma[dir] = sigma[dir];
        ((struct vessel*)curr_vessel)->Cbar[dir] = Cbar[dir];
    }

    //Save updated active radius
    ((struct vessel*) curr_vessel)->a_act[sn] = a_act;

    ((struct vessel*)curr_vessel)->sigma_inv = sigma[0]+sigma[1]+sigma[2];

}

void update_sigma_handshake(void* curr_vessel) {


    //Get current time index
    double s = ((struct vessel*)curr_vessel)->s;
    int sn = ((struct vessel*)curr_vessel)->sn;
    int nts = ((struct vessel*)curr_vessel)->nts;
    double dt = ((struct vessel*)curr_vessel)->dt;
    int taun_min = 0;

    double tau_max = 10000 * (1 / ((struct vessel*)curr_vessel)->k_alpha_h[3]); //max time of 10 half-lives

    //Find the current deformation gradient
    vector<double> F_s(9,0);
    vector<double> R_s = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    vector<double> R_tau = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    //Load in current deformation gradient
    for (int i = 0; i < 9; i++) {
         F_s[i] = ((struct vessel*)curr_vessel)->F_s[sn*9+i];
    }

    vector<double> C_s = mat_mul(mat_trans(F_s), F_s);
    double J_s = mat_det(F_s);

    //Calculate constituent specific stretches for evolving constituents at the current time
    int n_alpha = ((struct vessel*)curr_vessel)->n_alpha;
    vector<double> lambda_alpha_s(n_alpha, 0);
    double eta_alpha = 0;
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Check to see if constituent is isotropic
        eta_alpha = ((struct vessel*)curr_vessel)->eta_alpha_h[alpha];
        if (eta_alpha >= 0) {

            //Stretch is equal to the sqrt of I4
            lambda_alpha_s[alpha] = sqrt(C_s[8] * pow(cos(eta_alpha), 2)
                + C_s[4] * pow(sin(eta_alpha), 2));

            ((struct vessel*)curr_vessel)->lambda_alpha_tau[nts * alpha + sn] = lambda_alpha_s[alpha];
        }

    }

    //Find the mechanical contributions of each constituent for each direction
    //double a, h;
    double J_tau = 1;
    vector<double> F_tau = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    vector<double> G_alpha_h_N = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    vector<double> G_alpha_h = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    vector<double> A_alpha_act = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double lambda_alpha_ntau_s = 0;
    vector<double> F_alpha_ntau_s = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    double hat_S_alpha = 0;
    vector<double> sigma(9,0);
    double lagrange = 0;
    double pol_mod = 0;
    double g_alpha = 1.0;

    //Local active variables
    double C = 0;
    double lambda_act = 0;
    double lambda_act_1 = 0;
    double lambda_act_2 = 0;
    double parab_act = 0;
    double hat_sigma_act = 0, sigma_act = 0, S_act = 0;
    vector<double> sigma_act_mat(9,0);
    double Cbar_act = 0;

    ((struct vessel*) curr_vessel)->lambda_act[sn] = sqrt(C_s[4]);

    //Stiffness variables
    double hat_dS_dlambda2_alpha = 0;
    double dSdC_act = 0;
    vector<double> constitutive_return = { 0, 0 };

    //Full material stiffness tensor
    double CC[3][3][3][3] = {0};
    double CC_act_mat[3][3][3][3] = {0};

    double CC_sym[36] = {0};
    double hat_CC_1[3][3][3][3] = {0};
    double hat_CC_2[3][3][3][3] = {0};

    //Integration variables
    //For mass
    double mq_1 = 0, mq_2 = 0;
    double q_1 = 1.0, q_2 = 1.0;
    double k_1 = 0, k_2 = 0;

    int n = 0; //number of pts in integration interval

    //For stress
    vector<double> hat_sigma_1 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    vector<double> hat_sigma_2 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    //For active stress
    double k_act = ((struct vessel*) curr_vessel)->k_act;
    double q_act_1 = 0, q_act_2 = 0;


    //Boolean for checks
    bool deg_check = 0;

    //Determine if beyond initial time history
    if (s > tau_max) {
        taun_min = sn - int(tau_max / dt);
    }
    else {
        taun_min = 0;
    }

    n = (sn - taun_min) + 1; //number of integration pts

    //Similar Simpson's integration to that used for kinematics
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Trapz rule allows for fast heredity integral evaluation
        k_2 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + sn];
        q_2 = 1.0;
        mq_2 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + sn];

        //Find active radius from current cohort
        if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
            lambda_act_2 = sqrt(C_s[4])/((struct vessel*) curr_vessel)->lambda_act[sn];
            q_act_2 = 1.0;
        }

        //Kinematics
        F_tau = F_s;
        J_tau = J_s;

        eta_alpha = ((struct vessel*)curr_vessel)->eta_alpha_h[alpha];
        g_alpha = ((struct vessel*)curr_vessel)->g_alpha_h[alpha];

        for (int i = 0; i < 3; i++) {
            G_alpha_h[i*3 + i] = ((struct vessel*) curr_vessel)->G_alpha_h[3 * alpha + i];
        }
        R_tau = mat_pol_dec(F_tau);
        G_alpha_h_N = mat_mul(G_alpha_h,R_tau);

        constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, sn);
        hat_S_alpha = constitutive_return[0];
        hat_dS_dlambda2_alpha = constitutive_return[1];
        F_alpha_ntau_s = mat_mul(mat_mul(F_s, mat_inv(F_tau)),G_alpha_h_N);
        hat_sigma_2 = scl_mul(mat_mul(scl_mul(F_alpha_ntau_s, hat_S_alpha),mat_trans(F_alpha_ntau_s)), 1.0/J_s);

        get_material_stiffness(F_alpha_ntau_s, hat_dS_dlambda2_alpha,J_s,hat_CC_2);

        //Boolean for whether the constituent increases ref mass density
        deg_check = ((struct vessel*)curr_vessel)->mR_alpha_h[alpha] > 0;

        //Check if during G&R or at initial time point
        if (sn > 0 && deg_check) {

            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) {

                //Find the 1st intermediate deformation gradient
                for (int i = 0; i < 9; i++) {
                     F_tau[i] = ((struct vessel*)curr_vessel)->F_s[taun*9+i];
                }
                J_tau = mat_det(F_tau);
                R_tau = mat_pol_dec(F_tau);
                G_alpha_h_N = mat_mul(G_alpha_h,R_tau);

                //Find 1st intermediate kinetics
                k_1 = ((struct vessel*)curr_vessel)->k_alpha[nts * alpha + taun];
                q_1 = exp(-(k_2 + k_1) * dt / 2) * q_2;
                mq_1 = ((struct vessel*)curr_vessel)->mR_alpha[nts * alpha + taun] * q_1;

                //Find intermediate active state
                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    lambda_act_1 = sqrt(C_s[4])/((struct vessel*) curr_vessel)->lambda_act[taun];
                    q_act_1 = exp(-k_act * dt) * q_act_2;
                }

                constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, taun);
                hat_S_alpha = constitutive_return[0];
                hat_dS_dlambda2_alpha = constitutive_return[1];
                F_alpha_ntau_s = mat_mul(mat_mul(F_s, mat_inv(F_tau)),G_alpha_h_N);
                hat_sigma_1 = scl_mul(mat_mul(scl_mul(F_alpha_ntau_s, hat_S_alpha),mat_trans(F_alpha_ntau_s)), 1.0/J_s);
                get_material_stiffness(F_alpha_ntau_s, hat_dS_dlambda2_alpha,J_s,hat_CC_1);


                for(int i = 0; i<3; i++){
                    for(int j = 0; j<3; j++){

                        sigma[(i*3)+j] += (mq_2 * hat_sigma_2[(i*3)+j] + mq_1 * hat_sigma_1[(i*3)+j])
                            / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * dt / 2;

                        for(int k = 0; k<3; k++){
                            for(int l = 0; l<3; l++){

                                CC[i][j][k][l] += (mq_2 * hat_CC_2[i][j][k][l] + mq_1 * hat_CC_1[i][j][k][l])
                                    / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * dt / 2;
                            }
                        }
                    }
                }

                //Store active vars for next iteration
                //Find intermediate active state

                if (((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
                    lambda_act += k_act * (q_act_2 * lambda_act_2 + q_act_1 * lambda_act_1) * dt / 2;
                    lambda_act_2 = lambda_act_1;
                    q_act_2 = q_act_1;
                }

                //Store intermediate kinetics for next iteration
                k_2 = k_1;
                q_2 = q_1;
                mq_2 = mq_1;

                //Store intermediate stress and stiffness for next iteration
                hat_sigma_2 = hat_sigma_1;

                for(int i = 0; i<3; i++){
                    for(int j = 0; j<3; j++){
                        for(int k = 0; k<3; k++){
                            for(int l = 0; l<3; l++){
                                hat_CC_2[i][j][k][l] = hat_CC_1[i][j][k][l];
                            }
                        }
                    }
                }


            }

            //Add in the stress and stiffness contributions of the initial material
            if (taun_min == 0) {
                for(int i = 0; i<3; i++){
                    for(int j = 0; j<3; j++){
                        sigma[(i*3)+j] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + 0]
                            / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * q_1 * hat_sigma_1[(i*3)+j];
                        for(int k = 0; k<3; k++){
                            for(int l = 0; l<3; l++){
                                CC[i][j][k][l] +=((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + 0]
                                    / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * q_1 * hat_CC_1[i][j][k][l];
                            }
                        }
                    }
                }

            }

        }
        //Initial time point and constituents with prescribed degradation profiles
        else {
            //Find stress from initial cohort
            constitutive_return = constitutive(curr_vessel, lambda_alpha_s[alpha], alpha, 0);
            hat_S_alpha = constitutive_return[0];
            hat_dS_dlambda2_alpha = constitutive_return[1];

            F_alpha_ntau_s = mat_mul(F_s,G_alpha_h);

            //Check if anisotropic
            hat_sigma_2 = scl_mul(mat_mul(scl_mul(F_alpha_ntau_s, hat_S_alpha),mat_trans(F_alpha_ntau_s)), 1.0/J_s);

            get_material_stiffness(F_alpha_ntau_s, hat_dS_dlambda2_alpha,J_s,hat_CC_2);

            for(int i = 0; i<3; i++){
                for(int j = 0; j<3; j++){

                    sigma[(i*3)+j] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + sn] /
                        ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * hat_sigma_2[(i*3)+j];

                    for(int k = 0; k<3; k++){
                        for(int l = 0; l<3; l++){                 

                            CC[i][j][k][l] += ((struct vessel*)curr_vessel)->rhoR_alpha[nts * alpha + sn] /
                                ((struct vessel*)curr_vessel)->rho_hat_alpha_h[alpha] * hat_CC_2[i][j][k][l];

                        }
                    }
                }
            }
        }


        if (taun_min == 0 && ((struct vessel*) curr_vessel)->alpha_active[alpha] == 1) {
            if (!deg_check){
                q_act_1 = exp(-k_act * dt) * q_act_2;
            }
            lambda_act += sqrt(C_s[4])/((struct vessel*) curr_vessel)->lambda_act[0] * q_act_1;
        }

    }
        
    //Find active stress contribtion
    //add in initial active stress radius contribution
    C = ((struct vessel*) curr_vessel)->CB -
        ((struct vessel*) curr_vessel)->CS * (((struct vessel*) curr_vessel)->bar_tauw /
        ((struct vessel*) curr_vessel)->bar_tauw_h - 1);

    if (sn == 0) {
        C = ((struct vessel*) curr_vessel)->CB;
        lambda_act = 1.0;
    }

    parab_act = 1 - pow((((struct vessel*) curr_vessel)->lambda_m - lambda_act) /
        (((struct vessel*) curr_vessel)->lambda_m - ((struct vessel*) curr_vessel)->lambda_0), 2);

    //hat_sigma_act = ((struct vessel*) curr_vessel)->T_act * (1 - exp(-pow(C, 2))) * lambda_act * parab_act;

    S_act = ((struct vessel*) curr_vessel)->T_act * (1 - exp(-pow(C, 2))) * pow(lambda_act, -1) * parab_act *
            ((struct vessel*)curr_vessel)->rhoR_alpha[nts * 1 + sn] / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[1];


    dSdC_act = ((struct vessel*) curr_vessel)->T_act *  (1 - exp(-pow(C, 2))) * (1/lambda_act) *
               ((-1/pow(lambda_act,2))*parab_act + 2*(1/lambda_act)*(((struct vessel*) curr_vessel)->lambda_m - lambda_act)/ pow(((struct vessel*) curr_vessel)->lambda_m - ((struct vessel*) curr_vessel)->lambda_0, 2)) *
               ((struct vessel*)curr_vessel)->rhoR_alpha[nts * 1 + sn] / ((struct vessel*)curr_vessel)->rho_hat_alpha_h[1];

    //Radial stress estimated to be (p_int - p_ext)/2
    lagrange =  0.0; //-(((struct vessel*) curr_vessel)->P)/2.0 - mat_ddot(sigma, mat_mul(mat_mul(R_s,{1.0,0,0,0,0,0,0,0,0}),mat_trans(R_s)));

    int iv[6] = {0,1,2,0,1,2};
    int jv[6] = {0,1,2,1,2,0};

    sigma_act_mat = scl_mul(mat_mul(mat_mul(F_s,{0,0,0,0,S_act,0,0,0,0}),mat_trans(F_s)),1.0/J_s);
    get_active_material_stiffness(dSdC_act,CC_act_mat,F_s,J_s);

    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            sigma[i*3 + j] += sigma_act_mat[i*3 + j];
        }
    }

    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    CC[i][j][k][l] += CC_act_mat[i][j][k][l];
                }
            }
        }
    }

    //Enforcing minor symmetry and subtract lagrange
    for(int i = 0; i<3; i++){
        sigma[(i*3)+i] = sigma[(i*3)+i] + lagrange;
        for(int j = 0; j<3; j++){
            for(int k = 0; k<3; k++){
                for(int l = 0; l<3; l++){
                    CC[i][j][k][l] = 0.5*(CC[i][j][k][l]+CC[i][j][l][k]);
                    CC[i][j][l][k] = CC[i][j][k][l];
                }
            }
        }
    }

    for(int i = 0; i<9; i++){
        ((struct vessel*)curr_vessel)->sigma[i] = sigma[i];
    }

    for(int i = 0; i<6; i++){
        for(int j = 0; j<6; j++){
            CC_sym[i*6 + j] = CC[iv[i]] [jv[i]] [iv[j]] [jv[j]];
            ((struct vessel*)curr_vessel)->CC[i*6 + j] = CC_sym[i*6 + j];
        }
    }

}

vector<double> constitutive(void* curr_vessel, double lambda_alpha_s, int alpha, int ts) {

    double lambda_alpha_ntau_s = 0;
    double Q1 = 0;
    double Q2 = 0;
    double hat_S_alpha = 0;
    double hat_dS_dlambda2_alpha = 0;
    double pol_mod = 0;
    double epsilon_curr = 0;
    int nts = ((struct vessel*)curr_vessel)->nts;
    int sn = ((struct vessel*)curr_vessel)->sn;
    vector<double> return_constitutive = { 0, 0 };

    //Check if ansisotropic
    if (((struct vessel*)curr_vessel)->eta_alpha_h[alpha] >= 0) {

        lambda_alpha_ntau_s = ((struct vessel*)curr_vessel)->g_alpha_h[alpha] *
                              lambda_alpha_s / ((struct vessel*)curr_vessel)->lambda_alpha_tau[nts * alpha + ts];

        if (lambda_alpha_ntau_s < 1) {
            lambda_alpha_ntau_s = 1;
        }

        Q1 = (pow(lambda_alpha_ntau_s, 2) - 1);
        Q2 = ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha + 1] * pow(Q1, 2);
        hat_S_alpha = ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha] * Q1 * exp(Q2);
        hat_dS_dlambda2_alpha = ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha] * exp(Q2) * (1 + 2 * Q2);

    }
    else {

        if (alpha < ((struct vessel*)curr_vessel)->n_pol_alpha) {

            if (((struct vessel*)curr_vessel)->epsilon_alpha[nts * alpha + sn] <
                ((struct vessel*)curr_vessel)->epsilon_pol_min[alpha]) {
                ((struct vessel*)curr_vessel)->epsilon_pol_min[alpha] = ((struct vessel*)curr_vessel)->epsilon_alpha[nts * alpha + sn];
            }
            epsilon_curr = ((struct vessel*)curr_vessel)->epsilon_pol_min[alpha];
            pol_mod = 0.03 * pow(epsilon_curr, 2);
        }
        else {
            pol_mod = 1;
        }

        hat_S_alpha = pol_mod * ((struct vessel*)curr_vessel)->c_alpha_h[2 * alpha];
    }

    return_constitutive = { hat_S_alpha , hat_dS_dlambda2_alpha };
    return return_constitutive;

}

double get_app_visc(void* curr_vessel, int sn){
    //Returns apparent viscosity if diameter dependent factors from Secomb 2017 are in effect
    //Otherwise returns default for blood
    double d = 0.0;
    double mu = 0.00;

    if (((struct vessel*)curr_vessel)->app_visc_flag == 1){
        d = ((struct vessel*)curr_vessel)->a[sn] * 2 * 1000000;
	    mu = (1+(6*exp(-0.0858*d)+3.2-2.44*exp(-0.06*pow(d,0.645))-1)*pow(d/(d-1.1),2)*pow(d/(d-1.1),2)) * 0.0124;
    }
    else{
        mu = 0.04;
    }

    return mu;
}
