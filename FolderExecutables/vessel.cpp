// vessel.h
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>


#include "vessel.h"
#include "functions.h"


using std::string;
using std::vector;
using std::cout;


template <typename T>
// C++ template to print vessel contents
std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
    for (int i = 0; i < v.size(); ++i) { 
        os << v[i]; 
        if (i != v.size() - 1) 
            os << " "; 
    }
    return os; 
}

template <typename T>
// C++ template to read vector container elements
std::ifstream& operator>>(std::ifstream& in, vector<T>& v) {
    T x;
    char next;
    int i=0;
    while(in.get(next))
    {
        if (next == '\n')  
        {    break; }
    }
    while ((in.peek()!='\n') && (in>>x))
    {
        v[i]=x;
        i++;
    }
    return in; 
}

// C++ template to print vector container elements 
std::ostream& operator<<(std::ostream& os, vessel const &v) {
    os << v.nts << "\n";
    os << v.dt << "\n";
    os << v.sn << "\n";
    os << v.s << "\n";
    os << v.A_h << "\n";
    os << v.B_h << "\n";
    os << v.H_h << "\n";
    os << v.A_mid_h << "\n";
    os << v.a_h << "\n";
    os << v.b_h << "\n";
    os << v.h_h << "\n";
    os << v.a_mid_h << "\n";
    os << v.lambda_z_h << "\n";
    os << v.a << "\n";
    os << v.a_mid << "\n";
    os << v.b << "\n";
    os << v.h << "\n";
    os << v.A << "\n";
    os << v.A_mid << "\n";
    os << v.B << "\n";
    os << v.H << "\n";
    os << v.lambda_z_pre << "\n";
    os << v.n_alpha << "\n";
    os << v.n_pol_alpha << "\n";
    os << v.n_native_alpha << "\n";
    os << v.alpha_infl << "\n";
    os << v.c_alpha_h << "\n";
    os << v.eta_alpha_h << "\n";
    os << v.g_alpha_h << "\n";
    os << v.G_alpha_h << "\n";
    os << v.phi_alpha_h << "\n";
    os << v.rhoR_alpha_h << "\n";
    os << v.mR_alpha_h << "\n";
    os << v.k_alpha_h << "\n";
    os << v.K_sigma_p_alpha_h << "\n";
    os << v.K_sigma_d_alpha_h << "\n";
    os << v.K_tauw_p_alpha_h << "\n";
    os << v.K_tauw_d_alpha_h << "\n";
    os << v.rho_hat_alpha_h << "\n";
    os << v.epsilonR_alpha_0 << "\n";
    os << v.rhoR_h << "\n";
    return os;
}

vessel::vessel() { //Default constructor
    vessel_name = "default";
    file_name = "Vs_out";
    gnr_name = "GnR_out";
    equil_gnr_name = "Equil_GnR_out";
    exp_name = "Exp_out";
    hns_name = "HnS_out";

    //Time variables
    nts = 0; //total number of time steps
    dt = 0; //time increment
    sn = 0; //current time step index
    s = 0; //actual current time

    //Geometric quantities
    A_h = 0, B_h = 0, H_h = 0; //Traction-free reference
    A_mid_h = 0; //Tf Midpoint reference
    a_h = 0, b_h = 0, h_h = 0; //In-vivo reference
    a_mid_h = 0; //In vivo midpoint
    lambda_z_h = 0; //Reference in vivo stretch

    a = { 0 }, a_mid = { 0 }, b = { 0 }, h = { 0 }; //Evolving in vivo reference
    a_pas = { 0 }, a_mid_pas = { 0 }, h_pas = { 0 }; //Evolving passive configuration
    A = { 0 }, A_mid = { 0 }, B = { 0 }, H = { 0 }, lambda_z_pre = { 0 }; //Evolving traction free

    //Number of constituents
    n_alpha = 0;
    n_pol_alpha = 0;
    n_native_alpha = 0;

    //Constituent inflammatory status
    alpha_infl = { 0 };

    //Material properties
    c_alpha_h = { 0 }, eta_alpha_h = { 0 }, g_alpha_h = { 0 }, G_alpha_h = { 0 };
    epsilon_pol_min = { 0 };

    //Mass fractions, referential apparent densities, kinetic quantities
    phi_alpha_h = { 0 }, rhoR_alpha_h = { 0 }, mR_alpha_h = { 0 }, k_alpha_h = { 0 };
    K_sigma_p_alpha_h = { 0 }, K_sigma_d_alpha_h = { 0 }, K_tauw_p_alpha_h = { 0 }, K_tauw_d_alpha_h = { 0 };

    //True mass densities, initial volume fractions
    rho_hat_alpha_h = { 0 }, epsilonR_alpha_0 = { 0 };

    //Reference true mass density
    rhoR_h = 0;

    //Histories
    rhoR = { 0 }, rho = { 0 }, rhoR_alpha = { 0 }, mR_alpha = { 0 }, k_alpha = { 0 };
    epsilonR_alpha = { 0 }, epsilon_alpha = { 0 };
    ups_infl_p = { 0 }, ups_infl_d = { 0 };

    //Reference loading quantities
    P_h = 0, f_h = 0, bar_tauw_h = 0, Q_h = 0, P_prev=0, sigma_inv_h = 0;
    sigma_h = { 0 };

    //Current loading quantities
    lambda_th_curr = 0, lambda_z_curr = 0;
    P = 0, f = 0, bar_tauw = 0, Q = 0, sigma_inv = 0;
    sigma = { 0 }, Cbar = { 0 }, CC = { 0 }, lambda_alpha_tau = { 0 }, lambda_z_tau = { 0 };
    mb_equil = 0; //Current mechanobiological equil. state

    //Active stress quantities
    alpha_active = { 0 }; //boolean for active constituents
    a_act = { 0 }; //active radius history
    lambda_act = { 0 }; //active radius history
    T_act = 0, T_act_h = 0; //Max active stress, homeostatic max active stress
    k_act = 0; //Active remodeling parameters
    lambda_0 = 0; //Min contractile stretch
    lambda_m = 0; //Max contractile stretch
    CB = 0; //Basal VC to VD ratio
    CS = 0; //Scaling factor for VC to VD ratio

    //Mechanobiologically equilibrated quantities
    a_e = 0; //equilibrated radius
    h_e = 0; //equilibrated thickness
    rho_c_e = 0; //equilibrated collagen density
    rho_m_e = 0; //equilibrated smc density
    f_z_e = 0; //equilibrated axial force
    mb_equil_e = 0; //Current mechanobiological equil. state

    //Phenomenologic immune parameters
    Ki_p_h = 0, Ki_d_h = 0;

    //Flags
    num_exp_flag = 0; //indicates whether doing reg G&R step or a numerical experiment
    pol_only_flag = 0; //indicates whether other constituents are produced
    wss_calc_flag = 0; //indicates if GnR should update its own current WSS
    app_visc_flag = 0; //indicates whether to use the empirical correction for viscosity from Secomb 2017

}

//Initialize the reference vessel for the simulation    
void  vessel::initializeVesselHandshake(char* prefix_char, char* name_char, int num_days, double step_size,double anysm_arg, double tevg_arg){

    string prefix_arg = string(prefix_char);
    string name_arg = string(name_char);

    //Initialize the reference vessel for the simulation
    string native_file = "FolderVesselConfigurationFiles/Native_in_handshake_";// + name_arg;
    string immune_file = "FolderVesselConfigurationFiles/Immune_in_";// + name_arg;
    string scaffold_file = "FolderVesselConfigurationFiles/Scaffold_in_";// + name_arg;

    //Initialize TEVG files if necessary
    if(tevg_arg > 0.0) {
        //Initialize TEVG
        initializeNativeHandshake(native_file,num_days,step_size);
        initializeTEVGHandshake(scaffold_file,immune_file,tevg_arg,num_days,step_size);
    } else {
        initializeNativeHandshake(native_file,num_days,step_size);
    }

    //Get all other input arguements and apply to TEVG
    gnr_name = prefix_arg + "/" + gnr_name + "_" + name_arg;
    exp_name = prefix_arg + "/" + exp_name + "_" + name_arg;
    file_name = prefix_arg + "/" + file_name + "_" + name_arg;
    hns_name = prefix_arg + "/" + hns_name + "_" + name_arg;

    //------------------------------------------------------------------------

    //For elastin degradation 
    double Q_p1;

    if(tevg_arg <= 0.0) {
        for (int sn = 1; sn < nts; sn++) {
            //Calculate elastin degradation
            //s = sn * simulation_vessel.dt;
            
            Q_p1 = 1.0;
            
            epsilonR_alpha[0 * nts + sn] = Q_p1 * epsilonR_alpha[0 * nts + 0];

            rhoR_alpha[0 * nts + sn] = epsilonR_alpha[0 * nts + sn] * 
                                                                rho_hat_alpha_h[0];
        }
    }
    if(anysm_arg > 0.0) {
        //Change endothelial functioning to be proportional to damage
        for (int alpha = 0; alpha < n_alpha; alpha++) {
            K_tauw_p_alpha_h[alpha] = std::max(K_tauw_p_alpha_h[alpha] - anysm_arg, 0.0);
        }
        
        c_alpha_h[0] = c_alpha_h[0]*(1.0 - anysm_arg);
    }

}


int vessel::updateVesselHandshake(int restart_arg, int iter_arg, double sigma_arg, double tauw_arg, double * F, double * out_array) {


    try {
        if (sn == 0){
            //Write initial state to file
            update_sigma_handshake(this);
            update_kinetics(*this);
        }

        sigma_inv = sigma_arg;
        bar_tauw = tauw_arg;
        F_curr[0] = F[0];
        F_curr[1] = F[1];
        F_curr[2] = F[2];
        F_curr[3] = F[3];
        F_curr[4] = F[4];
        F_curr[5] = F[5];
        F_curr[6] = F[6];
        F_curr[7] = F[7];
        F_curr[8] = F[8];

        update_time_step_handshake(*this, iter_arg);

        out_array[0] = s;   
        for (int i = 0; i < 36; i++) {
            out_array[i+1] = CC[i];
        }
        for (int i = 0; i < 9; i++){
            out_array[i+37] = sigma[i];
        }
        out_array[46] = rhoR[sn];
        out_array[47] = rho[sn];
        for (int i = 0; i < 9; i++){
            out_array[i+48] = F_s[9*sn+i];
        }
        out_array[57] = sigma_inv;
        out_array[58] = bar_tauw;


    } catch(std::exception& e) {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;
}

//Initialize the reference vessel for the simulation    
void vessel::initializeNative(string native_name, double n_days_inp, double dt_inp) {

    //Set native vessel time parameters
    double n_days = n_days_inp;
    dt = dt_inp;
    nts = int(n_days / dt);
    sn = 0;
    s = 0.0;

    double mu = 0; //apparent viscosity

    //Load the expereimentally determined and prescribed properties of the vessel from file
    //Input arguments for scaffold input file (ELS)
    std::ifstream native_in(native_name);

    native_in >> vessel_name; //Type of vessel simulated

    //Initilize the parameters for the reference native vessel
    //Geometric parameters
    //Homeostatic parameters are those for a NATIVE vessel
    native_in >> a_h; //in vivo reference inner radius
    a_h = a_h * mm_to_m;
    native_in >> h_h; //in vivo reference medial thickness
    h_h = h_h * mm_to_m;
    a_mid_h = a_h + h_h / 2; //in vivo reference wall mid-point
    //curr_vessel.b_h = curr_vessel.a_h + curr_vessel.h_h; //in vivo reference outer diameter

    //Histories for geometric parameters
    a_mid.resize(nts);
    a.resize(nts);
    h.resize(nts);
    a_mid_pas.resize(nts);
    a_pas.resize(nts);
    h_pas.resize(nts);
    a_mid[0] = a_mid_h;
    a[0] = a_h;
    h[0] = h_h;
    a_mid_pas[0] = a_mid_h; 
    a_pas[0] = a_h;
    h_pas[0] = h_h;

    //Axial stretch with in vivo reference
    native_in >> lambda_z_h; //Reference in vivo stretch    

    //Constituent material properties
    n_alpha = 6; //number of constituents alpha
    n_pol_alpha = 0; //number of polymeric constituents
    alpha_infl = { 0, 0, 0, 0, 0, 0 }; //flags for inflammation

    //Medial and adventitial constituents have the same material behavior
    native_in >> c1_e >> c2_e >> c1_m >> c2_m >> c1_ct >> c2_ct;
    native_in >> c1_cz >> c2_cz >> c1_cd1 >> c2_cd1 >> c1_cd2 >> c2_cd2;
    c_alpha_h = { c1_e, c2_e, c1_m, c2_m, c1_ct, c2_ct,
                                c1_cz, c2_cz, c1_cd1, c2_cd1, c1_cd2, c2_cd2 };

    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents
    native_in >> eta_e_h >> eta_m_h >> eta_ct_h >> eta_cz_h >> eta_cd1_h >> eta_cd2_h;
    eta_alpha_h = { eta_e_h * M_PI / 180.0, eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0,
                                  eta_cz_h * M_PI / 180.0, eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0 };

    //Pre-stretch parameters
    native_in >> G_e_rh >> G_e_th >> G_e_zh;
    native_in >> g_e_h >> g_m_h >> g_ct_h >> g_cz_h >> g_cd1_h >> g_cd2_h;
    g_alpha_h = { g_e_h, g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h };

    //Mass density parameters
    native_in >> rho_hat_h; //density of tissue as a whole kg/m^3
    rhoR_h = rho_hat_h;
    rho[0] = rho_hat_h;
    rhoR[0] = rho_hat_h;
    rho_hat_alpha_h = { rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h };

    //Homeostatic mass fractions of constituents
    native_in >> phi_e_h >> phi_m_h >> phi_ct_h >> phi_cz_h >> phi_cd1_h >> phi_cd2_h;
    phi_alpha_h = { phi_e_h, phi_m_h, phi_ct_h, phi_cz_h, phi_cd1_h, phi_cd2_h };

    //Homeostatic mass densities
    rhoR_e_h = phi_e_h * rho_hat_h, rhoR_m_h = phi_m_h * rho_hat_h,
        rhoR_ct_h = phi_ct_h * rho_hat_h, rhoR_cz_h = phi_cz_h * rho_hat_h,
        rhoR_cd1_h = phi_cd1_h * rho_hat_h, rhoR_cd2_h = phi_cd2_h * rho_hat_h;
    rhoR_alpha_h = { rhoR_e_h, rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h };

    //Kinetic parameters
    //Degradation parameters
    native_in >> k_e_h >> k_m_h >> k_ct_h >> k_cz_h >> k_cd1_h >> k_cd2_h;
    k_alpha_h = { k_e_h, k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h };

    //Initialize kinetic arrays
    rho.resize(nts); //Current mass density history
    rhoR.resize(nts); //Referential mass density history
    rhoR_alpha.resize(nts * n_alpha); //referential apparent mass densities (time history)
    epsilon_alpha.resize(nts * n_alpha); //current volume fractions
    epsilonR_alpha.resize(nts * n_alpha); //referential volume fractions
    mR_alpha.resize(nts * n_alpha); //referential mass production rate (time history)
    k_alpha.resize(nts * n_alpha); //mass removal decay (time history)
    mR_alpha_h.resize(n_alpha);
    rhoR_alpha_h.resize(n_alpha);
    lambda_alpha_tau.resize(nts * n_alpha);

    //Gain parameters
    //Stress mediated production
    native_in >> K_sigma_p_e_h >> K_sigma_p_m_h >> K_sigma_p_ct_h >> K_sigma_p_cz_h >> K_sigma_p_cd1_h >> K_sigma_p_cd2_h;
    K_sigma_p_alpha_h = { K_sigma_p_e_h, K_sigma_p_m_h, K_sigma_p_ct_h,
                                        K_sigma_p_cz_h, K_sigma_p_cd1_h, K_sigma_p_cd2_h };

    //Stress mediated degradation
    native_in >> K_sigma_d_e_h >> K_sigma_d_m_h >> K_sigma_d_ct_h >> K_sigma_d_cz_h >> K_sigma_d_cd1_h >> K_sigma_d_cd2_h;
    K_sigma_d_alpha_h = { K_sigma_d_e_h, K_sigma_d_m_h, K_sigma_d_ct_h,
                                        K_sigma_d_cz_h, K_sigma_d_cd1_h, K_sigma_d_cd2_h };

    //Wall Shear Stress mediated production
    native_in >> K_tauw_p_e_h >> K_tauw_p_m_h >> K_tauw_p_ct_h >> K_tauw_p_cz_h >> K_tauw_p_cd1_h >> K_tauw_p_cd2_h;
    K_tauw_p_alpha_h = { K_tauw_p_e_h, K_tauw_p_m_h, K_tauw_p_ct_h,
                                        K_tauw_p_cz_h, K_tauw_p_cd1_h, K_tauw_p_cd2_h };

    //Wall Shear Stress mediated degradation
    native_in >> K_tauw_d_e_h >> K_tauw_d_m_h >> K_tauw_d_ct_h >> K_tauw_d_cz_h >> K_tauw_d_cd1_h >> K_tauw_d_cd2_h;
    K_tauw_d_alpha_h = { K_tauw_d_e_h, K_tauw_d_m_h, K_tauw_d_ct_h,
                                        K_tauw_d_cz_h, K_tauw_d_cd1_h, K_tauw_d_cd2_h };

    //Initialize homeostatic loading variables
    //Pressure and Flow
    native_in >> P_h >> Q_h;
    //Set WSS reference in dynes/cm2
    mu = get_app_visc(this, sn);
    bar_tauw_h = 4*mu*Q_h/(3.14159265*pow(a_h*100, 3));

    //Initial stresses in each direction
    sigma_h.resize(3);
    sigma.resize(3);
    Cbar.resize(3);
    CC.resize(36);

    //Active stress parameters
    alpha_active = { 0, 1, 0, 0, 0, 0 };
    a_act.resize(nts);
    a_act[0] = a_h;
    lambda_act.resize(nts);
    lambda_act[0] = 1.0;
    native_in >> k_act >> lambda_0 >> lambda_m; //active remodelling time, min active stretch, max active stretch
    native_in >> CB >> CS; //vasodilator ratios
    native_in >> T_act_h; //homeostatic active stress magnitude
    T_act = T_act_h; //current active stress 

    //Initializion with looping through constituents
    double rhoR_0 = 0.0;
    double eta = 0;
    double g_alpha = 0;
    G_alpha_h.resize(3 * n_alpha);
    
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Initialize the deposition tensor
        eta = eta_alpha_h[alpha];
        g_alpha = g_alpha_h[alpha];
        if (eta >= 0) { //for anisotropic constituents
            G_alpha_h[3 * alpha] = 0.0;
            G_alpha_h[3 * alpha + 1] = g_alpha * sin(eta);
            G_alpha_h[3 * alpha + 2] = g_alpha * cos(eta);
        }
        else { //for isotropic constituents (i.e. elastin)
            G_alpha_h[3 * alpha] = G_e_rh;
            G_alpha_h[3 * alpha + 1] = G_e_th;
            G_alpha_h[3 * alpha + 2] = G_e_zh;
        }

        //Initialize homeostatic mass productions
        mR_alpha_h[alpha] = k_alpha_h[alpha] * rhoR_alpha_h[alpha];

        //Initialize simulation mass kinetics to hemeostatic values
        mR_alpha[nts * alpha] = mR_alpha_h[alpha];
        k_alpha[nts * alpha] = k_alpha_h[alpha];

        //Initialize volume fractions, which are equal to mass fractions if all constituents have same true mass density
        epsilonR_alpha[nts * alpha] = phi_alpha_h[alpha];
        epsilon_alpha[nts * alpha] = phi_alpha_h[alpha];

        //Initialize density time histories to their homeostatic values for native vessel
        rhoR_alpha[nts * alpha] = rhoR_alpha_h[alpha];

        //Initilize stretch histories
        lambda_alpha_tau[nts * alpha] = 1.0;
    }

    //Solve for native stress state
    //Initialize variables for load history
    P = P_h;
    bar_tauw = bar_tauw_h;
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0; 
    lambda_z_tau.resize(nts);
    lambda_z_tau[0] = 1.0;

    //NOTE: If input homeostatic configuration is not close to actual homeostatic configuration
    //this will give an erroneous 1st passive config. 
    //Solve for native passive configuration
    int equil_check = 0;
    int num_act_incr = 10;
    num_exp_flag = 1;
    for (int i = 1; i < num_act_incr; i++){
        T_act = T_act_h * (1 - static_cast<double>(i) / num_act_incr);
        equil_check = find_iv_geom(this);
    }
    num_exp_flag = 0;
    //printf("%s %f %s %f\n", "pas inner radius: ", a[0], "pas thickness: ", h[0]);
    //fflush(stdout);
    
    //Assign passive geometry parameters
    h_pas[0] = h[0];
    a_pas[0] = a[0];
    a_mid_pas[0] = a_mid[0];
    //Reinitialize to input geometry parameters for solving
    a[0] = a_h;
    h[0] = h_h;
    a_mid[0] = a_mid_h;

    //Solve for native active (homeostatic) configuration
    T_act = T_act_h;
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;
    num_exp_flag = 0;
    equil_check = find_iv_geom(this);
    //printf("%s %f %s %f\n", "act inner radius: ", a[0], "act thickness: ", h[0]);
    //fflush(stdout);

    //Assign in vivo (active geometry parameters)
    h_h = h[0];
    a_h = a[0];
    a_mid_h = a_mid[0];

    //Update stress to in vivo stress for homeostatic stress calc
    update_sigma(this);

    //Set reference values
    sigma_h = sigma;
    sigma_inv_h = sigma_h[0]+sigma_h[1]+sigma_h[2];
    mu = get_app_visc(this, sn);
    bar_tauw_h = 4 * mu * Q_h / (3.14159265 * pow(a_h * 100, 3));

    P = P_h;
    bar_tauw = bar_tauw_h;
    Q = Q_h;

    //For pressure ramping
    P_prev = P;
    native_in.close();

}


//Initialize the reference vessel for the simulation    
void vessel::initializeNativeExplicit(string native_name, double n_days_inp, double dt_inp) {

    //Set native vessel time parameters
    double n_days = n_days_inp;
    dt = dt_inp;
    nts = int(n_days / dt);
    sn = 0;
    s = 0.0;

    double mu = 0; //apparent viscosity

    //Load the expereimentally determined and prescribed properties of the vessel from file
    //Input arguments for scaffold input file (ELS)
    std::ifstream native_in(native_name);

    native_in >> vessel_name; //Type of vessel simulated

    //Initilize the parameters for the reference native vessel
    //Geometric parameters
    //Homeostatic parameters are those for a NATIVE vessel
    native_in >> a_h; //in vivo reference inner radius
    a_h = a_h * mm_to_m;
    native_in >> h_h; //in vivo reference medial thickness
    h_h = h_h * mm_to_m;
    a_mid_h = a_h + h_h / 2; //in vivo reference wall mid-point
    //curr_vessel.b_h = curr_vessel.a_h + curr_vessel.h_h; //in vivo reference outer diameter

    //Histories for geometric parameters
    a_mid.resize(nts);
    a.resize(nts);
    h.resize(nts);
    a_mid_pas.resize(nts);
    a_pas.resize(nts);
    h_pas.resize(nts);
    a_mid[0] = a_mid_h;
    a[0] = a_h;
    h[0] = h_h;
    a_mid_pas[0] = a_mid_h; 
    a_pas[0] = a_h;
    h_pas[0] = h_h;

    //Axial stretch with in vivo reference
    native_in >> lambda_z_h; //Reference in vivo stretch    

    //Constituent material properties
    n_alpha = 6; //number of constituents alpha
    n_pol_alpha = 0; //number of polymeric constituents
    alpha_infl = { 0, 0, 0, 0, 0, 0 }; //flags for inflammation

    //Medial and adventitial constituents have the same material behavior
    native_in >> c1_e >> c2_e >> c1_m >> c2_m >> c1_ct >> c2_ct;
    native_in >> c1_cz >> c2_cz >> c1_cd1 >> c2_cd1 >> c1_cd2 >> c2_cd2;
    c_alpha_h = { c1_e, c2_e, c1_m, c2_m, c1_ct, c2_ct,
                                c1_cz, c2_cz, c1_cd1, c2_cd1, c1_cd2, c2_cd2 };

    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents
    native_in >> eta_e_h >> eta_m_h >> eta_ct_h >> eta_cz_h >> eta_cd1_h >> eta_cd2_h;
    eta_alpha_h = { eta_e_h * M_PI / 180.0, eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0,
                                  eta_cz_h * M_PI / 180.0, eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0 };

    //Pre-stretch parameters
    native_in >> G_e_rh >> G_e_th >> G_e_zh;
    native_in >> g_e_h >> g_m_h >> g_ct_h >> g_cz_h >> g_cd1_h >> g_cd2_h;
    g_alpha_h = { g_e_h, g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h };

    //Mass density parameters
    native_in >> rho_hat_h; //density of tissue as a whole kg/m^3
    rhoR_h = rho_hat_h;
    rho[0] = rho_hat_h;
    rhoR[0] = rho_hat_h;
    rho_hat_alpha_h = { rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h };

    //Homeostatic mass fractions of constituents
    native_in >> phi_e_h >> phi_m_h >> phi_ct_h >> phi_cz_h >> phi_cd1_h >> phi_cd2_h;
    phi_alpha_h = { phi_e_h, phi_m_h, phi_ct_h, phi_cz_h, phi_cd1_h, phi_cd2_h };

    //Homeostatic mass densities
    rhoR_e_h = phi_e_h * rho_hat_h, rhoR_m_h = phi_m_h * rho_hat_h,
        rhoR_ct_h = phi_ct_h * rho_hat_h, rhoR_cz_h = phi_cz_h * rho_hat_h,
        rhoR_cd1_h = phi_cd1_h * rho_hat_h, rhoR_cd2_h = phi_cd2_h * rho_hat_h;
    rhoR_alpha_h = { rhoR_e_h, rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h };

    //Kinetic parameters
    //Degradation parameters
    native_in >> k_e_h >> k_m_h >> k_ct_h >> k_cz_h >> k_cd1_h >> k_cd2_h;
    k_alpha_h = { k_e_h, k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h };

    //Initialize kinetic arrays
    rho.resize(nts); //Current mass density history
    rhoR.resize(nts); //Referential mass density history
    rhoR_alpha.resize(nts * n_alpha); //referential apparent mass densities (time history)
    epsilon_alpha.resize(nts * n_alpha); //current volume fractions
    epsilonR_alpha.resize(nts * n_alpha); //referential volume fractions
    mR_alpha.resize(nts * n_alpha); //referential mass production rate (time history)
    k_alpha.resize(nts * n_alpha); //mass removal decay (time history)
    mR_alpha_h.resize(n_alpha);
    rhoR_alpha_h.resize(n_alpha);
    lambda_alpha_tau.resize(nts * n_alpha);

    //Gain parameters
    //Stress mediated production
    native_in >> K_sigma_p_e_h >> K_sigma_p_m_h >> K_sigma_p_ct_h >> K_sigma_p_cz_h >> K_sigma_p_cd1_h >> K_sigma_p_cd2_h;
    K_sigma_p_alpha_h = { K_sigma_p_e_h, K_sigma_p_m_h, K_sigma_p_ct_h,
                                        K_sigma_p_cz_h, K_sigma_p_cd1_h, K_sigma_p_cd2_h };

    //Stress mediated degradation
    native_in >> K_sigma_d_e_h >> K_sigma_d_m_h >> K_sigma_d_ct_h >> K_sigma_d_cz_h >> K_sigma_d_cd1_h >> K_sigma_d_cd2_h;
    K_sigma_d_alpha_h = { K_sigma_d_e_h, K_sigma_d_m_h, K_sigma_d_ct_h,
                                        K_sigma_d_cz_h, K_sigma_d_cd1_h, K_sigma_d_cd2_h };

    //Wall Shear Stress mediated production
    native_in >> K_tauw_p_e_h >> K_tauw_p_m_h >> K_tauw_p_ct_h >> K_tauw_p_cz_h >> K_tauw_p_cd1_h >> K_tauw_p_cd2_h;
    K_tauw_p_alpha_h = { K_tauw_p_e_h, K_tauw_p_m_h, K_tauw_p_ct_h,
                                        K_tauw_p_cz_h, K_tauw_p_cd1_h, K_tauw_p_cd2_h };

    //Wall Shear Stress mediated degradation
    native_in >> K_tauw_d_e_h >> K_tauw_d_m_h >> K_tauw_d_ct_h >> K_tauw_d_cz_h >> K_tauw_d_cd1_h >> K_tauw_d_cd2_h;
    K_tauw_d_alpha_h = { K_tauw_d_e_h, K_tauw_d_m_h, K_tauw_d_ct_h,
                                        K_tauw_d_cz_h, K_tauw_d_cd1_h, K_tauw_d_cd2_h };

    //Initialize homeostatic loading variables
    //Pressure and Flow
    native_in >> P_h >> Q_h;
    //Set WSS reference in dynes/cm2
    mu = 0.04; //get_app_visc(this, sn);
    bar_tauw_h = 4*mu*Q_h/(3.14159265*pow(a_h*100, 3));

    //Initial stresses in each direction
    sigma_h.resize(9);
    sigma.resize(9);
    Cbar.resize(3);
    CC.resize(36);

    //Active stress parameters
    alpha_active = { 0, 1, 0, 0, 0, 0 };
    a_act.resize(nts);
    a_act[0] = a_h;
    lambda_act.resize(nts);
    lambda_act[0] = 1.0;
    native_in >> k_act >> lambda_0 >> lambda_m; //active remodelling time, min active stretch, max active stretch
    native_in >> CB >> CS; //vasodilator ratios
    native_in >> T_act_h; //homeostatic active stress magnitude
    T_act = T_act_h; //current active stress 

    //Initializion with looping through constituents
    double rhoR_0 = 0.0;
    double eta = 0;
    double g_alpha = 0;
    G_alpha_h.resize(3 * n_alpha);
    
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Initialize the deposition tensor
        eta = eta_alpha_h[alpha];
        g_alpha = g_alpha_h[alpha];
        if (eta >= 0) { //for anisotropic constituents
            G_alpha_h[3 * alpha] = 0.0;
            G_alpha_h[3 * alpha + 1] = g_alpha * sin(eta);
            G_alpha_h[3 * alpha + 2] = g_alpha * cos(eta);
        }
        else { //for isotropic constituents (i.e. elastin)
            G_alpha_h[3 * alpha] = G_e_rh;
            G_alpha_h[3 * alpha + 1] = G_e_th;
            G_alpha_h[3 * alpha + 2] = G_e_zh;
        }

        //Initialize homeostatic mass productions
        mR_alpha_h[alpha] = k_alpha_h[alpha] * rhoR_alpha_h[alpha];

        //Initialize simulation mass kinetics to hemeostatic values
        mR_alpha[nts * alpha] = mR_alpha_h[alpha];
        k_alpha[nts * alpha] = k_alpha_h[alpha];

        //Initialize volume fractions, which are equal to mass fractions if all constituents have same true mass density
        epsilonR_alpha[nts * alpha] = phi_alpha_h[alpha];
        epsilon_alpha[nts * alpha] = phi_alpha_h[alpha];

        //Initialize density time histories to their homeostatic values for native vessel
        rhoR_alpha[nts * alpha] = rhoR_alpha_h[alpha];

        //Initilize stretch histories
        lambda_alpha_tau[nts * alpha] = 1.0;
    }

    //Solve for native stress state
    //Initialize variables for load history
    P = P_h;
    bar_tauw = bar_tauw_h;
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0; 
    lambda_z_tau.resize(nts);
    lambda_z_tau[0] = 1.0;

    //NOTE: If input homeostatic configuration is not close to actual homeostatic configuration
    //this will give an erroneous 1st passive config. 
    //Solve for native passive configuration
    int equil_check = 0;
    int num_act_incr = 10;
    num_exp_flag = 1;
    T_act = T_act_h;
    num_exp_flag = 0;

    //Assign passive geometry parameters
    h_pas[0] = h[0];
    a_pas[0] = a[0];
    a_mid_pas[0] = a_mid[0];
    //Reinitialize to input geometry parameters for solving

    //Solve for native active (homeostatic) configuration
    T_act = T_act_h;
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;
    num_exp_flag = 0;

    native_in >> sigma_inv_h;
    sigma_inv = sigma_inv_h;
    mu = get_app_visc(this, sn);
    bar_tauw_h = 4 * mu * Q_h / (3.14159265 * pow(a_h * 100, 3));

    native_in >> P;
    bar_tauw = 4 * mu * Q_h / (3.14159265 * pow(a[0] * 100, 3));;
    native_in >> Q;

    //For pressure ramping
    P_prev = P;
    native_in.close();

}




//Initialize the reference vessel for the simulation    
void vessel::initializeNativeHandshake(string native_name, double n_days_inp, double dt_inp) {

    //Set native vessel time parameters
    double n_days = n_days_inp;
    dt = dt_inp;
    nts = int(n_days / dt);
    sn = 0;
    s = 0.0;

    double mu = 0; //apparent viscosity

    //Load the expereimentally determined and prescribed properties of the vessel from file
    //Input arguments for scaffold input file (ELS)
    std::ifstream native_in(native_name);

    native_in >> vessel_name; //Type of vessel simulated

    //Initilize the parameters for the reference native vessel
    //Geometric parameters
    //Homeostatic parameters are those for a NATIVE vessel
    native_in >> a_h; //in vivo reference inner radius
    a_h = a_h * mm_to_m;
    native_in >> h_h; //in vivo reference medial thickness
    h_h = h_h * mm_to_m;

    //Axial stretch with in vivo reference
    native_in >> lambda_z_h; //Reference in vivo stretch    

    //Constituent material properties
    n_alpha = 6; //number of constituents alpha
    n_pol_alpha = 0; //number of polymeric constituents
    alpha_infl = { 0, 0, 0, 0, 0, 0 }; //flags for inflammation

    //Medial and adventitial constituents have the same material behavior
    native_in >> c1_e >> c2_e >> c1_m >> c2_m >> c1_ct >> c2_ct;
    native_in >> c1_cz >> c2_cz >> c1_cd1 >> c2_cd1 >> c1_cd2 >> c2_cd2;
    c_alpha_h = { c1_e, c2_e, c1_m, c2_m, c1_ct, c2_ct,
                                c1_cz, c2_cz, c1_cd1, c2_cd1, c1_cd2, c2_cd2 };

    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents
    native_in >> eta_e_h >> eta_m_h >> eta_ct_h >> eta_cz_h >> eta_cd1_h >> eta_cd2_h;
    eta_alpha_h = { eta_e_h * M_PI / 180.0, eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0,
                                  eta_cz_h * M_PI / 180.0, eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0 };

    //Pre-stretch parameters
    native_in >> G_e_rh >> G_e_th >> G_e_zh;
    native_in >> g_e_h >> g_m_h >> g_ct_h >> g_cz_h >> g_cd1_h >> g_cd2_h;
    g_alpha_h = { g_e_h, g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h };

    //Mass density parameters
    native_in >> rho_hat_h; //density of tissue as a whole kg/m^3
    rhoR_h = rho_hat_h;
    rho[0] = rho_hat_h;
    rhoR[0] = rho_hat_h;
    rho_hat_alpha_h = { rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h };

    //Homeostatic mass fractions of constituents
    native_in >> phi_e_h >> phi_m_h >> phi_ct_h >> phi_cz_h >> phi_cd1_h >> phi_cd2_h;
    phi_alpha_h = { phi_e_h, phi_m_h, phi_ct_h, phi_cz_h, phi_cd1_h, phi_cd2_h };

    //Homeostatic mass densities
    rhoR_e_h = phi_e_h * rho_hat_h, rhoR_m_h = phi_m_h * rho_hat_h,
        rhoR_ct_h = phi_ct_h * rho_hat_h, rhoR_cz_h = phi_cz_h * rho_hat_h,
        rhoR_cd1_h = phi_cd1_h * rho_hat_h, rhoR_cd2_h = phi_cd2_h * rho_hat_h;
    rhoR_alpha_h = { rhoR_e_h, rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h };

    //Kinetic parameters
    //Degradation parameters
    native_in >> k_e_h >> k_m_h >> k_ct_h >> k_cz_h >> k_cd1_h >> k_cd2_h;
    k_alpha_h = { k_e_h, k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h };

    //Initialize kinetic arrays
    rho.resize(nts); //Current mass density history
    rhoR.resize(nts); //Referential mass density history
    rhoR_alpha.resize(nts * n_alpha); //referential apparent mass densities (time history)
    epsilon_alpha.resize(nts * n_alpha); //current volume fractions
    epsilonR_alpha.resize(nts * n_alpha); //referential volume fractions
    mR_alpha.resize(nts * n_alpha); //referential mass production rate (time history)
    k_alpha.resize(nts * n_alpha); //mass removal decay (time history)
    mR_alpha_h.resize(n_alpha);
    rhoR_alpha_h.resize(n_alpha);
    lambda_alpha_tau.resize(nts * n_alpha);

    //Gain parameters
    //Stress mediated production
    native_in >> K_sigma_p_e_h >> K_sigma_p_m_h >> K_sigma_p_ct_h >> K_sigma_p_cz_h >> K_sigma_p_cd1_h >> K_sigma_p_cd2_h;
    K_sigma_p_alpha_h = { K_sigma_p_e_h, K_sigma_p_m_h, K_sigma_p_ct_h,
                                        K_sigma_p_cz_h, K_sigma_p_cd1_h, K_sigma_p_cd2_h };

    //Stress mediated degradation
    native_in >> K_sigma_d_e_h >> K_sigma_d_m_h >> K_sigma_d_ct_h >> K_sigma_d_cz_h >> K_sigma_d_cd1_h >> K_sigma_d_cd2_h;
    K_sigma_d_alpha_h = { K_sigma_d_e_h, K_sigma_d_m_h, K_sigma_d_ct_h,
                                        K_sigma_d_cz_h, K_sigma_d_cd1_h, K_sigma_d_cd2_h };

    //Wall Shear Stress mediated production
    native_in >> K_tauw_p_e_h >> K_tauw_p_m_h >> K_tauw_p_ct_h >> K_tauw_p_cz_h >> K_tauw_p_cd1_h >> K_tauw_p_cd2_h;
    K_tauw_p_alpha_h = { K_tauw_p_e_h, K_tauw_p_m_h, K_tauw_p_ct_h,
                                        K_tauw_p_cz_h, K_tauw_p_cd1_h, K_tauw_p_cd2_h };

    //Wall Shear Stress mediated degradation
    native_in >> K_tauw_d_e_h >> K_tauw_d_m_h >> K_tauw_d_ct_h >> K_tauw_d_cz_h >> K_tauw_d_cd1_h >> K_tauw_d_cd2_h;
    K_tauw_d_alpha_h = { K_tauw_d_e_h, K_tauw_d_m_h, K_tauw_d_ct_h,
                                        K_tauw_d_cz_h, K_tauw_d_cd1_h, K_tauw_d_cd2_h };

    //Initialize homeostatic loading variables
    //Stress invariant and WSS
    native_in >> sigma_inv_h >> bar_tauw_h >> P_h;
    sigma_inv = sigma_inv_h;
    bar_tauw = bar_tauw_h;
    P = P_h;

    //Material stiffness matrices
    Cbar.resize(3);
    CC.resize(36);
    sigma.resize(9);

    //Active stress parameters
    alpha_active = { 0, 1, 0, 0, 0, 0 };
    a_act.resize(nts);
    a_act[0] = a_h;
    lambda_act.resize(nts);
    lambda_act[0] = 1.0;
    native_in >> k_act >> lambda_0 >> lambda_m; //active remodelling time, min active stretch, max active stretch
    native_in >> CB >> CS; //vasodilator ratios
    native_in >> T_act_h; //homeostatic active stress magnitude
    T_act = T_act_h; //current active stress 

    F_s.resize(9*nts);
    Fi_curr.resize(9);
    F_curr.resize(9);
    //Initialize to no stretch
    F_s[0] = 1.0;
    F_s[4] = 1.0;
    F_s[8] = 1.0;
    F_curr[0] = 1.0;
    F_curr[4] = 1.0;
    F_curr[8] = 1.0;
    Fi_curr[0] = 1.0;
    Fi_curr[4] = 1.0;
    Fi_curr[8] = 1.0;
    J_di = 1.0;

    //Initializion with looping through constituents
    double rhoR_0 = 0.0;
    double eta = 0;
    double g_alpha = 0;
    G_alpha_h.resize(3 * n_alpha);
    
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Initialize the deposition tensor
        eta = eta_alpha_h[alpha];
        g_alpha = g_alpha_h[alpha];
        if (eta >= 0) { //for anisotropic constituents
            G_alpha_h[3 * alpha] = 0.0;
            G_alpha_h[3 * alpha + 1] = g_alpha * sin(eta);
            G_alpha_h[3 * alpha + 2] = g_alpha * cos(eta);
        }
        else { //for isotropic constituents (i.e. elastin)
            G_alpha_h[3 * alpha] = G_e_rh;
            G_alpha_h[3 * alpha + 1] = G_e_th;
            G_alpha_h[3 * alpha + 2] = G_e_zh;
        }

        //Initialize homeostatic mass productions
        mR_alpha_h[alpha] = k_alpha_h[alpha] * rhoR_alpha_h[alpha];

        //Initialize simulation mass kinetics to hemeostatic values
        mR_alpha[nts * alpha] = mR_alpha_h[alpha];
        k_alpha[nts * alpha] = k_alpha_h[alpha];

        //Initialize volume fractions, which are equal to mass fractions if all constituents have same true mass density
        epsilonR_alpha[nts * alpha] = phi_alpha_h[alpha];
        epsilon_alpha[nts * alpha] = phi_alpha_h[alpha];

        //Initialize density time histories to their homeostatic values for native vessel
        rhoR_alpha[nts * alpha] = rhoR_alpha_h[alpha];

        //Initilize stretch histories
        lambda_alpha_tau[nts * alpha] = 1.0;
    }



    //Update stress to in vivo stress for homeostatic stress calc
    update_sigma_handshake(this);
    GnR_out.precision(12);
    HnS_out.precision(12);

    native_in.close();

}


void vessel::initializeTEVG(string scaffold_name, string immune_name, vessel const &native_vessel, double infl_scale_trans, double n_days_inp, double dt_inp) {
    //Copy initialization variables over from native vessel counterpart
    //Initialization parameters
    //Initializing constituents
    c1_e = native_vessel.c1_e;
    c2_e = native_vessel.c2_e;
    c1_m = native_vessel.c1_m;
    c2_m = native_vessel.c2_m;
    c1_ct = native_vessel.c1_ct;
    c2_ct = native_vessel.c2_ct;
    c1_cz = native_vessel.c1_cz;
    c2_cz = native_vessel.c2_cz;
    c1_cd1 = native_vessel.c1_cd1;
    c2_cd1 = native_vessel.c2_cd1;
    c1_cd2 = native_vessel.c1_cd2;
    c2_cd2 = native_vessel.c2_cd2;
    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents
    eta_e_h = native_vessel.eta_e_h;
    eta_m_h = native_vessel.eta_m_h;
    eta_ct_h = native_vessel.eta_ct_h;
    eta_cz_h = native_vessel.eta_cz_h;
    eta_cd1_h = native_vessel.eta_cd1_h;
    eta_cd2_h = native_vessel.eta_cd2_h;
    //Pre-stretch parameters
    g_e_h = native_vessel.g_e_h;
    g_m_h = native_vessel.g_m_h;
    g_ct_h = native_vessel.g_ct_h;
    g_cz_h = native_vessel.g_cz_h;
    g_cd1_h = native_vessel.g_cd1_h;
    g_cd2_h = native_vessel.g_cd2_h;
    //Mass density parameters
    rho_hat_h = native_vessel.rho_hat_h;
    //Homeostatic mass fractions of constituents
    phi_e_h = native_vessel.phi_e_h;
    phi_m_h = native_vessel.phi_m_h;
    phi_ct_h = native_vessel.phi_ct_h;
    phi_cz_h = native_vessel.phi_cz_h;
    phi_cd1_h = native_vessel.phi_cd1_h;
    phi_cd2_h = native_vessel.phi_cd2_h;
    //Homeostatic mass densities
    rhoR_e_h = native_vessel.rhoR_e_h;
    rhoR_m_h = native_vessel.rhoR_m_h;
    rhoR_ct_h = native_vessel.rhoR_ct_h;
    rhoR_cz_h = native_vessel.rhoR_cz_h;
    rhoR_cd1_h = native_vessel.rhoR_cd1_h;
    rhoR_cd2_h = native_vessel.rhoR_cd2_h;
    //Degradation parameters
    k_e_h = native_vessel.k_e_h;
    k_m_h = native_vessel.k_m_h;
    k_ct_h = native_vessel.k_ct_h;
    k_cz_h = native_vessel.k_cz_h;
    k_cd1_h = native_vessel.k_cd1_h;
    k_cd2_h = native_vessel.k_cd2_h;
    //Stress mediated production
    K_sigma_p_e_h = native_vessel.K_sigma_p_e_h;
    K_sigma_p_m_h = native_vessel.K_sigma_p_m_h;
    K_sigma_p_ct_h = native_vessel.K_sigma_p_ct_h;
    K_sigma_p_cz_h = native_vessel.K_sigma_p_cz_h;
    K_sigma_p_cd1_h = native_vessel.K_sigma_p_cd1_h;
    K_sigma_p_cd2_h = native_vessel.K_sigma_p_cd2_h;
    //Stress mediated degradation
    K_sigma_d_e_h = native_vessel.K_sigma_d_e_h;
    K_sigma_d_m_h = native_vessel.K_sigma_d_m_h;
    K_sigma_d_ct_h = native_vessel.K_sigma_d_ct_h;
    K_sigma_d_cz_h = native_vessel.K_sigma_d_cz_h;
    K_sigma_d_cd1_h = native_vessel.K_sigma_d_cd1_h;
    K_sigma_d_cd2_h = native_vessel.K_sigma_d_cd2_h;
    //Wall Shear Stress mediated production
    K_tauw_p_e_h = native_vessel.K_tauw_p_e_h;
    K_tauw_p_m_h = native_vessel.K_tauw_p_m_h;
    K_tauw_p_ct_h = native_vessel.K_tauw_p_ct_h;
    K_tauw_p_cz_h = native_vessel.K_tauw_p_cz_h;
    K_tauw_p_cd1_h = native_vessel.K_tauw_p_cd1_h;
    K_tauw_p_cd2_h = native_vessel.K_tauw_p_cd2_h;
    //Wall Shear Stress mediated degradation
    K_tauw_d_e_h = native_vessel.K_tauw_d_e_h;
    K_tauw_d_m_h = native_vessel.K_tauw_d_m_h;
    K_tauw_d_ct_h = native_vessel.K_tauw_d_ct_h;
    K_tauw_d_cz_h = native_vessel.K_tauw_d_cz_h;
    K_tauw_d_cd1_h = native_vessel.K_tauw_d_cd1_h;
    K_tauw_d_cd2_h = native_vessel.K_tauw_d_cd2_h;
    //Active Stress Parameters
    k_act = native_vessel.k_act;
    lambda_0 = native_vessel.lambda_0;
    lambda_m = native_vessel.lambda_m;
     //active remodelling time, min active stretch, max active stretch
    CB = native_vessel.CB; //vasodilator ratios
    CS = native_vessel.CS; //vasodilator ratios
    T_act_h = native_vessel.T_act_h; //homeostatic active stress magnitude
    T_act = native_vessel.T_act; //current active stress 

    //Read in scaffold properties/immune properties and store
    //Input arguments for scaffold input file (ELS)
    std::ifstream Scaffold_in(scaffold_name);
    std::ifstream Immune_in(immune_name);
    //std::ifstream Fit_in("FitParams_In.txt");

    Scaffold_in >> vessel_name;

    //Time parameters
    double n_days = n_days_inp; //days simulated
    dt = dt_inp; //time step size
    nts = int(n_days / dt); //number of G&R time steps
    sn = 0; //Initialize current time index to zero;
    s = 0; //Initialize the actual current time to zero

    double mu = 0; //apparent viscosity

    //Geometric parameters
    //Homeostatic parameters are those for a NATIVE vessel
    Scaffold_in >> A_h; //unloaded inner radius
    A_h = A_h * mm_to_m;
    Scaffold_in >> H_h; //unloaded medial thickness
    H_h = H_h * mm_to_m;
    A_mid_h = A_h + H_h / 2;

    //Initialize the loaded geometry history
    a.resize(nts); //loaded inner radius history
    a[0] = A_h;
    a_h = A_h;
    //curr_vessel.b.resize(curr_vessel.nts); //loaded outer radius history
    //curr_vessel.b[0] = curr_vessel.b_h;
    h.resize(nts); //loaded thickness history
    h[0] = H_h;
    h_h = h[0];
    a_mid.resize(nts); //loaded mid-radial history
    a_mid[0] = A_mid_h;
    a_mid_h = a_mid[0];

    //Axial stretch with in vivo reference
    lambda_z_h = 1.0; //Reference in vivo stretch    

    //Constituent material properties
    n_alpha = 13; //number of constituents alpha
    n_pol_alpha = 2; //number of polmyer constituents
    n_native_alpha = 5; //number of constituents from native vessel
    alpha_infl = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0 };

    //PGA, PCLA, mech SMC, 4 mech col ff, infl SMC, 4 infl col ff 
    double c1_p1, c2_p1, c1_p2, c2_p2, c1_gnd, c2_gnd;
    c2_gnd = 0.0;
    c1_gnd = 10;
    Scaffold_in >> c1_p1 >> c2_p1 >> c1_p2 >> c2_p2;
    //Fit_in >> c1_p1 >> c1_p2;;
    double gamma_i_1, gamma_i_2;
    Immune_in >> gamma_i_1 >> gamma_i_2;
    c_alpha_h = { c1_p1, c2_p1, c1_p2, c2_p2,
        c1_m, c2_m, c1_ct, c2_ct, c1_cz, c2_cz, c1_cd1, c2_cd1, c1_cd2, c2_cd2,
        c1_m * gamma_i_1, c2_m * gamma_i_2, c1_ct * gamma_i_1, c2_ct * gamma_i_2,
        c1_cz * gamma_i_1, c2_cz * gamma_i_2, c1_cd1 * gamma_i_1, c2_cd1 * gamma_i_2,
        c1_cd2 * gamma_i_1, c2_cd2 * gamma_i_2, c1_gnd, c2_gnd };

    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents, for polymer and elastin
    double eta_p1_h, eta_p2_h, eta_gnd_h;
    Scaffold_in >> eta_p1_h >> eta_p2_h;
    eta_gnd_h = eta_e_h;
    eta_alpha_h = { eta_p1_h * M_PI / 180.0, eta_p2_h * M_PI / 180.0,
                         eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0, eta_cz_h * M_PI / 180.0,
                         eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0,
                         eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0, eta_cz_h * M_PI / 180.0,
                         eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0, eta_gnd_h * M_PI / 180.0 };

    //Pre-stretch parameters
    double g_p1_h, g_p2_h, g_gnd_h;
    Scaffold_in >> g_p1_h >> g_p2_h;
    g_gnd_h = 1.0;
    g_alpha_h = { g_p1_h, g_p2_h,
                       g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h,
                       g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h,
                       g_gnd_h };

    //Mass density parameters
    //True mass densities
    double rho_hat_p1, rho_hat_p2;
    Scaffold_in >> rho_hat_p1 >> rho_hat_p2;
    rho_hat_alpha_h = { rho_hat_p1, rho_hat_p2, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h,
                            rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h };

    //Volume fractions only for polymer constituents
    double epsilon_p1_0, epsilon_p2_0, epsilon_gnd_0;
    Scaffold_in >> epsilon_p1_0 >> epsilon_p2_0;
    epsilon_gnd_0 = 1.0 - epsilon_p1_0 - epsilon_p2_0;
    epsilonR_alpha_0 = { epsilon_p1_0, epsilon_p2_0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, epsilon_gnd_0 };

    //Minimum volume fractions for damage-like behavior
    epsilon_pol_min = epsilonR_alpha_0;

    //Initialize the traction free geometry history
    A.resize(nts); //TF inner radius history
    //curr_vessel.B.resize(curr_vessel.nts); //TF outer radius history
    A_mid.resize(nts); //TF mid-radial history
    H.resize(nts); //TF Thickness history
    lambda_z_pre.resize(nts); //Axial pre-stretch

    //Kinetic parameters
    //Degradation parameters
    double k_gnd_h;
    k_gnd_h = k_m_h;
    k_alpha_h = { 0.0, 0.0, k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h,
                            k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h, 0.0 };

    //Gains for mechano-mediated kinetics
    K_sigma_p_alpha_h = { 0.0, 0.0,
                               K_sigma_p_m_h, K_sigma_p_ct_h, K_sigma_p_cz_h,
                               K_sigma_p_cd1_h, K_sigma_p_cd2_h,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_sigma_d_alpha_h = { 0.0, 0.0,
                              K_sigma_d_m_h, K_sigma_d_ct_h, K_sigma_d_cz_h,
                              K_sigma_d_cd1_h, K_sigma_d_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_tauw_p_alpha_h = { 0.0, 0.0,
                              K_tauw_p_m_h, K_tauw_p_ct_h, K_tauw_p_cz_h,
                              K_tauw_p_cd1_h, K_tauw_p_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_tauw_d_alpha_h = { 0.0, 0.0,
                              K_tauw_d_m_h, K_tauw_d_ct_h, K_tauw_d_cz_h,
                              K_tauw_d_cd1_h, K_tauw_d_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    //Initialize homeostatic apparent mass densities
    rhoR_alpha_h = { 0.0, 0.0, rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h,
                          rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h, 0.0 };

    rhoR_h = native_vessel.rhoR_h;

    //Initializing time history storage
    rho.resize(nts); //Current mass density history
    rhoR.resize(nts); //Referential mass density history
    rhoR_alpha.resize(nts * n_alpha); //referential apparent mass densities (time history)
    epsilon_alpha.resize(nts * n_alpha); //current volume fractions
    epsilonR_alpha.resize(nts * n_alpha); //referential volume fractions
    mR_alpha.resize(nts * n_alpha); //referential mass production rate (time history)
    k_alpha.resize(nts * n_alpha); //mass removal decay (time history)
    mR_alpha_h.resize(n_alpha);
    rhoR_alpha_h.resize(n_alpha);
    lambda_alpha_tau.resize(nts * n_alpha);

    //Initializion with looping through constituents
    double rhoR_0 = 0.0;
    double eta = 0;
    double g_alpha = 0;
    G_alpha_h.resize(3 * n_alpha);
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Initialize the deposition tensor
        eta = eta_alpha_h[alpha];
        g_alpha = g_alpha_h[alpha];
        if (eta >= 0) { //for anisotropic constituents
            G_alpha_h[3 * alpha] = 0.0;
            G_alpha_h[3 * alpha + 1] = g_alpha * sin(eta);
            G_alpha_h[3 * alpha + 2] = g_alpha * cos(eta);
        }
        else { //for isotropic constituents
            G_alpha_h[3 * alpha] = 1.0 / pow(g_alpha, 2);
            G_alpha_h[3 * alpha + 1] = g_alpha;
            G_alpha_h[3 * alpha + 2] = g_alpha;
        }

        //Initialize homeostatic mass productions
        mR_alpha_h[alpha] = k_alpha_h[alpha] * rhoR_alpha_h[alpha];

        //Initialize time histories to their homeostatic values for native vessel
        rhoR_alpha[nts * alpha] = epsilonR_alpha_0[alpha] * rho_hat_alpha_h[alpha];
        rhoR_0 += rhoR_alpha[nts * alpha];

        mR_alpha[nts * alpha] = 0; // mR_alpha_h[alpha];
        k_alpha[nts * alpha] = 0; // k_alpha_h[alpha];
        epsilonR_alpha[nts * alpha] = epsilonR_alpha_0[alpha];
        epsilon_alpha[nts * alpha] = epsilonR_alpha_0[alpha];

        //Initilize stretch histories
        lambda_alpha_tau[nts * alpha] = 1.0;

    }
    rhoR[0] = rhoR_0;
    rho[0] = rhoR_0;

    //Pre-calculating kinetic quanities for polymer degradation and immuno-driven production/removal
    double k_p1, zeta_p1, k_p2, zeta_p2, gamma_p_d1, gamma_p_d2;

    //In vitro degradation characteristics
    Scaffold_in >> k_p1 >> zeta_p1 >> k_p2 >> zeta_p2;
    //Fit_in >> k_p1 >> zeta_p1 >> k_p2 >> zeta_p2;

    //Enhanced degradation in vivo
    Immune_in >> gamma_p_d1;
    gamma_p_d2 = 1 / gamma_p_d1;

    //Scaffold properties affecting inflammation
    double fd_p_1, fd_p_2, fd_p, ps_p, epsilon_p;
    Scaffold_in >> fd_p_1 >> fd_p_2;

    //Determine bulk scaffold properties
    epsilon_p = epsilon_p1_0 + epsilon_p2_0;
    fd_p = (fd_p_1 < fd_p_2) * fd_p_1 + (fd_p_2 < fd_p_1) * fd_p_2 + (fd_p_1 == fd_p_2) * fd_p_1;
    ps_p = -sqrt(M_PI) / 4 * (1 + M_PI / (2 * log(1 - epsilon_p))) * fd_p;

    //Immunological production and degradation parameters
    double K_i_p_mic, K_i_p_wound, K_i_d_max, delta_i_p, beta_i_p;
    Immune_in >> K_i_p_mic >> K_i_p_wound >> K_i_d_max >> delta_i_p >> beta_i_p;

    double ps_norm, fd_norm;
    Immune_in >> ps_norm >> fd_norm;

    double K_i_p_trans = infl_scale_trans * K_i_p_mic * (ps_p / ps_norm * fd_p / fd_norm) + K_i_p_wound;
    double K_i_p_steady = K_i_p_mic * (ps_norm / ps_p + fd_p / fd_norm - 2);
    K_i_p_steady = (K_i_p_steady > 0)* K_i_p_steady;
    double K_i_d = K_i_p_trans / K_i_d_max;

    //Polymer degradation initialization
    double Q_p1, Q_p2, Q_gnd;
    double s_gnd_off = 14.0, epsilonR_gnd_min = 0.10, s = 0.0;

    //Inflammation initialization
    double gamma_fun, steady_fun;
    ups_infl_p.resize(nts);
    ups_infl_d.resize(nts);

    for (int sn = 1; sn < nts; sn++) {
        //Calculate polymer/ground degradation
        s = sn * dt;
        Q_p1 = (1 + exp(-k_p1 * zeta_p1)) / (1 + exp(k_p1 * gamma_p_d1 * (s - zeta_p1 * gamma_p_d2)));
        Q_p2 = (1 + exp(-k_p2 * zeta_p2)) / (1 + exp(k_p2 * gamma_p_d1 * (s - zeta_p2 * gamma_p_d2)));
        Q_gnd = (s > s_gnd_off)* ((1 - epsilonR_gnd_min) * exp(-k_gnd_h * (s - s_gnd_off)) + epsilonR_gnd_min)
            + (s <= s_gnd_off) * 1;

        epsilonR_alpha[0 * nts + sn] = Q_p1 * epsilonR_alpha_0[0];
        epsilonR_alpha[1 * nts + sn] = Q_p2 * epsilonR_alpha_0[1];
        epsilonR_alpha[12 * nts + sn] = Q_gnd * epsilonR_alpha_0[12];

        rhoR_alpha[0 * nts + sn] = epsilonR_alpha[0 * nts + sn] * rho_hat_alpha_h[0];
        rhoR_alpha[1 * nts + sn] = epsilonR_alpha[1 * nts + sn] * rho_hat_alpha_h[1];
        rhoR_alpha[12 * nts + sn] = epsilonR_alpha[12 * nts + sn] * rho_hat_alpha_h[12];

        //Calculate immunological stimulus
        gamma_fun = (pow(delta_i_p, beta_i_p)) * pow(s, beta_i_p - 1) * exp(-delta_i_p * s)
            / (delta_i_p * pow((beta_i_p - 1), (beta_i_p - 1)) * exp(1 - beta_i_p));
        steady_fun = (1 - exp(-delta_i_p * s));
        ups_infl_p[sn] = K_i_p_trans* gamma_fun + K_i_p_steady * steady_fun;
        ups_infl_d[sn] = K_i_d;

    }

    //Set initial loading conditions
    P =  native_vessel.P_h;
    Q_h = native_vessel.Q_h;
    Q = native_vessel.Q_h;
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;
    lambda_z_tau.resize(nts);
    lambda_z_tau[0] = 1.0;

    //Intialize stresses
    sigma_h.resize(3);
    sigma.resize(3);
    Cbar.resize(3);
    CC.resize(36);

    //Active stress parameters
    //NOTE: currently assume TEVG has no active stress contribution
    alpha_active = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    a_act.resize(nts);
    a_act[0] = a_h;
    lambda_act.resize(nts);
    lambda_act[0] = 1.0;

    //Find the initial loaded state
    update_sigma(this);
    num_exp_flag = 1;
    int equil_check = find_iv_geom(this);
    num_exp_flag = 0;

    //printf("%s %f %s %f\n", "Inner radius: ", a[0],
    //    "Thickness: ", h[0]);
    //fflush(stdout);

    //Need to add pre-stretch to the polymeric constituents to give them and in vivo reference
    vector<double> G_invivo_ref = { A_mid_h / a_mid[0], a_mid[0] / A_mid_h, 1.0 };
    for (int alpha = 0; alpha < n_pol_alpha + 1; alpha++) {

        //Account for ground matrix the last constituent
        if (alpha == n_pol_alpha) {
            G_alpha_h[3 * (n_alpha - 1)] = G_invivo_ref[0];
            G_alpha_h[3 * (n_alpha - 1) + 1] = G_invivo_ref[1];
            G_alpha_h[3 * (n_alpha - 1) + 2] = G_invivo_ref[2];
        }
        else {
            G_alpha_h[3 * alpha] = G_invivo_ref[0];
            G_alpha_h[3 * alpha + 1] = G_invivo_ref[1];
            G_alpha_h[3 * alpha + 2] = G_invivo_ref[2];
        }
    }

    //Reset the initial stretch conditions
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;

    //Reinitialize the reference loaded configuration 
    //to be the initial loaded configuration
    a_mid_h = a_mid[0];
    a_h = a[0];
    h_h = h[0];

    //Assign passive geometry parameters to be the same as active
    h_pas[0] = h[0];
    a_pas[0] = a[0];
    a_mid_pas[0] = a_mid[0];

    //Store homeostatic stresses from native vessel
    bar_tauw_h = native_vessel.bar_tauw_h;
    f_h = native_vessel.f_h;
    sigma_h = native_vessel.sigma_h;
    sigma_inv_h = native_vessel.sigma_inv_h;

    //Initialize current pressure and WSS
    mu = get_app_visc(this, sn);
    bar_tauw = 4*mu*Q_h/(3.14159265*pow(a[0]*100, 3));

    //To only consider initial constituents and prescribed mass changes
    pol_only_flag = 0;

    //For pressure ramping
    P_prev = P;

    Scaffold_in.close();
    Immune_in.close();
}

void vessel::initializeTEVGExplicit(string scaffold_name, string immune_name, vessel const &native_vessel, double infl_scale_trans, double n_days_inp, double dt_inp) {
    //Copy initialization variables over from native vessel counterpart
    //Initialization parameters
    //Initializing constituents
    c1_e = native_vessel.c1_e;
    c2_e = native_vessel.c2_e;
    c1_m = native_vessel.c1_m;
    c2_m = native_vessel.c2_m;
    c1_ct = native_vessel.c1_ct;
    c2_ct = native_vessel.c2_ct;
    c1_cz = native_vessel.c1_cz;
    c2_cz = native_vessel.c2_cz;
    c1_cd1 = native_vessel.c1_cd1;
    c2_cd1 = native_vessel.c2_cd1;
    c1_cd2 = native_vessel.c1_cd2;
    c2_cd2 = native_vessel.c2_cd2;
    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents
    eta_e_h = native_vessel.eta_e_h;
    eta_m_h = native_vessel.eta_m_h;
    eta_ct_h = native_vessel.eta_ct_h;
    eta_cz_h = native_vessel.eta_cz_h;
    eta_cd1_h = native_vessel.eta_cd1_h;
    eta_cd2_h = native_vessel.eta_cd2_h;
    //Pre-stretch parameters
    g_e_h = native_vessel.g_e_h;
    g_m_h = native_vessel.g_m_h;
    g_ct_h = native_vessel.g_ct_h;
    g_cz_h = native_vessel.g_cz_h;
    g_cd1_h = native_vessel.g_cd1_h;
    g_cd2_h = native_vessel.g_cd2_h;
    //Mass density parameters
    rho_hat_h = native_vessel.rho_hat_h;
    //Homeostatic mass fractions of constituents
    phi_e_h = native_vessel.phi_e_h;
    phi_m_h = native_vessel.phi_m_h;
    phi_ct_h = native_vessel.phi_ct_h;
    phi_cz_h = native_vessel.phi_cz_h;
    phi_cd1_h = native_vessel.phi_cd1_h;
    phi_cd2_h = native_vessel.phi_cd2_h;
    //Homeostatic mass densities
    rhoR_e_h = native_vessel.rhoR_e_h;
    rhoR_m_h = native_vessel.rhoR_m_h;
    rhoR_ct_h = native_vessel.rhoR_ct_h;
    rhoR_cz_h = native_vessel.rhoR_cz_h;
    rhoR_cd1_h = native_vessel.rhoR_cd1_h;
    rhoR_cd2_h = native_vessel.rhoR_cd2_h;
    //Degradation parameters
    k_e_h = native_vessel.k_e_h;
    k_m_h = native_vessel.k_m_h;
    k_ct_h = native_vessel.k_ct_h;
    k_cz_h = native_vessel.k_cz_h;
    k_cd1_h = native_vessel.k_cd1_h;
    k_cd2_h = native_vessel.k_cd2_h;
    //Stress mediated production
    K_sigma_p_e_h = native_vessel.K_sigma_p_e_h;
    K_sigma_p_m_h = native_vessel.K_sigma_p_m_h;
    K_sigma_p_ct_h = native_vessel.K_sigma_p_ct_h;
    K_sigma_p_cz_h = native_vessel.K_sigma_p_cz_h;
    K_sigma_p_cd1_h = native_vessel.K_sigma_p_cd1_h;
    K_sigma_p_cd2_h = native_vessel.K_sigma_p_cd2_h;
    //Stress mediated degradation
    K_sigma_d_e_h = native_vessel.K_sigma_d_e_h;
    K_sigma_d_m_h = native_vessel.K_sigma_d_m_h;
    K_sigma_d_ct_h = native_vessel.K_sigma_d_ct_h;
    K_sigma_d_cz_h = native_vessel.K_sigma_d_cz_h;
    K_sigma_d_cd1_h = native_vessel.K_sigma_d_cd1_h;
    K_sigma_d_cd2_h = native_vessel.K_sigma_d_cd2_h;
    //Wall Shear Stress mediated production
    K_tauw_p_e_h = native_vessel.K_tauw_p_e_h;
    K_tauw_p_m_h = native_vessel.K_tauw_p_m_h;
    K_tauw_p_ct_h = native_vessel.K_tauw_p_ct_h;
    K_tauw_p_cz_h = native_vessel.K_tauw_p_cz_h;
    K_tauw_p_cd1_h = native_vessel.K_tauw_p_cd1_h;
    K_tauw_p_cd2_h = native_vessel.K_tauw_p_cd2_h;
    //Wall Shear Stress mediated degradation
    K_tauw_d_e_h = native_vessel.K_tauw_d_e_h;
    K_tauw_d_m_h = native_vessel.K_tauw_d_m_h;
    K_tauw_d_ct_h = native_vessel.K_tauw_d_ct_h;
    K_tauw_d_cz_h = native_vessel.K_tauw_d_cz_h;
    K_tauw_d_cd1_h = native_vessel.K_tauw_d_cd1_h;
    K_tauw_d_cd2_h = native_vessel.K_tauw_d_cd2_h;
    //Active Stress Parameters
    k_act = native_vessel.k_act;
    lambda_0 = native_vessel.lambda_0;
    lambda_m = native_vessel.lambda_m;
     //active remodelling time, min active stretch, max active stretch
    CB = native_vessel.CB; //vasodilator ratios
    CS = native_vessel.CS; //vasodilator ratios
    T_act_h = native_vessel.T_act_h; //homeostatic active stress magnitude
    T_act = native_vessel.T_act; //current active stress 

    //Read in scaffold properties/immune properties and store
    //Input arguments for scaffold input file (ELS)
    std::ifstream Scaffold_in(scaffold_name);
    std::ifstream Immune_in(immune_name);
    //std::ifstream Fit_in("FitParams_In.txt");

    Scaffold_in >> vessel_name;

    //Time parameters
    double n_days = n_days_inp; //days simulated
    dt = dt_inp; //time step size
    nts = int(n_days / dt); //number of G&R time steps
    sn = 0; //Initialize current time index to zero;
    s = 0; //Initialize the actual current time to zero

    double mu = 0; //apparent viscosity

    //Geometric parameters
    //Homeostatic parameters are those for a NATIVE vessel
    Scaffold_in >> A_h; //unloaded inner radius
    A_h = A_h * mm_to_m;
    Scaffold_in >> H_h; //unloaded medial thickness
    H_h = H_h * mm_to_m;
    A_mid_h = A_h + H_h / 2;

    //Initialize the loaded geometry history
    a.resize(nts); //loaded inner radius history
    a[0] = A_h;
    a_h = A_h;
    //curr_vessel.b.resize(curr_vessel.nts); //loaded outer radius history
    //curr_vessel.b[0] = curr_vessel.b_h;
    h.resize(nts); //loaded thickness history
    h[0] = H_h;
    h_h = h[0];
    a_mid.resize(nts); //loaded mid-radial history
    a_mid[0] = A_mid_h;
    a_mid_h = a_mid[0];

    //Axial stretch with in vivo reference
    lambda_z_h = 1.0; //Reference in vivo stretch    

    //Constituent material properties
    n_alpha = 13; //number of constituents alpha
    n_pol_alpha = 2; //number of polmyer constituents
    n_native_alpha = 5; //number of constituents from native vessel
    alpha_infl = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0 };

    //PGA, PCLA, mech SMC, 4 mech col ff, infl SMC, 4 infl col ff 
    double c1_p1, c2_p1, c1_p2, c2_p2, c1_gnd, c2_gnd;
    c2_gnd = 0.0;
    c1_gnd = 10;
    Scaffold_in >> c1_p1 >> c2_p1 >> c1_p2 >> c2_p2;
    //Fit_in >> c1_p1 >> c1_p2;;
    double gamma_i_1, gamma_i_2;
    Immune_in >> gamma_i_1 >> gamma_i_2;
    c_alpha_h = { c1_p1, c2_p1, c1_p2, c2_p2,
        c1_m, c2_m, c1_ct, c2_ct, c1_cz, c2_cz, c1_cd1, c2_cd1, c1_cd2, c2_cd2,
        c1_m * gamma_i_1, c2_m * gamma_i_2, c1_ct * gamma_i_1, c2_ct * gamma_i_2,
        c1_cz * gamma_i_1, c2_cz * gamma_i_2, c1_cd1 * gamma_i_1, c2_cd1 * gamma_i_2,
        c1_cd2 * gamma_i_1, c2_cd2 * gamma_i_2, c1_gnd, c2_gnd };

    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents, for polymer and elastin
    double eta_p1_h, eta_p2_h, eta_gnd_h;
    Scaffold_in >> eta_p1_h >> eta_p2_h;
    eta_gnd_h = eta_e_h;
    eta_alpha_h = { eta_p1_h * M_PI / 180.0, eta_p2_h * M_PI / 180.0,
                         eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0, eta_cz_h * M_PI / 180.0,
                         eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0,
                         eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0, eta_cz_h * M_PI / 180.0,
                         eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0, eta_gnd_h * M_PI / 180.0 };

    //Pre-stretch parameters
    double g_p1_h, g_p2_h, g_gnd_h;
    Scaffold_in >> g_p1_h >> g_p2_h;
    g_gnd_h = 1.0;
    g_alpha_h = { g_p1_h, g_p2_h,
                       g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h,
                       g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h,
                       g_gnd_h };

    //Mass density parameters
    //True mass densities
    double rho_hat_p1, rho_hat_p2;
    Scaffold_in >> rho_hat_p1 >> rho_hat_p2;
    rho_hat_alpha_h = { rho_hat_p1, rho_hat_p2, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h,
                            rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h };

    //Volume fractions only for polymer constituents
    double epsilon_p1_0, epsilon_p2_0, epsilon_gnd_0;
    Scaffold_in >> epsilon_p1_0 >> epsilon_p2_0;
    epsilon_gnd_0 = 1.0 - epsilon_p1_0 - epsilon_p2_0;
    epsilonR_alpha_0 = { epsilon_p1_0, epsilon_p2_0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, epsilon_gnd_0 };

    //Minimum volume fractions for damage-like behavior
    epsilon_pol_min = epsilonR_alpha_0;

    //Initialize the traction free geometry history
    A.resize(nts); //TF inner radius history
    //curr_vessel.B.resize(curr_vessel.nts); //TF outer radius history
    A_mid.resize(nts); //TF mid-radial history
    H.resize(nts); //TF Thickness history
    lambda_z_pre.resize(nts); //Axial pre-stretch

    //Kinetic parameters
    //Degradation parameters
    double k_gnd_h;
    k_gnd_h = k_m_h;
    k_alpha_h = { 0.0, 0.0, k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h,
                       k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h, 0.0 };

    //Gains for mechano-mediated kinetics
    K_sigma_p_alpha_h = { 0.0, 0.0,
                               K_sigma_p_m_h, K_sigma_p_ct_h, K_sigma_p_cz_h,
                               K_sigma_p_cd1_h, K_sigma_p_cd2_h,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_sigma_d_alpha_h = { 0.0, 0.0,
                              K_sigma_d_m_h, K_sigma_d_ct_h, K_sigma_d_cz_h,
                              K_sigma_d_cd1_h, K_sigma_d_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_tauw_p_alpha_h = { 0.0, 0.0,
                              K_tauw_p_m_h, K_tauw_p_ct_h, K_tauw_p_cz_h,
                              K_tauw_p_cd1_h, K_tauw_p_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_tauw_d_alpha_h = { 0.0, 0.0,
                              K_tauw_d_m_h, K_tauw_d_ct_h, K_tauw_d_cz_h,
                              K_tauw_d_cd1_h, K_tauw_d_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    //Initialize homeostatic apparent mass densities
    rhoR_alpha_h = { 0.0, 0.0, rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h,
                          rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h, 0.0 };

    rhoR_h = native_vessel.rhoR_h;

    //Initializing time history storage
    rho.resize(nts); //Current mass density history
    rhoR.resize(nts); //Referential mass density history
    rhoR_alpha.resize(nts * n_alpha); //referential apparent mass densities (time history)
    epsilon_alpha.resize(nts * n_alpha); //current volume fractions
    epsilonR_alpha.resize(nts * n_alpha); //referential volume fractions
    mR_alpha.resize(nts * n_alpha); //referential mass production rate (time history)
    k_alpha.resize(nts * n_alpha); //mass removal decay (time history)
    mR_alpha_h.resize(n_alpha);
    rhoR_alpha_h.resize(n_alpha);
    lambda_alpha_tau.resize(nts * n_alpha);

    //Initializion with looping through constituents
    double rhoR_0 = 0.0;
    double eta = 0;
    double g_alpha = 0;
    G_alpha_h.resize(3 * n_alpha);
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Initialize the deposition tensor
        eta = eta_alpha_h[alpha];
        g_alpha = g_alpha_h[alpha];
        if (eta >= 0) { //for anisotropic constituents
            G_alpha_h[3 * alpha] = 0.0;
            G_alpha_h[3 * alpha + 1] = g_alpha * sin(eta);
            G_alpha_h[3 * alpha + 2] = g_alpha * cos(eta);
        }
        else { //for isotropic constituents
            G_alpha_h[3 * alpha] = 1.0 / pow(g_alpha, 2);
            G_alpha_h[3 * alpha + 1] = g_alpha;
            G_alpha_h[3 * alpha + 2] = g_alpha;
        }

        //Initialize homeostatic mass productions
        mR_alpha_h[alpha] = k_alpha_h[alpha] * rhoR_alpha_h[alpha];

        //Initialize time histories to their homeostatic values for native vessel
        rhoR_alpha[nts * alpha] = epsilonR_alpha_0[alpha] * rho_hat_alpha_h[alpha];
        rhoR_0 += rhoR_alpha[nts * alpha];

        mR_alpha[nts * alpha] = 0; // mR_alpha_h[alpha];
        k_alpha[nts * alpha] = 0; // k_alpha_h[alpha];
        epsilonR_alpha[nts * alpha] = epsilonR_alpha_0[alpha];
        epsilon_alpha[nts * alpha] = epsilonR_alpha_0[alpha];

        //Initilize stretch histories
        lambda_alpha_tau[nts * alpha] = 1.0;

    }
    rhoR[0] = rhoR_0;
    rho[0] = rhoR_0;

    //Pre-calculating kinetic quanities for polymer degradation and immuno-driven production/removal
    double k_p1, zeta_p1, k_p2, zeta_p2, gamma_p_d1, gamma_p_d2;

    //In vitro degradation characteristics
    Scaffold_in >> k_p1 >> zeta_p1 >> k_p2 >> zeta_p2;
    //Fit_in >> k_p1 >> zeta_p1 >> k_p2 >> zeta_p2;

    //Enhanced degradation in vivo
    Immune_in >> gamma_p_d1;
    gamma_p_d2 = 1 / gamma_p_d1;

    //Scaffold properties affecting inflammation
    double fd_p_1, fd_p_2, fd_p, ps_p, epsilon_p;
    Scaffold_in >> fd_p_1 >> fd_p_2;

    //Determine bulk scaffold properties
    epsilon_p = epsilon_p1_0 + epsilon_p2_0;
    fd_p = (fd_p_1 < fd_p_2) * fd_p_1 + (fd_p_2 < fd_p_1) * fd_p_2 + (fd_p_1 == fd_p_2) * fd_p_1;
    ps_p = -sqrt(M_PI) / 4 * (1 + M_PI / (2 * log(1 - epsilon_p))) * fd_p;

    //Immunological production and degradation parameters
    double K_i_p_mic, K_i_p_wound, K_i_d_max, delta_i_p, beta_i_p;
    Immune_in >> K_i_p_mic >> K_i_p_wound >> K_i_d_max >> delta_i_p >> beta_i_p;

    double ps_norm, fd_norm;
    Immune_in >> ps_norm >> fd_norm;

    double K_i_p_trans = infl_scale_trans * K_i_p_mic * (ps_p / ps_norm * fd_p / fd_norm) + K_i_p_wound;
    double K_i_p_steady = K_i_p_mic * (ps_norm / ps_p + fd_p / fd_norm - 2);
    K_i_p_steady = (K_i_p_steady > 0)* K_i_p_steady;
    double K_i_d = K_i_p_trans / K_i_d_max;

    //Polymer degradation initialization
    double Q_p1, Q_p2, Q_gnd;
    double s_gnd_off = 14.0, epsilonR_gnd_min = 0.10, s = 0.0;

    //Inflammation initialization
    double gamma_fun, steady_fun;
    ups_infl_p.resize(nts);
    ups_infl_d.resize(nts);

    for (int sn = 1; sn < nts; sn++) {
        //Calculate polymer/ground degradation
        s = sn * dt;
        Q_p1 = (1 + exp(-k_p1 * zeta_p1)) / (1 + exp(k_p1 * gamma_p_d1 * (s - zeta_p1 * gamma_p_d2)));
        Q_p2 = (1 + exp(-k_p2 * zeta_p2)) / (1 + exp(k_p2 * gamma_p_d1 * (s - zeta_p2 * gamma_p_d2)));
        Q_gnd = (s > s_gnd_off)* ((1 - epsilonR_gnd_min) * exp(-k_gnd_h * (s - s_gnd_off)) + epsilonR_gnd_min)
            + (s <= s_gnd_off) * 1;

        epsilonR_alpha[0 * nts + sn] = Q_p1 * epsilonR_alpha_0[0];
        epsilonR_alpha[1 * nts + sn] = Q_p2 * epsilonR_alpha_0[1];
        epsilonR_alpha[12 * nts + sn] = Q_gnd * epsilonR_alpha_0[12];

        rhoR_alpha[0 * nts + sn] = epsilonR_alpha[0 * nts + sn] * rho_hat_alpha_h[0];
        rhoR_alpha[1 * nts + sn] = epsilonR_alpha[1 * nts + sn] * rho_hat_alpha_h[1];
        rhoR_alpha[12 * nts + sn] = epsilonR_alpha[12 * nts + sn] * rho_hat_alpha_h[12];

        //Calculate immunological stimulus
        gamma_fun = (pow(delta_i_p, beta_i_p)) * pow(s, beta_i_p - 1) * exp(-delta_i_p * s)
            / (delta_i_p * pow((beta_i_p - 1), (beta_i_p - 1)) * exp(1 - beta_i_p));
        steady_fun = (1 - exp(-delta_i_p * s));
        ups_infl_p[sn] = K_i_p_trans* gamma_fun + K_i_p_steady * steady_fun;
        ups_infl_d[sn] = K_i_d;

    }

    //Set initial loading conditions
    P =  native_vessel.P_h;
    Q_h = native_vessel.Q_h;
    Q = native_vessel.Q_h;
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;
    lambda_z_tau.resize(nts);
    lambda_z_tau[0] = 1.0;

    //Intialize stresses
    sigma_h.resize(3);
    sigma.resize(3);
    Cbar.resize(3);
    CC.resize(36);

    //Active stress parameters
    //NOTE: currently assume TEVG has no active stress contribution
    alpha_active = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    a_act.resize(nts);
    a_act[0] = a_h;
    lambda_act.resize(nts);
    lambda_act[0] = 1.0;

    //Find the initial loaded state
    update_sigma(this);
    //num_exp_flag = 1;
    //int equil_check = find_iv_geom(this);
    num_exp_flag = 0;

    //printf("%s %f %s %f\n", "Inner radius: ", a[0],
    //    "Thickness: ", h[0]);
    //fflush(stdout);

    //Need to add pre-stretch to the polymeric constituents to give them and in vivo reference
    vector<double> G_invivo_ref = { A_mid_h / a_mid[0], a_mid[0] / A_mid_h, 1.0 };
    for (int alpha = 0; alpha < n_pol_alpha + 1; alpha++) {

        //Account for ground matrix the last constituent
        if (alpha == n_pol_alpha) {
            G_alpha_h[3 * (n_alpha - 1)] = G_invivo_ref[0];
            G_alpha_h[3 * (n_alpha - 1) + 1] = G_invivo_ref[1];
            G_alpha_h[3 * (n_alpha - 1) + 2] = G_invivo_ref[2];
        }
        else {
            G_alpha_h[3 * alpha] = G_invivo_ref[0];
            G_alpha_h[3 * alpha + 1] = G_invivo_ref[1];
            G_alpha_h[3 * alpha + 2] = G_invivo_ref[2];
        }
    }

    //Reset the initial stretch conditions
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;

    //Reinitialize the reference loaded configuration 
    //to be the initial loaded configuration
    a_mid_h = a_mid[0];
    a_h = a[0];
    h_h = h[0];

    //Assign passive geometry parameters to be the same as active
    h_pas[0] = h[0];
    a_pas[0] = a[0];
    a_mid_pas[0] = a_mid[0];

    //Store homeostatic stresses from native vessel
    bar_tauw_h = native_vessel.bar_tauw_h;
    f_h = native_vessel.f_h;
    sigma_h = native_vessel.sigma_h;
    sigma_inv_h = native_vessel.sigma_inv_h;

    //Initialize current pressure and WSS
    mu = get_app_visc(this, sn);
    bar_tauw = 4*mu*Q_h/(3.14159265*pow(a[0]*100, 3));

    //To only consider initial constituents and prescribed mass changes
    pol_only_flag = 0;

    //For pressure ramping
    P_prev = P;

    Scaffold_in.close();
    Immune_in.close();
}

void vessel::initializeTEVGHandshake(string scaffold_name, string immune_name, double infl_scale_trans, double n_days_inp, double dt_inp) {

    //Read in scaffold properties/immune properties and store
    //Input arguments for scaffold input file (ELS)
    std::ifstream Scaffold_in(scaffold_name);
    std::ifstream Immune_in(immune_name);
    //std::ifstream Fit_in("FitParams_In.txt");

    Scaffold_in >> vessel_name;

    //Time parameters
    double n_days = n_days_inp; //days simulated
    dt = dt_inp; //time step size
    nts = int(n_days / dt); //number of G&R time steps
    sn = 0; //Initialize current time index to zero;
    s = 0; //Initialize the actual current time to zero

    double mu = 0; //apparent viscosity

    //Geometric parameters
    //Homeostatic parameters are those for a NATIVE vessel
    Scaffold_in >> A_h; //unloaded inner radius
    A_h = A_h * mm_to_m;
    Scaffold_in >> H_h; //unloaded medial thickness
    H_h = H_h * mm_to_m;
    A_mid_h = A_h + H_h / 2;

    //Initialize the loaded geometry history
    a.resize(nts); //loaded inner radius history
    a[0] = A_h;
    a_h = A_h;
    //curr_vessel.b.resize(curr_vessel.nts); //loaded outer radius history
    //curr_vessel.b[0] = curr_vessel.b_h;
    h.resize(nts); //loaded thickness history
    h[0] = H_h;
    h_h = h[0];
    a_mid.resize(nts); //loaded mid-radial history
    a_mid[0] = A_mid_h;
    a_mid_h = a_mid[0];

    //Axial stretch with in vivo reference
    lambda_z_h = 1.0; //Reference in vivo stretch    

    //Constituent material properties
    n_alpha = 13; //number of constituents alpha
    n_pol_alpha = 2; //number of polmyer constituents
    n_native_alpha = 5; //number of constituents from native vessel
    alpha_infl = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0 };

    //PGA, PCLA, mech SMC, 4 mech col ff, infl SMC, 4 infl col ff 
    double c1_p1, c2_p1, c1_p2, c2_p2, c1_gnd, c2_gnd;
    c2_gnd = 0.0;
    c1_gnd = 10;
    Scaffold_in >> c1_p1 >> c2_p1 >> c1_p2 >> c2_p2;
    //Fit_in >> c1_p1 >> c1_p2;;
    double gamma_i_1, gamma_i_2;
    Immune_in >> gamma_i_1 >> gamma_i_2;
    c_alpha_h = { c1_p1, c2_p1, c1_p2, c2_p2,
        c1_m, c2_m, c1_ct, c2_ct, c1_cz, c2_cz, c1_cd1, c2_cd1, c1_cd2, c2_cd2,
        c1_m * gamma_i_1, c2_m * gamma_i_2, c1_ct * gamma_i_1, c2_ct * gamma_i_2,
        c1_cz * gamma_i_1, c2_cz * gamma_i_2, c1_cd1 * gamma_i_1, c2_cd1 * gamma_i_2,
        c1_cd2 * gamma_i_1, c2_cd2 * gamma_i_2, c1_gnd, c2_gnd };

    //Constituent orientations
    //Orientations in the reference configuration (the in vivo state for the DTA)
    //orientation < 1 for isotropic constituents, for polymer and elastin
    double eta_p1_h, eta_p2_h, eta_gnd_h;
    Scaffold_in >> eta_p1_h >> eta_p2_h;
    eta_gnd_h = eta_e_h;
    eta_alpha_h = { eta_p1_h * M_PI / 180.0, eta_p2_h * M_PI / 180.0,
                         eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0, eta_cz_h * M_PI / 180.0,
                         eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0,
                         eta_m_h * M_PI / 180.0, eta_ct_h * M_PI / 180.0, eta_cz_h * M_PI / 180.0,
                         eta_cd1_h * M_PI / 180.0, eta_cd2_h * M_PI / 180.0, eta_gnd_h * M_PI / 180.0 };

    //Pre-stretch parameters
    double g_p1_h, g_p2_h, g_gnd_h;
    Scaffold_in >> g_p1_h >> g_p2_h;
    g_gnd_h = 1.0;
    g_alpha_h = { g_p1_h, g_p2_h,
                       g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h,
                       g_m_h, g_ct_h, g_cz_h, g_cd1_h, g_cd2_h,
                       g_gnd_h };

    //Mass density parameters
    //True mass densities
    double rho_hat_p1, rho_hat_p2;
    Scaffold_in >> rho_hat_p1 >> rho_hat_p2;
    rho_hat_alpha_h = { rho_hat_p1, rho_hat_p2, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h,
                            rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h, rho_hat_h };

    //Volume fractions only for polymer constituents
    double epsilon_p1_0, epsilon_p2_0, epsilon_gnd_0;
    Scaffold_in >> epsilon_p1_0 >> epsilon_p2_0;
    epsilon_gnd_0 = 1.0 - epsilon_p1_0 - epsilon_p2_0;
    epsilonR_alpha_0 = { epsilon_p1_0, epsilon_p2_0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, epsilon_gnd_0 };

    //Minimum volume fractions for damage-like behavior
    epsilon_pol_min = epsilonR_alpha_0;

    //Initialize the traction free geometry history
    A.resize(nts); //TF inner radius history
    //curr_vessel.B.resize(curr_vessel.nts); //TF outer radius history
    A_mid.resize(nts); //TF mid-radial history
    H.resize(nts); //TF Thickness history
    lambda_z_pre.resize(nts); //Axial pre-stretch

    //Kinetic parameters
    //Degradation parameters
    double k_gnd_h;
    k_gnd_h = k_m_h;
    k_alpha_h = { 0.0, 0.0, k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h,
                       k_m_h, k_ct_h, k_cz_h, k_cd1_h, k_cd2_h, 0.0 };

    //Gains for mechano-mediated kinetics
    K_sigma_p_alpha_h = { 0.0, 0.0,
                               K_sigma_p_m_h, K_sigma_p_ct_h, K_sigma_p_cz_h,
                               K_sigma_p_cd1_h, K_sigma_p_cd2_h,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_sigma_d_alpha_h = { 0.0, 0.0,
                              K_sigma_d_m_h, K_sigma_d_ct_h, K_sigma_d_cz_h,
                              K_sigma_d_cd1_h, K_sigma_d_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_tauw_p_alpha_h = { 0.0, 0.0,
                              K_tauw_p_m_h, K_tauw_p_ct_h, K_tauw_p_cz_h,
                              K_tauw_p_cd1_h, K_tauw_p_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    K_tauw_d_alpha_h = { 0.0, 0.0,
                              K_tauw_d_m_h, K_tauw_d_ct_h, K_tauw_d_cz_h,
                              K_tauw_d_cd1_h, K_tauw_d_cd2_h,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    //Initialize homeostatic apparent mass densities
    rhoR_alpha_h = { 0.0, 0.0, rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h,
                          rhoR_m_h, rhoR_ct_h, rhoR_cz_h, rhoR_cd1_h, rhoR_cd2_h, 0.0 };

    //Initializing time history storage
    rho.resize(nts); //Current mass density history
    rhoR.resize(nts); //Referential mass density history
    rhoR_alpha.resize(nts * n_alpha); //referential apparent mass densities (time history)
    epsilon_alpha.resize(nts * n_alpha); //current volume fractions
    epsilonR_alpha.resize(nts * n_alpha); //referential volume fractions
    mR_alpha.resize(nts * n_alpha); //referential mass production rate (time history)
    k_alpha.resize(nts * n_alpha); //mass removal decay (time history)
    mR_alpha_h.resize(n_alpha);
    rhoR_alpha_h.resize(n_alpha);
    lambda_alpha_tau.resize(nts * n_alpha);

    //Initializion with looping through constituents
    double rhoR_0 = 0.0;
    double eta = 0;
    double g_alpha = 0;
    G_alpha_h.resize(3 * n_alpha);
    for (int alpha = 0; alpha < n_alpha; alpha++) {

        //Initialize the deposition tensor
        eta = eta_alpha_h[alpha];
        g_alpha = g_alpha_h[alpha];
        if (eta >= 0) { //for anisotropic constituents
            G_alpha_h[3 * alpha] = 0.0;
            G_alpha_h[3 * alpha + 1] = g_alpha * sin(eta);
            G_alpha_h[3 * alpha + 2] = g_alpha * cos(eta);
        }
        else { //for isotropic constituents
            G_alpha_h[3 * alpha] = 1.0 / pow(g_alpha, 2);
            G_alpha_h[3 * alpha + 1] = g_alpha;
            G_alpha_h[3 * alpha + 2] = g_alpha;
        }

        //Initialize homeostatic mass productions
        mR_alpha_h[alpha] = k_alpha_h[alpha] * rhoR_alpha_h[alpha];

        //Initialize time histories to their homeostatic values for native vessel
        rhoR_alpha[nts * alpha] = epsilonR_alpha_0[alpha] * rho_hat_alpha_h[alpha];
        rhoR_0 += rhoR_alpha[nts * alpha];

        mR_alpha[nts * alpha] = 0; // mR_alpha_h[alpha];
        k_alpha[nts * alpha] = 0; // k_alpha_h[alpha];
        epsilonR_alpha[nts * alpha] = epsilonR_alpha_0[alpha];
        epsilon_alpha[nts * alpha] = epsilonR_alpha_0[alpha];

        //Initilize stretch histories
        lambda_alpha_tau[nts * alpha] = 1.0;

    }
    rhoR[0] = rhoR_0;
    rho[0] = rhoR_0;

    //Pre-calculating kinetic quanities for polymer degradation and immuno-driven production/removal
    double k_p1, zeta_p1, k_p2, zeta_p2, gamma_p_d1, gamma_p_d2;

    //In vitro degradation characteristics
    Scaffold_in >> k_p1 >> zeta_p1 >> k_p2 >> zeta_p2;
    //Fit_in >> k_p1 >> zeta_p1 >> k_p2 >> zeta_p2;

    //Enhanced degradation in vivo
    Immune_in >> gamma_p_d1;
    gamma_p_d2 = 1 / gamma_p_d1;

    //Scaffold properties affecting inflammation
    double fd_p_1, fd_p_2, fd_p, ps_p, epsilon_p;
    Scaffold_in >> fd_p_1 >> fd_p_2;

    //Determine bulk scaffold properties
    epsilon_p = epsilon_p1_0 + epsilon_p2_0;
    fd_p = (fd_p_1 < fd_p_2) * fd_p_1 + (fd_p_2 < fd_p_1) * fd_p_2 + (fd_p_1 == fd_p_2) * fd_p_1;
    ps_p = -sqrt(M_PI) / 4 * (1 + M_PI / (2 * log(1 - epsilon_p))) * fd_p;

    //Immunological production and degradation parameters
    double K_i_p_mic, K_i_p_wound, K_i_d_max, delta_i_p, beta_i_p;
    Immune_in >> K_i_p_mic >> K_i_p_wound >> K_i_d_max >> delta_i_p >> beta_i_p;

    double ps_norm, fd_norm;
    Immune_in >> ps_norm >> fd_norm;

    double K_i_p_trans = infl_scale_trans * K_i_p_mic * (ps_p / ps_norm * fd_p / fd_norm) + K_i_p_wound;
    double K_i_p_steady = K_i_p_mic * (ps_norm / ps_p + fd_p / fd_norm - 2);
    K_i_p_steady = (K_i_p_steady > 0)* K_i_p_steady;
    double K_i_d = K_i_p_trans / K_i_d_max;

    //Polymer degradation initialization
    double Q_p1, Q_p2, Q_gnd;
    double s_gnd_off = 14.0, epsilonR_gnd_min = 0.10, s = 0.0;

    //Inflammation initialization
    double gamma_fun, steady_fun;
    ups_infl_p.resize(nts);
    ups_infl_d.resize(nts);

    for (int sn = 1; sn < nts; sn++) {
        //Calculate polymer/ground degradation
        s = sn * dt;
        Q_p1 = (1 + exp(-k_p1 * zeta_p1)) / (1 + exp(k_p1 * gamma_p_d1 * (s - zeta_p1 * gamma_p_d2)));
        Q_p2 = (1 + exp(-k_p2 * zeta_p2)) / (1 + exp(k_p2 * gamma_p_d1 * (s - zeta_p2 * gamma_p_d2)));
        Q_gnd = (s > s_gnd_off)* ((1 - epsilonR_gnd_min) * exp(-k_gnd_h * (s - s_gnd_off)) + epsilonR_gnd_min)
            + (s <= s_gnd_off) * 1;

        epsilonR_alpha[0 * nts + sn] = Q_p1 * epsilonR_alpha_0[0];
        epsilonR_alpha[1 * nts + sn] = Q_p2 * epsilonR_alpha_0[1];
        epsilonR_alpha[12 * nts + sn] = Q_gnd * epsilonR_alpha_0[12];

        rhoR_alpha[0 * nts + sn] = epsilonR_alpha[0 * nts + sn] * rho_hat_alpha_h[0];
        rhoR_alpha[1 * nts + sn] = epsilonR_alpha[1 * nts + sn] * rho_hat_alpha_h[1];
        rhoR_alpha[12 * nts + sn] = epsilonR_alpha[12 * nts + sn] * rho_hat_alpha_h[12];

        //Calculate immunological stimulus
        gamma_fun = (pow(delta_i_p, beta_i_p)) * pow(s, beta_i_p - 1) * exp(-delta_i_p * s)
            / (delta_i_p * pow((beta_i_p - 1), (beta_i_p - 1)) * exp(1 - beta_i_p));
        steady_fun = (1 - exp(-delta_i_p * s));
        ups_infl_p[sn] = K_i_p_trans* gamma_fun + K_i_p_steady * steady_fun;
        ups_infl_d[sn] = K_i_d;

    }

    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;
    lambda_z_tau.resize(nts);
    lambda_z_tau[0] = 1.0;

    //Intialize stresses
    sigma.resize(9);
    Cbar.resize(3);
    CC.resize(36);

    //Active stress parameters
    //NOTE: currently assume TEVG has no active stress contribution
    alpha_active = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    a_act.resize(nts);
    a_act[0] = a_h;
    lambda_act.resize(nts);
    lambda_act[0] = 1.0;


    F_s.resize(9*nts);
    F_curr.resize(9);
    Fi_curr.resize(9);
    //Initialize to no stretch
    F_s[0] = 1.0;
    F_s[4] = 1.0;
    F_s[8] = 1.0;
    F_curr[0] = 1.0;
    F_curr[4] = 1.0;
    F_curr[8] = 1.0;
    Fi_curr[0] = 1.0;
    Fi_curr[4] = 1.0;
    Fi_curr[8] = 1.0;
    J_di = 1.0;


    //Need to add pre-stretch to the polymeric constituents to give them and in vivo reference
    for (int alpha = 0; alpha < n_pol_alpha + 1; alpha++) {

        //Account for ground matrix the last constituent
        if (alpha == n_pol_alpha) {
            G_alpha_h[3 * (n_alpha - 1)] = 1.0;
            G_alpha_h[3 * (n_alpha - 1) + 1] = 1.0;
            G_alpha_h[3 * (n_alpha - 1) + 2] = 1.0;
        }
        else {
            G_alpha_h[3 * alpha] = 1.0;
            G_alpha_h[3 * alpha + 1] = 1.0;
            G_alpha_h[3 * alpha + 2] = 1.0;
        }
    }

    //Reset the initial stretch conditions
    lambda_th_curr = 1.0;
    lambda_z_curr = 1.0;

    sigma_inv = sigma_inv_h;
    bar_tauw = bar_tauw_h;

    //To only consider initial constituents and prescribed mass changes
    pol_only_flag = 0;

    //Update stress to in vivo stress for homeostatic stress calc
    update_sigma_handshake(this);

    GnR_out.precision(12);
    HnS_out.precision(12);

    Scaffold_in.close();
    Immune_in.close();
}



void vessel::print() {
    std::cout << "nts : " << nts << std::endl;
    std::cout << "dt : " << dt << std::endl;
    std::cout << "sn : " << sn << std::endl;
    std::cout << "s : " << s << std::endl;
    std::cout << "A_h : " << A_h << std::endl;
    std::cout << "B_h : " << B_h << std::endl;
    std::cout << "H_h : " << H_h << std::endl;
    std::cout << "A_mid_h : " << A_mid_h << std::endl;
    std::cout << "a_h : " << a_h << std::endl;
    std::cout << "b_h : " << b_h << std::endl;
    std::cout << "h_h : " << h_h << std::endl;
    std::cout << "a_mid_h : " << a_mid_h << std::endl;
    std::cout << "lambda_z_h : " << lambda_z_h << std::endl;
    std::cout << "a : " << a << std::endl;
    std::cout << "a_mid : " << a_mid << std::endl;
    std::cout << "b : " << b << std::endl;
    std::cout << "h : " << h << std::endl;
    std::cout << "A : " << A << std::endl;
    std::cout << "A_mid : " << A_mid << std::endl;
    std::cout << "B : " << B << std::endl;
    std::cout << "H : " << H << std::endl;
    std::cout << "lambda_z_pre : " << lambda_z_pre << std::endl;
    std::cout << "n_alpha : " << n_alpha << std::endl;
    std::cout << "n_pol_alpha : " << n_pol_alpha << std::endl;
    std::cout << "n_native_alpha : " << n_native_alpha << std::endl;
    std::cout << "alpha_infl : " << alpha_infl << std::endl;
    std::cout << "c_alpha_h : " << c_alpha_h << std::endl;
    std::cout << "eta_alpha_h : " << eta_alpha_h << std::endl;
    std::cout << "g_alpha_h : " << g_alpha_h << std::endl;
    std::cout << "G_alpha_h : " << G_alpha_h << std::endl;
    std::cout << "phi_alpha_h : " << phi_alpha_h << std::endl;
    std::cout << "rhoR_alpha_h : " << rhoR_alpha_h << std::endl;
    std::cout << "mR_alpha_h : " << mR_alpha_h << std::endl;
    std::cout << "k_alpha_h : " << k_alpha_h << std::endl;
    std::cout << "K_sigma_p_alpha_h : " << K_sigma_p_alpha_h << std::endl;
    std::cout << "K_sigma_d_alpha_h : " << K_sigma_d_alpha_h << std::endl;
    std::cout << "K_tauw_p_alpha_h : " << K_tauw_p_alpha_h << std::endl;
    std::cout << "K_tauw_d_alpha_h : " << K_tauw_d_alpha_h << std::endl;
    std::cout << "rho_hat_alpha_h : " << rho_hat_alpha_h << std::endl;
    std::cout << "epsilonR_alpha_0 : " << epsilonR_alpha_0 << std::endl;
    std::cout << "rhoR_h : " << rhoR_h << std::endl;
    std::cout << "ups_infl_p : " << ups_infl_p << std::endl;
    std::cout << "ups_infl_d : " << ups_infl_d << std::endl;
    return;
}

void vessel::printNativeOutputs() {
    GnR_out << s << "\t";
    GnR_out << a[sn] << "\t" << h[sn] << "\t" << rhoR[sn] << "\t" << rho[sn] << "\t" << rhoR_alpha[0 * nts + sn] << "\t"
        << rhoR_alpha[1 * nts + sn] << "\t" << rhoR_alpha[2 * nts + sn]<< "\t"
        << bar_tauw << "\t" << bar_tauw_h << "\t" << P << "\t" << P_h << "\t" << f << "\t" << f_h
        << "\t" << Q << "\t" << Q_h << "\t";
    for (int i = 0; i < 3; i++){
        GnR_out << sigma[i] << " ";
    }
    GnR_out << sigma_inv;
    GnR_out << "\n";
    GnR_out.flush();

    return;

}

void vessel::printOutputsHandshake() {
    /*
    GnR_out << s << "\t";
    GnR_out << rhoR[sn] << "\t";
    GnR_out << rho[sn] << "\t";
    for (int alpha = 0; alpha < 6; alpha++) {
        GnR_out << rhoR_alpha[alpha * nts + sn] << "\t";
    }

    GnR_out << bar_tauw << "\t" << sigma_inv << "\t";
    GnR_out << J_di << "\t";
    for (int i = 0; i < 9; i++){
        GnR_out << F_s[9*sn+i] << " ";
    }
    GnR_out << sigma << "\t";
    GnR_out << "\n";
    GnR_out.flush();
    */

    //HnS_out <<  CC << " " << J_di << " " << sigma << "\n";
    //HnS_out.flush();

    HnS_out<< std::setprecision(9)
    << s << " " << CC  << " " << sigma << " " << rhoR[sn]<< " " << rho[sn]
    << " " << F_s[9*sn+0]
    << " " << F_s[9*sn+1]
    << " " << F_s[9*sn+2]
    << " " << F_s[9*sn+3]
    << " " << F_s[9*sn+4]
    << " " << F_s[9*sn+5]
    << " " << F_s[9*sn+6]
    << " " << F_s[9*sn+7]
    << " " << F_s[9*sn+8]
    << " " << sigma_inv << " " << bar_tauw << "\n";
    HnS_out.flush();

    return;
}


void vessel::printNativeEquilibratedOutputs() {
    Equil_GnR_out << a_e << "\t" << h_e << "\t" << rho_m_e << "\t" << rho_c_e << "\t"
        << f_z_e << "\t" << mb_equil_e << "\n";
    Equil_GnR_out.flush();

    return;

}

void vessel::printTEVGOutputs() {
    GnR_out << s << "\t";
    GnR_out << a[sn] << "\t" << h[sn] << "\t" << rhoR[sn] << "\t" << rho[sn] << "\t" << rhoR_alpha[0 * nts + sn] << "\t"
        << rhoR_alpha[1 * nts + sn] << "\t" << rhoR_alpha[2 * nts + sn]<< "\t"
        << bar_tauw << "\t" << bar_tauw_h << "\t" << P << "\t" << P_h << "\t" << f << "\t" << f_h
        << "\t" << Q << "\t" << Q_h << "\t";
    for (int i = 0; i < 3; i++){
        GnR_out << sigma[i] << " ";
    }
    GnR_out << sigma_inv << " ";
    GnR_out << sigma_inv_h << " ";
    GnR_out << "\n";

    //GnR_out << a[sn] << "\t" << h[sn] << "\t" << rhoR[sn] << "\t" << rhoR_alpha[0 * nts + sn] << "\t"
    //    << rhoR_alpha[1 * nts + sn] << "\t" << rhoR_alpha[8 * nts + sn] << "\t"
    //    << bar_tauw << "\t" << bar_tauw_h << "\t" << sigma_inv << "\t" << sigma_inv_h << "\t" << f << "\t" << f_h
    //    << "\t" << Q << "\t" << Q_h << "\n";
    GnR_out.flush();

    return;
}

void vessel::load() {
    std::ifstream ar;
    ar.open(file_name);
    ar >> nts;
    ar >> dt;
    ar >> sn;
    ar >> s;
    ar >> A_h;
    ar >> B_h;
    ar >> H_h;
    ar >> A_mid_h;
    ar >> a_h;
    ar >> b_h;
    ar >> h_h;
    ar >> a_mid_h;
    ar >> lambda_z_h;
    ar >> a;
    ar >> a_mid;
    ar >> b;
    ar >> h;
    ar >> A;
    ar >> A_mid;
    ar >> B;
    ar >> H;
    ar >> a_pas;
    ar >> a_mid_pas;
    ar >> h_pas;
    ar >> lambda_z_pre;
    ar >> n_alpha;
    ar >> n_pol_alpha;
    ar >> n_native_alpha;
    ar >> alpha_infl;
    ar >> c_alpha_h;
    ar >> eta_alpha_h;
    ar >> g_alpha_h;
    ar >> G_alpha_h;
    ar >> phi_alpha_h;
    ar >> rhoR_alpha_h;
    ar >> mR_alpha_h;
    ar >> k_alpha_h;
    ar >> K_sigma_p_alpha_h;
    ar >> K_sigma_d_alpha_h;
    ar >> K_tauw_p_alpha_h;
    ar >> K_tauw_d_alpha_h;
    ar >> rho_hat_alpha_h;
    ar >> epsilonR_alpha_0;
    ar >> rhoR_h;
    ar >> rhoR;
    ar >> rho;
    ar >> rhoR_alpha;
    ar >> mR_alpha;
    ar >> k_alpha;
    ar >> epsilonR_alpha;
    ar >> epsilon_alpha;
    ar >> epsilon_pol_min;
    ar >> ups_infl_p;
    ar >> ups_infl_d;
    ar >> P_h;
    ar >> f_h;
    ar >> bar_tauw_h;
    ar >> Q_h;
    ar >> sigma_h;
    ar >> sigma_inv_h;
    ar >> lambda_th_curr;
    ar >> lambda_z_curr;
    ar >> P;
    ar >> f;
    ar >> bar_tauw;
    ar >> Q;
    ar >> sigma;
    ar >> sigma_inv;
    ar >> Cbar;
    ar >> CC;
    ar >> lambda_alpha_tau;
    ar >> lambda_z_tau;
    ar >> mb_equil;
    ar >> alpha_active;
    ar >> a_act;
    ar >> lambda_act;
    ar >> T_act;
    ar >> T_act_h;
    ar >> k_act;
    ar >> lambda_0;
    ar >> lambda_m;
    ar >> CB;
    ar >> CS;
    ar >> a_e;
    ar >> h_e;
    ar >> rho_c_e;
    ar >> rho_m_e;
    ar >> f_z_e;
    ar >> mb_equil_e;
    ar >> Ki_p_h;
    ar >> Ki_d_h;
    ar >> num_exp_flag;
    ar >> pol_only_flag;
    ar >> wss_calc_flag;
    ar >> app_visc_flag;
    ar >> F_s;
    ar >> F_curr;
    ar >> Fi_curr;
    ar >> J_di;
    ar >> P_prev;
    ar.close();

    //std::cout << v;
    //std::cout << "Loaded vessel." << "\n";

}

void vessel::save() {
    std::ofstream ar;
    ar.open(file_name);
    ar.precision(20);

    ar << nts << "\n";
    ar << dt << "\n";
    ar << sn << "\n";
    ar << s << "\n";
    ar << A_h << "\n";
    ar << B_h << "\n";
    ar << H_h << "\n";
    ar << A_mid_h << "\n";
    ar << a_h << "\n";
    ar << b_h << "\n";
    ar << h_h << "\n";
    ar << a_mid_h << "\n";
    ar << lambda_z_h << "\n";
    ar << a << "\n";
    ar << a_mid << "\n";
    ar << b << "\n";
    ar << h << "\n";
    ar << A << "\n";
    ar << A_mid << "\n";
    ar << B << "\n";
    ar << H << "\n";
    ar << a_pas << "\n";;
    ar << a_mid_pas << "\n";;
    ar << h_pas << "\n";;
    ar << lambda_z_pre << "\n";
    ar << n_alpha << "\n";
    ar << n_pol_alpha << "\n";
    ar << n_native_alpha << "\n";
    ar << alpha_infl << "\n";
    ar << c_alpha_h << "\n";
    ar << eta_alpha_h << "\n";
    ar << g_alpha_h << "\n";
    ar << G_alpha_h << "\n";
    ar << phi_alpha_h << "\n";
    ar << rhoR_alpha_h << "\n";
    ar << mR_alpha_h << "\n";
    ar << k_alpha_h << "\n";
    ar << K_sigma_p_alpha_h << "\n";
    ar << K_sigma_d_alpha_h << "\n";
    ar << K_tauw_p_alpha_h << "\n";
    ar << K_tauw_d_alpha_h << "\n";
    ar << rho_hat_alpha_h << "\n";
    ar << epsilonR_alpha_0 << "\n";
    ar << rhoR_h << "\n";
    ar << rhoR << "\n";
    ar << rho << "\n";
    ar << rhoR_alpha << "\n";
    ar << mR_alpha << "\n";
    ar << k_alpha << "\n";
    ar << epsilonR_alpha << "\n";
    ar << epsilon_alpha << "\n";
    ar << epsilon_pol_min << "\n";
    ar << ups_infl_p << "\n";
    ar << ups_infl_d << "\n";
    ar << P_h << "\n";
    ar << f_h << "\n";
    ar << bar_tauw_h << "\n";
    ar << Q_h << "\n";
    ar << sigma_h << "\n";
    ar << sigma_inv_h << "\n";
    ar << lambda_th_curr << "\n";
    ar << lambda_z_curr << "\n";
    ar << P << "\n";
    ar << f << "\n";
    ar << bar_tauw << "\n";
    ar << Q << "\n";
    ar << sigma << "\n";
    ar << sigma_inv << "\n";
    ar << Cbar << "\n";
    ar << CC << "\n";
    ar << lambda_alpha_tau << "\n";
    ar << lambda_z_tau << "\n";
    ar << mb_equil << "\n";
    ar << alpha_active << "\n";
    ar << a_act << "\n";
    ar << lambda_act << "\n";
    ar << T_act << "\n";
    ar << T_act_h << "\n";
    ar << k_act << "\n";
    ar << lambda_0 << "\n";
    ar << lambda_m << "\n";
    ar << CB << "\n";
    ar << CS << "\n";
    ar << a_e << "\n";
    ar << h_e << "\n";
    ar << rho_c_e << "\n";
    ar << rho_m_e << "\n";
    ar << f_z_e << "\n";
    ar << mb_equil_e << "\n";
    ar << Ki_p_h << "\n";
    ar << Ki_d_h << "\n";
    ar << num_exp_flag << "\n";
    ar << pol_only_flag << "\n";
    ar << wss_calc_flag << "\n";
    ar << app_visc_flag << "\n";
    ar << F_s << "\n";
    ar << F_curr << "\n";
    ar << Fi_curr << "\n";
    ar << J_di << "\n";
    ar << P_prev << "\n";

    ar.close();
    //std::cout << "Saved vessel." << "\n";
};

