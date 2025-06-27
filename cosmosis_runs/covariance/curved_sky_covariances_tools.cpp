#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_poly.h>
#include <algorithm> 
#include <cassert>

namespace constants {
	
	double pi = 4.0*atan(1.0);
	double arcmin = pi/180.0/60.0;
	double degree = pi/180.0;
	
};


using namespace std;
using namespace constants;

bool LSZ_correction_on = 1;
int Npair_modus = 0;
/* Npair_modus = 0: everything as usual ; Npair_modus = 1: masking wtheta ignored ; Npair_modus = 2: optimized for FLASK. */


/*
 * set_angular_bins:
 * Returns the edges of angular bins that are log-spaced and in radians.
 *
 * Params:
 * - double theta_min : minimal angular scale
 * - double theta_max : maximal angular scale
 * - int N_bin : number of angular bins
 *
 */
vector<double> set_angular_bins_in_radians(double theta_min, double theta_max, int N_bin){
  vector<double> bins(N_bin+1, theta_min*arcmin);
  double alpha = pow(theta_max/theta_min, 1.0/double(N_bin));
  
  for(int i = 1; i < N_bin+1; i++){
    bins[i] = bins[i-1]*alpha;
  }
  
  return bins;
  
}

vector<double> set_linear_angular_bins_in_radians(double theta_min, double theta_max, int N_bin){
  vector<double> bins(N_bin+1, theta_min*arcmin);
  double dtheta = (theta_max-theta_min)/double(N_bin);
  
  for(int i = 1; i < N_bin+1; i++){
    bins[i] = bins[i-1]+dtheta;
  }
  
  return bins;
  
}


/*
 * return_Legendres_at_one_theta:
 * Returns Legendre polynomials and their derivatives at x = cos(theta).
 *
 * Params:
 * - int n_ell : (maximum value of ell) + 1
 * - double cos_theta : argument of Legendre polynomials
 * - vector<double> *Pl : pointer to vector in which polynomial values will be stored
 * - vector<double> *Pl_deriv : pointer to vector in which derivative values will be stored
 *
 */
void return_Legendres_at_one_theta(int n_ell, double cos_theta, vector<double> *Pl,  vector<double> *Pl_deriv){
  
  (*Pl) = vector<double>(n_ell, 0.0);
  (*Pl_deriv) = vector<double>(n_ell, 0.0);
  
  
  double *Pl_th_aux = new double[n_ell+3];
  double *Pl_th_derivs_aux = new double[n_ell+3];
  
  gsl_sf_legendre_Pl_deriv_array(n_ell-1, cos_theta, Pl_th_aux, Pl_th_derivs_aux);
  
  for(int l = 0; l < n_ell; l++){
    (*Pl)[l] = Pl_th_aux[l];
    (*Pl_deriv)[l] = Pl_th_derivs_aux[l];
  }
  
  delete[] Pl_th_aux;
  delete[] Pl_th_derivs_aux;
  
}





/*
 * return_associated_Legendres_at_one_theta:
 * Returns associated Legendre polynomials and their derivatives at x = cos(theta).
 *
 * Params:
 * - int n_ell : (maximum value of ell) + 1
 * - int m: degree of association (is that even a thing?)
 * - double cos_theta : argument of Legendre polynomials
 * - vector<double> *Pl : pointer to vector in which polynomial values will be stored
 * - vector<double> *Pl_deriv : pointer to vector in which derivative values will be stored
 *
 */
void return_associated_Legendres_at_one_theta(int n_ell, int m, double cos_theta, vector<double> *Pl2){
  
  (*Pl2) = vector<double>(n_ell, 0.0);
  
  long int n = gsl_sf_legendre_array_index(n_ell, n_ell);
  double *Plm_th_aux = new double[n+1];
  
  gsl_sf_legendre_array(GSL_SF_LEGENDRE_NONE, n_ell-1, cos_theta, Plm_th_aux);
  
  for(int l = m; l < n_ell; l++){
    int i = gsl_sf_legendre_array_index(l, m);
    (*Pl2)[l] = Plm_th_aux[i];
  }
  
  delete[] Plm_th_aux;
  
}






vector<vector<double> > return_Legendres_small_bins(int nu, int shear, double ell_max, vector<double> bins){

  int n_ell = ell_max+1;
  int n_th = bins.size()-1;
    
  vector<vector<double> > Legendres(n_th, vector<double>(n_ell, 0.0));
  
  double x;
  
  if(nu == 0){
    for(int th = 0; th < n_th; th++){
      x = 2.0/3.0*(pow(bins[th+1], 3)-pow(bins[th], 3))/(pow(bins[th+1], 2)-pow(bins[th], 2));
      for(int l = 0; l < n_ell; l++){
        Legendres[th][l] = gsl_sf_bessel_J0(x*double(l));
      }
    }
  }
  
  if(nu == 2){
    for(int th = 0; th < n_th; th++){
      x = 2.0/3.0*(pow(bins[th+1], 3)-pow(bins[th], 3))/(pow(bins[th+1], 2)-pow(bins[th], 2));
      for(int l = 1; l < n_ell; l++){
        Legendres[th][l] = gsl_sf_bessel_Jn(2, x*double(l));
      }
    }
  }
  
  if(nu == 4){
    for(int th = 0; th < n_th; th++){
      x = 2.0/3.0*(pow(bins[th+1], 3)-pow(bins[th], 3))/(pow(bins[th+1], 2)-pow(bins[th], 2));
      for(int l = 1; l < n_ell; l++){
        Legendres[th][l] = gsl_sf_bessel_Jn(4, x*double(l));
      }
    }
  }
  
  
  return Legendres;

}


/*
 * return_Legendres_finite_bins:
 * Returns Legendre-related functions F_{bin}(\ell) such that
 * \xi(bin) = \sum_\ell (2.0*\ell + 1.0)/(4.0*\pi) * C_\ell * F_{bin}(\ell)
 * --> NOTE: this function assumes that the number of galaxy pairs
 * --> with separation \theta goes as ~ 2\pi(1-\cos\theta)
 * 
 * Params:
 * - int nu : == 0 for \xi_+ or w(\theta) ; == 4 for \xi_- ; == 2 for \gamma_t
 * - int shear : == 1 for cosmic shear correlation functions ; == 0 else
 * - int ell_max : maximum value of \ell
 * - vector<double> bins : vector of bin EDGES in radians
 * 
 * returns:
 * - vector<vector<double> > Legendres : 2-dim array storing F_{bin}(\ell) where first dimension labels the angular bins and second dimension labels \ell
 * 
 */

vector<vector<double> > return_Legendres_finite_bins(int nu, int shear, int ell_max, vector<double> bins){

  int n_ell = ell_max+1;
  int n_th = bins.size();
  
  double x;
  double alpha = pow(bins[n_th-1]/bins[0], 1.0/double(n_th-1));
  
  vector<vector<double> > Legendres_coefficients(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > Legendres(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > Legendres_derivs(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > Legendres_l2(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > Pl_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  vector<vector<double> > Pl_deriv_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  vector<vector<double> > Pl2_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  vector<vector<double> > Pl2_over_1_min_x2_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  vector<vector<double> > Pl2_times_x_over_1_min_x2_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  vector<vector<double> > G_minus_l2(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > G_plus_l2(n_th, vector<double>(n_ell, 0.0));
  
  vector<double> Pl_th_aux(n_ell, 0.0);
  vector<double> Pl2_th_aux(n_ell, 0.0);
  vector<double> Pl_th_derivs_aux(n_ell, 0.0);
  
  for(int th = 0; th < n_th; th++){

    return_Legendres_at_one_theta(n_ell, cos(bins[th]), &Pl_th_aux, &Pl_th_derivs_aux);
    //return_associated_Legendres_at_one_theta(n_ell, 2, cos(bins[th]), &Legendres_l2[th]);

    for(int l = 0; l < n_ell; l++){ 
      Legendres[th][l] = Pl_th_aux[l];
      Legendres_derivs[th][l] = Pl_th_derivs_aux[l];
    }
    
  }
  
  for(int th = 0; th < n_th-1; th++){
    
    double x1 = cos(bins[th]);
    double x2 = cos(bins[th+1]);
    double x1_sq = x1*x1;
    double x2_sq = x2*x2;
    double normalisation = x1-x2;
    double ell;
    for(int l = 1; l < n_ell-1; l++){
      ell = double(l);
      
      Pl_bin_average[th][l] = (Legendres[th][l+1] - Legendres[th][l-1])/(2.0*ell+1);
      Pl_bin_average[th][l] -= (Legendres[th+1][l+1] - Legendres[th+1][l-1])/(2.0*ell+1);
      Pl_bin_average[th][l] /= normalisation;
      
      Pl2_bin_average[th][l] = (1.0-x1_sq)*Legendres_derivs[th][l];
      Pl2_bin_average[th][l] -= (1.0-x2_sq)*Legendres_derivs[th+1][l];
      Pl2_bin_average[th][l] += 2.0*x1*Legendres[th][l];
      Pl2_bin_average[th][l] -= 2.0*x2*Legendres[th+1][l];
      Pl2_bin_average[th][l] -= 2.0/(2.0*ell+1)*(Legendres[th][l+1] - Legendres[th][l-1]);
      Pl2_bin_average[th][l] += 2.0/(2.0*ell+1)*(Legendres[th+1][l+1] - Legendres[th+1][l-1]);
      Pl2_bin_average[th][l] /= normalisation;
      
      // testing more efficient bin average for P_l^2:
      Pl2_bin_average[th][l] = 2.0*x1*Legendres[th][l];
      Pl2_bin_average[th][l] -= 2.0*x2*Legendres[th+1][l];
      Pl2_bin_average[th][l] /= normalisation;
      Pl2_bin_average[th][l] -= (2.0+ell*(ell+1))*Pl_bin_average[th][l];
      
      
      Pl2_over_1_min_x2_bin_average[th][l] = Legendres_derivs[th][l];
      Pl2_over_1_min_x2_bin_average[th][l] -= Legendres_derivs[th+1][l];
      Pl2_over_1_min_x2_bin_average[th][l] /= normalisation;
      
      Pl2_times_x_over_1_min_x2_bin_average[th][l] = x1*Legendres_derivs[th][l];
      Pl2_times_x_over_1_min_x2_bin_average[th][l] -= x2*Legendres_derivs[th+1][l];
      Pl2_times_x_over_1_min_x2_bin_average[th][l] -= Legendres[th][l];
      Pl2_times_x_over_1_min_x2_bin_average[th][l] += Legendres[th+1][l];
      Pl2_times_x_over_1_min_x2_bin_average[th][l] /= normalisation;

      
    }
    Pl2_bin_average[th][1] = 0.0;
    Pl2_over_1_min_x2_bin_average[th][1] = 0.0;
    Pl2_times_x_over_1_min_x2_bin_average[th][1] = 0.0;
  }
  
  double Legendre_new;
  
  // wtheta
  if(nu == 0 && shear == 0){
    for(int th = 0; th < n_th-1; th++){
      for(int l = 1; l < n_ell-1; l++){
        Legendres_coefficients[th][l] = Pl_bin_average[th][l];
      }
      Legendres_coefficients[th][0] = 1.0;
    }
  }
  
  // ggl
  if(nu == 2){
    for(int th = 0; th < n_th-1; th++){
      double ell;
      for(int l = 1; l < n_ell-1; l++){
        ell = double(l);
        Legendres_coefficients[th][l] = Pl2_bin_average[th][l]/(ell*(ell+1));
      }
      Legendres_coefficients[th][0] = 0.0;
      Legendres_coefficients[th][1] = 0.0;
    }
  }
  
  // xi_plus
  if(nu == 0 && shear == 1){
    for(int th = 0; th < n_th-1; th++){
      double ell;
      for(int l = 1; l < n_ell-1; l++){
        ell = double(l);
        
        Legendres_coefficients[th][l] = -ell*(ell-1.0)/2.0*Pl2_bin_average[th][l];
        Legendres_coefficients[th][l] += (4.0-ell)*Pl2_over_1_min_x2_bin_average[th][l];
        Legendres_coefficients[th][l] += (ell+2.0)*Pl2_times_x_over_1_min_x2_bin_average[th][l-1];
        
        // These contributions change sign between xi_plus and xi_minus:
        Legendres_coefficients[th][l] += 2.0*(ell-1.0)*Pl2_times_x_over_1_min_x2_bin_average[th][l];
        Legendres_coefficients[th][l] -= 2.0*(ell+2.0)*Pl2_over_1_min_x2_bin_average[th][l-1];

        Legendres_coefficients[th][l] /= 0.5*pow(ell*(ell+1.0), 2);
        
      }
      Legendres_coefficients[th][0] = 0.0;
      Legendres_coefficients[th][1] = 0.0;
      //Legendres_coefficients[th][2] = 0.0;
    }
  }
  
  // xi_minus
  if(nu == 4 && shear == 1){
    for(int th = 0; th < n_th-1; th++){
      double ell;
      
      for(int l = 1; l < n_ell-1; l++){
        ell = double(l);
        
        Legendres_coefficients[th][l] = -ell*(ell-1.0)/2.0*Pl2_bin_average[th][l];
        Legendres_coefficients[th][l] += (4.0-ell)*Pl2_over_1_min_x2_bin_average[th][l];
        Legendres_coefficients[th][l] += (ell+2.0)*Pl2_times_x_over_1_min_x2_bin_average[th][l-1];
        
        // These contributions change sign between xi_plus and xi_minus:
        Legendres_coefficients[th][l] -= 2.0*(ell-1.0)*Pl2_times_x_over_1_min_x2_bin_average[th][l];
        Legendres_coefficients[th][l] += 2.0*(ell+2.0)*Pl2_over_1_min_x2_bin_average[th][l-1];
        
        Legendres_coefficients[th][l] /= 0.5*pow(ell*(ell+1.0), 2);
        
      }
      Legendres_coefficients[th][0] = 0.0;
      Legendres_coefficients[th][1] = 0.0;
    }
  }
  
  return Legendres_coefficients;

}

vector<vector<double> > return_Legendres_finite_bins_flat_sky(int nu, int shear, int ell_max, vector<double> bins){

  int n_ell = ell_max+1;
  int n_th = bins.size();
  
  double x;
  double alpha = pow(bins[n_th-1]/bins[0], 1.0/double(n_th-1));
  
  vector<vector<double> > Bessel_coefficients(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > J0(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > J1(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > J2(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > J3(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > J4(n_th, vector<double>(n_ell, 0.0));
  vector<vector<double> > J0_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  vector<vector<double> > J2_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  vector<vector<double> > J4_bin_average(n_th-1, vector<double>(n_ell, 0.0));
  
  vector<double> Pl_th_aux(n_ell, 0.0);
  vector<double> Pl2_th_aux(n_ell, 0.0);
  vector<double> Pl_th_derivs_aux(n_ell, 0.0);
  
  for(int th = 0; th < n_th; th++){
    x = bins[th];
    for(int l = 0; l < n_ell; l++){ 
      J0[th][l] = gsl_sf_bessel_J0(x*double(l));
      J1[th][l] = gsl_sf_bessel_J1(x*double(l));
      J2[th][l] = gsl_sf_bessel_Jn(2, x*double(l));
      J3[th][l] = gsl_sf_bessel_Jn(3, x*double(l));
      J4[th][l] = gsl_sf_bessel_Jn(4, x*double(l));
    }
  }
  
  double x1;
  double x2;
  double x1_sq;
  double x2_sq;
  double normalisation;
  double ell;
  
  for(int th = 0; th < n_th-1; th++){
    
    double x1 = bins[th];
    double x2 = bins[th+1];
    double x1_sq = x1*x1;
    double x2_sq = x2*x2;
    double normalisation = 0.5*(x1_sq-x2_sq);
    double ell;
    for(int l = 1; l < n_ell-1; l++){
      double ell = double(l);
      x1 = ell*bins[th];
      x2 = ell*bins[th+1];
      x1_sq = x1*x1;
      x2_sq = x2*x2;
      normalisation = 0.5*(x1_sq-x2_sq);
      
      J0_bin_average[th][l] = x1*J1[th][l];
      J0_bin_average[th][l] -= x2*J1[th+1][l];
      J0_bin_average[th][l] /= normalisation;
      
      J2_bin_average[th][l] = -2.0*J0[th][l];
      J2_bin_average[th][l] -= -2.0*J0[th+1][l];
      J2_bin_average[th][l] /= normalisation;
      J2_bin_average[th][l] -= J0_bin_average[th][l];
      
      J4_bin_average[th][l] = 2.0*x1*(J1[th][l]-J3[th][l]);
      J4_bin_average[th][l] -= 2.0*x2*(J1[th+1][l]-J3[th+1][l]);
      J4_bin_average[th][l] += -2.0*J0[th][l];
      J4_bin_average[th][l] -= -2.0*J0[th+1][l];
      J4_bin_average[th][l] += -4.0*J2[th][l];
      J4_bin_average[th][l] -= -4.0*J2[th+1][l];
      J4_bin_average[th][l] /= normalisation;
      J4_bin_average[th][l] += J2_bin_average[th][l];
      J4_bin_average[th][l] -= 2.0*J0_bin_average[th][l];
      
      
    }
  }
  
  double Legendre_new;
  
  if(nu == 0 && shear == 0){
    for(int th = 0; th < n_th-1; th++){
      for(int l = 1; l < n_ell-1; l++){
        Bessel_coefficients[th][l] = J0_bin_average[th][l];
      }
      Bessel_coefficients[th][0] = 1.0;
    }
  }
  
  if(nu == 2){
    for(int th = 0; th < n_th-1; th++){
      for(int l = 1; l < n_ell-1; l++){
        Bessel_coefficients[th][l] = J2_bin_average[th][l];
      }
    }
  }
  
  if(nu == 0 && shear == 1){
    for(int th = 0; th < n_th-1; th++){
      for(int l = 1; l < n_ell-1; l++){        
        Bessel_coefficients[th][l] = J0_bin_average[th][l];
      }
      Bessel_coefficients[th][0] = 1.0;
    }
  }
  
  if(nu == 4 && shear == 1){
    for(int th = 0; th < n_th-1; th++){
      for(int l = 1; l < n_ell-1; l++){
        Bessel_coefficients[th][l] = J4_bin_average[th][l];
      }
    }
  }
  
  return Bessel_coefficients;

}

/*
 * xi_curved_sky:
 * Returns value of 2-point correlation functions.
 * 
 * Params:
 * - int i_th : index of angular bin
 * - vector<double> *Cell : pointer to vector of C_ell
 * - vector<vector<double> >* Legendres : pointer to the result of the function return_Legendres_finite_bins
 * 
 * returns:
 * - double xi : value of 2-point function at a particular angular bin
 * 
 */
double xi_curved_sky(int i_th, vector<double> *Cell, vector<vector<double> >* Legendres){
 
  int n_ell = std::min((*Legendres)[i_th].size(),(*Cell).size());
  double xi = 0.0;
  double ell;
  double Legendre_factor;
  
  for(int l = 0; l < n_ell; l++){
    ell = double(l);
    Legendre_factor = (*Legendres)[i_th][l];
    xi += (2.0*ell + 1.0)/(4.0*pi)*Legendre_factor*(*Cell)[l];
  }
    
  return xi;
  
}


/*
 * get_lognormal_shift:
 * Returns shift parameter of a zero-mean-shifted lognormal PDF from the variance and skewness of the PDF.
 * NOTE: Correlation function estimators do NOT have 0 mean, so the PDF still has to be shifted in the end.
 * 
 * Params:
 * - double twoPt : 2nd central moment of the PDF (Variance)
 * - double threePt : 3rd central moment of the PDF
 * 
 * returns:
 * - double x0 : shift parameter of a zero-mean-shifted lognormal PDF
 * 
 */
double get_lognormal_shift(double twoPt, double threePt){
  
  double S3 = threePt/twoPt/twoPt;
  double x0, x1, x2;
  gsl_poly_solve_cubic(-3.0/S3, 0.0, -twoPt/S3, &x0, &x1, &x2);
  return x0;
  
}




/*
 * return_Cells_from_file:
 * Read and return C_ell values from a file.
 * 
 * Params:
 * - string Cell_file : Path to your C_ell file ; File should be of the format: column 1 == ell, column 2 == C_ell
 * 
 * returns:
 * - vector<double> Cell : vector with C_ell values
 * 
 */
vector<double> return_Cells_from_file(string Cell_file){

  vector<double> Cell(0, 0.0);
  vector<double> dummy(2, 0.0);
  
  fstream in;
  in.open(Cell_file.c_str());
  
  int count = 0;
  while(in.good()){
    in >> dummy[0];
    if(in.good()){
      if(count == 0){
        for(int ell = 0; ell < dummy[0]; ell++){
          Cell.push_back(0.0);
        }
        count++;
      }
      in >> dummy[1];
      //if(dummy[0]<10.) dummy[1]=0.0;
      Cell.push_back(dummy[1]);
    }
  }
  in.close();
  
  return Cell;
  
}


enum TYPE_OF_FIELD {GALAXY_DENSITY, SHEAR_FIELD};
const char *TYPE_OF_FIELD_types[] = {"GALAXY_DENSITY", "SHEAR_FIELD"};

enum TYPE_OF_2POINT_FUNCTION {WTHETA, GAMMA_T, XI_PLUS, XI_MINUS};
const char *TYPE_OF_2POINT_FUNCTION_types[] = {"WTHETA", "GAMMA_T", "XI_PLUS", "XI_MINUS"};
// use: static_cast<PDF_MODUS>(PDF_modus)







/*
 * covariance_curved_sky:
 * Returns covariance matrix of 2-point correlation functions.
 * 
 * Params:
 * - SO MANY!!!
 * 
 * returns:
 * - double cov : value of covariance matrix for a particular bin pair
 * 
 */
double cosmic_variance_curved_sky(int i_th1, int i_th2, double f_sky, vector<double> *Cell_13, vector<double> *Cell_24, vector<double> *Cell_14, vector<double> *Cell_23, vector<vector<double> >* Legendres_1, vector<vector<double> >* Legendres_2){
 
  int n_ell = (*Cell_13).size();
  
  double cov = 0.0;
  double Legendre_1;
  double Legendre_2;
  double sigma_sq;
  
  for(int l = 0; l < n_ell-1; l++){
    Legendre_1 =  (*Legendres_1)[i_th1][l];
    Legendre_2 =  (*Legendres_2)[i_th2][l];
    sigma_sq  = (*Cell_13)[l]*(*Cell_24)[l];
    sigma_sq += (*Cell_14)[l]*(*Cell_23)[l];
    sigma_sq /= (2.0*double(l)+1.0);
    cov += Legendre_1*Legendre_2*sigma_sq*pow(2.0*double(l)+1.0, 2);
  }
  
  cov *= 1.0/(f_sky*pow(4.0*pi, 2));
  
  return cov;
  
}




/*
 * covariance_curved_sky:
 * Returns covariance matrix of 2-point correlation functions.
 * 
 * Params:
 * - SO MANY!!!
 * 
 * returns:
 * - double cov : value of covariance matrix for a particular bin pair
 * 
 */
void covariance_curved_sky(double *covariance, double *Gaussian_part, double *Cov_SS, double *Cov_SN, double *Cov_NN, int i_th1, int i_th2, int i, int j, double f_sky, double noise_1, double noise_2, double noise_3, double noise_4, int field_index_1, int field_index_2, int field_index_3, int field_index_4, TYPE_OF_FIELD field_type_1, TYPE_OF_FIELD field_type_2, TYPE_OF_FIELD field_type_3, TYPE_OF_FIELD field_type_4, vector<double> *bin_edges, vector<double> *Cell_13, vector<double> *Cell_24, vector<double> *Cell_14, vector<double> *Cell_23, vector<vector<double> >* Legendres_1, vector<vector<double> >* Legendres_2, vector<double> *data_vector, vector<vector<double> > *lognormal_terms, vector<double> *mask_wtheta){
 
  int n_bins = (*bin_edges).size()-1;
  int n_ell = (*Legendres_1)[0].size();
  int ell_max = (*Cell_23).size();

  double cov = 0.0;
  double Legendre_1;
  double Legendre_2;
  double sigma_sq = 0.0;
  
  double noise_13 = 0.0;
  double noise_24 = 0.0;
  double noise_14 = 0.0;
  double noise_23 = 0.0;
  
  double n_random = 1.36/(constants::arcmin*constants::arcmin);
  
  if(field_index_1 == field_index_3){
    noise_13 = noise_1*noise_3;
  }
  else if(field_type_1 == GALAXY_DENSITY && field_type_3 == GALAXY_DENSITY && LSZ_correction_on){
    noise_13 = 1.0/n_random; // we use the same random points in all bins!
  }
  if(field_index_2 == field_index_4){
    noise_24 = noise_2*noise_4;
  }
  else if(field_type_2 == GALAXY_DENSITY && field_type_4 == GALAXY_DENSITY && LSZ_correction_on){
    noise_24 = 1.0/n_random; // we use the same random points in all bins!
  }
  if(field_index_1 == field_index_4){
    noise_14 = noise_1*noise_4;
  }
  else if(field_type_1 == GALAXY_DENSITY && field_type_4 == GALAXY_DENSITY && LSZ_correction_on){
    noise_14 = 1.0/n_random; // we use the same random points in all bins!
  }
  if(field_index_2 == field_index_3){
    noise_23 = noise_2*noise_3;
  }
  else if(field_type_2 == GALAXY_DENSITY && field_type_3 == GALAXY_DENSITY && LSZ_correction_on){
    noise_23 = 1.0/n_random; // we use the same random points in all bins!
  }
  
  
  double cov_signal_signal = 0.0;
  double cov_noise_signal = 0.0;
  
  for(int l = 0; l < n_ell-1; l++){
    Legendre_1 =  (*Legendres_1)[i_th1][l];
    Legendre_2 =  (*Legendres_2)[i_th2][l];
    sigma_sq = 0.0;
    if(l < ell_max-1){
      sigma_sq = (*Cell_13)[l]*(*Cell_24)[l] + (*Cell_14)[l]*(*Cell_23)[l];
      sigma_sq /= (2.0*double(l)+1.0);
      cov_signal_signal += Legendre_1*Legendre_2*sigma_sq*pow(2.0*double(l)+1.0, 2);
      sigma_sq = (*Cell_13)[l]*noise_24 + noise_13*(*Cell_24)[l] + (*Cell_14)[l]*noise_23 + noise_14*(*Cell_23)[l];
      sigma_sq /= (2.0*double(l)+1.0);
      cov_noise_signal += Legendre_1*Legendre_2*sigma_sq*pow(2.0*double(l)+1.0, 2);
      //sigma_sq  += ((*Cell_13)[l] + noise_13)*((*Cell_24)[l] + noise_24) - noise_13*noise_24;
      //sigma_sq += ((*Cell_14)[l] + noise_14)*((*Cell_23)[l] + noise_23) - noise_14*noise_23;
    }
    //sigma_sq /= (2.0*double(l)+1.0);
    //cov += Legendre_1*Legendre_2*sigma_sq*pow(2.0*double(l)+1.0, 2);
  }
  
  
  cov_noise_signal *= 1.0/(f_sky*pow(4.0*pi, 2));
  cov_signal_signal *= 1.0/(f_sky*pow(4.0*pi, 2));
  cov = cov_noise_signal + cov_signal_signal;
  
  double A_sky = 4.0*pi;
  double A_bin = 2.0*pi*(cos((*bin_edges)[i_th1]) - cos((*bin_edges)[i_th1+1]));
  
  
  double N_pair_ratio = 1.0;
  //Npair_modus = 0;
  if(Npair_modus == 1){
    (*mask_wtheta)[i_th1] = f_sky;
  }
  //cout << "HAHAHA\n";
  //cov = 0.0;
  
  (*Cov_NN) = 0.0;
  
  if(i == j && field_index_1 == field_index_2 && field_type_1 == SHEAR_FIELD && field_type_2 == SHEAR_FIELD){
    double sigma_eps_complex_pow4_over_n_squared = noise_1*noise_2*noise_1*noise_2*2.0;
    double N_pair_for_unit_density = A_sky*A_bin/2.0;
    N_pair_for_unit_density *= (*mask_wtheta)[i_th1];
    N_pair_for_unit_density *= N_pair_ratio;
    (*Cov_NN) += sigma_eps_complex_pow4_over_n_squared/N_pair_for_unit_density;
  }
  if(i == j && field_index_1 != field_index_2 && field_type_1 == SHEAR_FIELD && field_type_2 == SHEAR_FIELD){
    double sigma_eps_complex_pow4_over_n_squared = noise_1*noise_1*noise_2*noise_2*2.0;
    double N_pair_for_unit_density = A_sky*A_bin;
    N_pair_for_unit_density *= (*mask_wtheta)[i_th1];
    N_pair_for_unit_density *= N_pair_ratio;
    (*Cov_NN) += sigma_eps_complex_pow4_over_n_squared/N_pair_for_unit_density;
  }
  if(i == j && field_index_1 == field_index_2 && field_type_1 == GALAXY_DENSITY && field_type_2 == GALAXY_DENSITY){
    double sigma_eps_complex_pow4_over_n_squared = noise_1*noise_2*noise_1*noise_2;
    double N_pair_for_unit_density = A_sky*A_bin/2.0;
    N_pair_for_unit_density *= (*mask_wtheta)[i_th1];
    N_pair_for_unit_density *= N_pair_ratio;
    (*Cov_NN) += sigma_eps_complex_pow4_over_n_squared/N_pair_for_unit_density;
  }
  if(i == j && field_index_1 != field_index_2 && field_type_1 == GALAXY_DENSITY && field_type_2 == GALAXY_DENSITY){
    double sigma_eps_complex_pow4_over_n_squared = noise_1*noise_2*noise_1*noise_2;
    double N_pair_for_unit_density = A_sky*A_bin;
    N_pair_for_unit_density *= (*mask_wtheta)[i_th1];
    N_pair_for_unit_density *= N_pair_ratio;
    (*Cov_NN) += sigma_eps_complex_pow4_over_n_squared/N_pair_for_unit_density;
  }
  if(i != j && LSZ_correction_on && i%n_bins == j%n_bins && field_type_1 == GALAXY_DENSITY && field_type_2 == GALAXY_DENSITY && field_type_3 == GALAXY_DENSITY && field_type_4 == GALAXY_DENSITY){
    // THIS IS NEEDED BECAUSE WE USE SAME RANDOMS IN EACH BIN
    double sigma_eps_complex_pow4_over_n_squared = 1.0/pow(n_random,2);
    double N_pair_for_unit_density = A_sky*A_bin/2.0;
    N_pair_for_unit_density *= (*mask_wtheta)[i_th1];
    N_pair_for_unit_density *= N_pair_ratio;
    (*Cov_NN) += sigma_eps_complex_pow4_over_n_squared/N_pair_for_unit_density;
  }
  if(i == j && field_type_1 == GALAXY_DENSITY && field_type_2 == SHEAR_FIELD){
    double sigma_eps_complex_pow4_over_n_squared = noise_1*noise_1*noise_2*noise_2;
    double N_pair_for_unit_density = A_sky*A_bin;
    N_pair_for_unit_density *= (*mask_wtheta)[i_th1];
    N_pair_for_unit_density *= N_pair_ratio;
    (*Cov_NN) += sigma_eps_complex_pow4_over_n_squared/N_pair_for_unit_density;
  }
  if(i != j && LSZ_correction_on && i%n_bins == j%n_bins && field_type_1 == GALAXY_DENSITY && field_type_2 == SHEAR_FIELD && field_type_3 == GALAXY_DENSITY && field_type_4 == SHEAR_FIELD && field_index_2 == field_index_4){
    // THIS IS NEEDED BECAUSE WE USE SAME RANDOMS IN EACH BIN
    double sigma_eps_complex_pow4_over_n_squared = noise_2*noise_4/n_random;
    double N_pair_for_unit_density = A_sky*A_bin;
    N_pair_for_unit_density *= (*mask_wtheta)[i_th1];
    N_pair_for_unit_density *= N_pair_ratio;
    (*Cov_NN) += sigma_eps_complex_pow4_over_n_squared/N_pair_for_unit_density;
  }
  
  cov += (*Cov_NN);
  (*Gaussian_part) = cov;
  (*Cov_SS) = cov_signal_signal;
  (*Cov_SN) = cov_noise_signal;
  
  cov += (*data_vector)[i]*(*data_vector)[j]*(*lognormal_terms)[field_index_1][field_index_3];
  cov += (*data_vector)[i]*(*data_vector)[j]*(*lognormal_terms)[field_index_2][field_index_4];
  cov += (*data_vector)[i]*(*data_vector)[j]*(*lognormal_terms)[field_index_1][field_index_4];
  cov += (*data_vector)[i]*(*data_vector)[j]*(*lognormal_terms)[field_index_2][field_index_3];
  
  (*covariance) = cov;
  
}




/*
 * variance_and_skewness:
 * Returns covariance matrix of 2-point correlation functions.
 * 
 * Params:
 * - SO MANY!!!
 * 
 * returns:
 * - double cov : value of covariance matrix for a particular bin pair
 * 
 */
void variance_and_skewness(double *variance_value, double *skewness_value, int i_th, double f_sky, double noise_1, double noise_2, double noise_3, double noise_4, int field_index_1, int field_index_2, int field_index_3, int field_index_4, TYPE_OF_FIELD field_type_1, TYPE_OF_FIELD field_type_2, TYPE_OF_FIELD field_type_3, TYPE_OF_FIELD field_type_4, vector<double> *bin_edges, vector<double> *Cell_13, vector<double> *Cell_24, vector<double> *Cell_14, vector<double> *Cell_23, vector<vector<double> >* Legendres, vector<double> *data_vector, vector<vector<double> > *lognormal_terms, vector<double> *mask_wtheta){
 
  int n_bins = (*bin_edges).size()-1;
  int n_ell = (*Legendres)[0].size();
  int ell_max = (*Cell_23).size();

  double Legendre;
  double skew = 0.0;
  double cov = 0.0;
  double sigma_sq = 0.0;
  double sigma_cu = 0.0;
  
  double noise_13 = 0.0;
  double noise_24 = 0.0;
  double noise_14 = 0.0;
  double noise_23 = 0.0;
  
  double n_random = 1.36/(constants::arcmin*constants::arcmin);
  
  if(field_index_1 == field_index_3){
    noise_13 = noise_1*noise_3;
  }
  else if(field_type_1 == GALAXY_DENSITY && field_type_3 == GALAXY_DENSITY && LSZ_correction_on){
    noise_13 = 1.0/n_random; // we use the same random points in all bins!
  }
  if(field_index_2 == field_index_4){
    noise_24 = noise_2*noise_4;
  }
  else if(field_type_2 == GALAXY_DENSITY && field_type_4 == GALAXY_DENSITY && LSZ_correction_on){
    noise_24 = 1.0/n_random; // we use the same random points in all bins!
  }
  if(field_index_1 == field_index_4){
    noise_14 = noise_1*noise_4;
  }
  else if(field_type_1 == GALAXY_DENSITY && field_type_4 == GALAXY_DENSITY && LSZ_correction_on){
    noise_14 = 1.0/n_random; // we use the same random points in all bins!
  }
  if(field_index_2 == field_index_3){
    noise_23 = noise_2*noise_3;
  }
  else if(field_type_2 == GALAXY_DENSITY && field_type_3 == GALAXY_DENSITY && LSZ_correction_on){
    noise_23 = 1.0/n_random; // we use the same random points in all bins!
  }
  
  
  for(int l = 0; l < n_ell-1; l++){
    Legendre =  (*Legendres)[i_th][l];
    sigma_sq = 0.0;
    if(l < ell_max-1){
      sigma_sq  += ((*Cell_13)[l] + noise_13)*((*Cell_24)[l] + noise_24);
      sigma_sq += ((*Cell_14)[l] + noise_14)*((*Cell_23)[l] + noise_23);
    }
    else{
      sigma_sq  += noise_13*noise_24;
      sigma_sq += noise_14*noise_23;
    }
    sigma_sq /= (2.0*double(l)+1.0);
    cov += pow((2.0*double(l)+1.0)*Legendre/(4.0*pi), 2)*sigma_sq;
  }
  
  cov /= f_sky;
  
  double C13, C24, C14, C23;
  
  for(int l = 0; l < n_ell-1; l++){
    Legendre =  (*Legendres)[i_th][l];
    sigma_cu = 0.0;
    if(l < ell_max-1){
    C13 = ((*Cell_13)[l] + noise_13);
    C24 = ((*Cell_24)[l] + noise_24);
    C14 = ((*Cell_14)[l] + noise_14);
    C23 = ((*Cell_23)[l] + noise_23);
    }
    else{
    C13 = (noise_13);
    C24 = (noise_24);
    C14 = (noise_14);
    C23 = (noise_23);
    }
    sigma_cu  += pow(C14, 3);
    sigma_cu  += 3.0*C14*C13*C24;
    sigma_cu *= 2.0/pow(2.0*double(l)+1.0, 2);
    skew += pow((2.0*double(l)+1.0)*Legendre/(4.0*pi), 3)*sigma_cu;
  }
  
  skew /= f_sky*f_sky;
  
  
  (*variance_value) = cov;
  (*skewness_value) = skew;
  
}




/* 
 * Variance of a field with power spectrum Cell, when averaged over a cicular patch of radius angular_radius.
 * 
 * 
 */
double return_variance(double angular_radius, vector<double> Cell){
  
  int n_ell = Cell.size();
  
  double variance = 0.0;
  double Legendre_factor;
  double prefactor;
  double cos_th = cos(angular_radius);
  double area = 2.0*pi*(1.0-cos_th);
  double area_sq = area*area;
  double ell;
  
  double Pl_th[n_ell+3]; 
    
  gsl_sf_legendre_Pl_array(n_ell+1, cos_th, Pl_th);

  for(int l = 1; l < n_ell-1; l++){
    ell = double(l);
    Legendre_factor = -(Pl_th[l+1]-Pl_th[l-1]);
    prefactor = pi/(area_sq*(2.0*ell +1.0));
    variance += Cell[l]*prefactor*Legendre_factor*Legendre_factor;
  }

  return variance;
  
}



/*
 * return_covariance:
 * Compute multi-probe 2-point covariance matrix.
 *
 * Params:
 * - double* covariance : pointer to output, see below
 * - const char *config_file_cstr : file containing various input quantities (see inline comments below)
 * - const char *Cell_files_cstr : file listing the file_names containing the various auto- and cross-power spectra
 * - double theta_min : minimal angular scale in arcmin
 * - double theta_max : maximal angular scale in arcmin
 * - int N_bins : number of log-spaced angular bins
 *
 * returns (in the form of pointers):
 * - double array covariance : array storing values of the covariance matrix
 *
 */
extern "C" void return_covariance(double* covariance, double* Gaussian_part, double *Cov_SS, double *Cov_SN, double *Cov_NN, const char *config_file_cstr, const char *Cell_files_cstr, double theta_min, double theta_max, int N_bin, double f_sky, int small_bin, int LSZ_correction_input, int Npair_modus_input){
  
  LSZ_correction_on = LSZ_correction_input;
  Npair_modus = Npair_modus_input;
  // The following inputs are read from config_file_cstr :
  int number_of_fields;
  int number_of_shear_fields = 0;
  int number_of_2pt_functions;
  vector<string> Cell_files(0);
  vector<TYPE_OF_FIELD> field_type(0);
  vector<double> number_density(0,0); // should be given in arcmin^-2
  vector<double> shape_noise(0,0); // set to 0 for galagy density field
  vector<double> log_normal_shift(0,0);
  
  /**** Reading input from configuration file. ****/
  fstream input_stream;
  string  input_string;
  int     input_int;
  double  input_double;
  input_stream.open(config_file_cstr);
  getline(input_stream, input_string); // first line is just head of table
  
  while(input_stream.good()){
    input_stream >> input_int;
    if(input_stream.good()){
      field_type.push_back(static_cast<TYPE_OF_FIELD>(input_int)); // 0 == galaxy density ; 1 == shear field
      if(input_int == 1) number_of_shear_fields++;
      
      input_stream >> input_double; number_density.push_back(input_double/constants::arcmin/constants::arcmin);
      
      input_stream >> input_double; shape_noise.push_back(input_double);
      
      input_stream >> input_double; log_normal_shift.push_back(input_double);
    }
  }
  input_stream.close();
  number_of_fields = field_type.size();
  
  input_stream.open(Cell_files_cstr);
  int i = 0;
  //cout << "Using the following C_ell files:\n";
  while(input_stream.good()){
    input_stream >> input_string;
    if(input_stream.good()){
      Cell_files.push_back(input_string);
      //cout << input_string << '\n';
    }
  }
  input_stream.close();
  number_of_2pt_functions = Cell_files.size();
  assert(number_of_2pt_functions == number_of_fields*(number_of_fields+1)/2);
  /**** Done reading input from configuration file. ****/
  
  /**** Reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  int ell_max = 0;
  int index = 0;
  vector<vector<vector<double> > > Cells(number_of_fields, vector<vector<double> >(number_of_fields, vector<double>(0, 0.0)));
  for(int i = 0; i < number_of_fields; i++){
    for(int j = i; j < number_of_fields; j++){
      Cells[i][j] = return_Cells_from_file(Cell_files[index]);
      /*if(field_type[i] == SHEAR_FIELD || field_type[j] == SHEAR_FIELD){
        int nn = Cells[i][j].size();
        Cells[i][j] = vector<double>(nn, 0.0);
      }*/
      Cells[j][i] = Cells[i][j];
      
      ell_max = max(ell_max, int(Cells[i][j].size()-1));
      index++;
    }
  }
  
  //cout << "ell_max = " << ell_max << '\n';
  //ell_max = min(ell_max, 8192);
  
  
  /**** Done reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  
  vector<double> bins = set_angular_bins_in_radians(theta_min, theta_max, N_bin);
  
  
  vector<vector<double> > Legendre_coefficients_wtheta;
  vector<vector<double> > Legendre_coefficients_gt;
  vector<vector<double> > Legendre_coefficients_xi_plus;
  vector<vector<double> > Legendre_coefficients_xi_minus;
  
  if(small_bin == 1){
    Legendre_coefficients_wtheta   = return_Legendres_small_bins(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_small_bins(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_small_bins(0, 1, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_small_bins(4, 1, ell_max, bins);
  }
  else if(small_bin == 2){ // This is flat sky!!!!
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins_flat_sky(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins_flat_sky(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins_flat_sky(0, 1, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins_flat_sky(4, 1, ell_max, bins);
  }
  else{
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins(0, 1, ell_max, bins);
    //Legendre_coefficients_xi_plus  = return_Legendres_finite_bins(0, 0, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins(4, 1, ell_max, bins);
  }
  
  vector<vector<double> > *pointer_to_Legendres_i;
  vector<vector<double> > *pointer_to_Legendres_j;
    
  //int n_data = (number_of_2pt_functions + number_of_shear_fields*(number_of_shear_fields+1)/2)*N_bin; // xi_minus makes everything more complicated...
  int n_data = (number_of_shear_fields*(number_of_shear_fields+1) + number_of_shear_fields*(number_of_fields - number_of_shear_fields) + (number_of_fields - number_of_shear_fields))*N_bin;
  int i_without_counting_2_cosmic_shear_2pt_functions = 0;
  int j_without_counting_2_cosmic_shear_2pt_functions = 0;
  int angular_index_i;
  int angular_index_j;
  int field_index_1;
  int field_index_2;
  int field_index_3;
  int field_index_4;
  
  /**** The array elements field_indeces_1[i] and field_indeces_2[i] store between which fields the i'th correlation function is formed. ****/
  index = 0;
  vector<int> field_indeces_1(n_data/N_bin, 0);
  vector<int> field_indeces_2(n_data/N_bin, 0);
  
  // Xi_plus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // Xi_minus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // gamma_t:
  for(int j = number_of_shear_fields; j < number_of_fields; j++){
    for(int i = 0; i < number_of_shear_fields; i++){
      field_indeces_1[index] = j;
      field_indeces_2[index] = i;
      index++;
    }
  }
  
  // wtheta:
  for(int i = number_of_shear_fields; i < number_of_fields; i++){
    field_indeces_1[index] = i;
    field_indeces_2[index] = i;
    index++;
  }
  /**** End of generating the arrays field_indeces_1 and field_indeces_2. ****/
  
  vector<double> noise_terms(number_of_fields, 0.0);
  double n_random = 1.36/constants::arcmin/constants::arcmin;
  for(int i = 0; i < number_of_fields; i++){
    switch(field_type[i]){
      case GALAXY_DENSITY:
        if(LSZ_correction_on)
          noise_terms[i] = sqrt(1.0/number_density[i]+1.0/n_random);
        else
          noise_terms[i] = sqrt(1.0/number_density[i]);
        break;
      case SHEAR_FIELD:
        noise_terms[i] = 1.0*shape_noise[i]/sqrt(number_density[i]);
        break;
    }
    
  }
  
  vector<double> data_vector(n_data, 0.0);
  for(int i = 0; i < n_data; i++){
    angular_index_i = i%N_bin;
    field_index_1 = field_indeces_1[i/N_bin];
    field_index_2 = field_indeces_2[i/N_bin];
    
    if(i < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_plus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_minus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_gt;
    }
    else{
      pointer_to_Legendres_i = &Legendre_coefficients_wtheta;
    }
    
    data_vector[i] = xi_curved_sky(angular_index_i, &Cells[field_index_1][field_index_2], pointer_to_Legendres_i);    
  }
  
  
  vector<vector<double> > lognormal_terms(number_of_fields, vector<double>(number_of_fields, 0.0));
  double survey_radius = sqrt(4.0*f_sky);
  for(int i = 0; i < number_of_fields; i++){
    for(int j = 0; j < number_of_fields; j++){
      lognormal_terms[i][j] = return_variance(survey_radius, Cells[i][j])/log_normal_shift[i]/log_normal_shift[j];
    }
  }
  
  
  index = 0;
  vector<double> mask_Cl = return_Cells_from_file("mask_cell_assuming_non-fractional_mask_and_Olivers_normalisation.dat");
  vector<double> mask_wtheta(N_bin, 0.0);
  vector<vector<double> > Legendre_coefficients_mask;
  if(small_bin == 0){
    Legendre_coefficients_mask = return_Legendres_finite_bins(0, 0, mask_Cl.size()-1, bins);
  }
  else if(small_bin == 1){ // This is flat sky!!!!
    Legendre_coefficients_mask = return_Legendres_finite_bins_flat_sky(0, 0, mask_Cl.size()-1, bins);
  }
  else if(small_bin == 2){ // This is flat sky!!!!
    Legendre_coefficients_mask = return_Legendres_finite_bins_flat_sky(0, 0, mask_Cl.size()-1, bins);
  }
  for(int i = 0; i < N_bin; i++){
    mask_wtheta[i] = xi_curved_sky(i, &mask_Cl, &Legendre_coefficients_mask);
  }

  
  for(int i = 0; i < n_data; i++){
    //cout << "\rcomputing row " << i+1 << "   ";
    //cout.flush();
    angular_index_i = i%N_bin;
    field_index_1 = field_indeces_1[i/N_bin];
    field_index_2 = field_indeces_2[i/N_bin];
    
    if(i < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_plus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_minus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_gt;
    }
    else{
      pointer_to_Legendres_i = &Legendre_coefficients_wtheta;
    }
    
    
    for(int j = 0; j < n_data; j++){
      angular_index_j = j%N_bin;
      field_index_3 = field_indeces_1[j/N_bin];
      field_index_4 = field_indeces_2[j/N_bin];
    
      if(j < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
        pointer_to_Legendres_j = &Legendre_coefficients_xi_plus;
      }
      else if(j < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
        pointer_to_Legendres_j = &Legendre_coefficients_xi_minus;
      }
      else if(j < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
        pointer_to_Legendres_j = &Legendre_coefficients_gt;
      }
      else{
        pointer_to_Legendres_j = &Legendre_coefficients_wtheta;
      }
      covariance_curved_sky(&covariance[index], &Gaussian_part[index], &Cov_SS[index], &Cov_SN[index], &Cov_NN[index], angular_index_i, angular_index_j, i, j, f_sky, noise_terms[field_index_1], noise_terms[field_index_2], noise_terms[field_index_3], noise_terms[field_index_4], field_index_1, field_index_2, field_index_3, field_index_4, field_type[field_index_1], field_type[field_index_2], field_type[field_index_3], field_type[field_index_4], &bins, &Cells[field_index_1][field_index_3], &Cells[field_index_2][field_index_4], &Cells[field_index_1][field_index_4], &Cells[field_index_2][field_index_3], pointer_to_Legendres_i, pointer_to_Legendres_j, &data_vector, &lognormal_terms, &mask_wtheta);
            
      index++;
    }
  }
  //cout << "\rcomputation done            \n";
  
  
}



/*
 * return_data_vector:
 * Compute multi-probe 2-point covariance matrix.
 *
 * Params:
 * - double* data_vector : pointer to output, see below
 * - const char *config_file_cstr : file containing various input quantities (see inline comments below)
 * - const char *Cell_files_cstr : file listing the file_names containing the various auto- and cross-power spectra
 * - double theta_min : minimal angular scale in arcmin
 * - double theta_max : maximal angular scale in arcmin
 * - int N_bins : number of log-spaced angular bins
 *
 * returns (in the form of pointers):
 * - double array data_vector : array storing values of the data vector
 *
 */
extern "C" void return_data_vector(double* data_vector, const char *config_file_cstr, const char *Cell_files_cstr, double theta_min, double theta_max, int N_bin, double f_sky, int small_bin, int LSZ_correction_input){
  
  LSZ_correction_on = LSZ_correction_input;
  // The following inputs are read from config_file_cstr :
  int number_of_fields;
  int number_of_shear_fields = 0;
  int number_of_2pt_functions;
  vector<string> Cell_files(0);
  vector<TYPE_OF_FIELD> field_type(0);
  vector<double> number_density(0,0); // should be given in arcmin^-2
  vector<double> shape_noise(0,0); // set to 0 for galagy density field
  vector<double> log_normal_shift(0,0);
  
  /**** Reading input from configuration file. ****/
  fstream input_stream;
  string  input_string;
  int     input_int;
  double  input_double;
  input_stream.open(config_file_cstr);
  getline(input_stream, input_string); // first line is just head of table
  
  while(input_stream.good()){
    input_stream >> input_int;
    if(input_stream.good()){
      field_type.push_back(static_cast<TYPE_OF_FIELD>(input_int)); // 0 == galaxy density ; 1 == shear field
      if(input_int == 1) number_of_shear_fields++;
      
      input_stream >> input_double; number_density.push_back(input_double/constants::arcmin/constants::arcmin);
      
      input_stream >> input_double; shape_noise.push_back(input_double);
      
      input_stream >> input_double; log_normal_shift.push_back(input_double);
    }
  }
  input_stream.close();
  number_of_fields = field_type.size();
  
  input_stream.open(Cell_files_cstr);
  int i = 0;
  while(input_stream.good()){
    input_stream >> input_string;
    if(input_stream.good()){
      Cell_files.push_back(input_string);
      //cout << input_string << '\n';
    }
  }
  input_stream.close();
  number_of_2pt_functions = Cell_files.size();
  assert(number_of_2pt_functions == number_of_fields*(number_of_fields+1)/2);
  /**** Done reading input from configuration file. ****/
  
  /**** Reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  int ell_max = 0;
  int index = 0;
  vector<vector<vector<double> > > Cells(number_of_fields, vector<vector<double> >(number_of_fields, vector<double>(0, 0.0)));
  for(int i = 0; i < number_of_fields; i++){
    for(int j = i; j < number_of_fields; j++){
      Cells[i][j] = return_Cells_from_file(Cell_files[index]);
      Cells[j][i] = Cells[i][j];
      ell_max = max(ell_max, int(Cells[i][j].size()-1));
      index++;
    }
  }
  
  ell_max = min(ell_max, 8192);
  //ell_max = min(ell_max, 19000);
  
  
  /**** Done reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  
  vector<double> bins = set_angular_bins_in_radians(theta_min, theta_max, N_bin);
  
  
  vector<vector<double> > Legendre_coefficients_wtheta;
  vector<vector<double> > Legendre_coefficients_gt;
  vector<vector<double> > Legendre_coefficients_xi_plus;
  vector<vector<double> > Legendre_coefficients_xi_minus;
  
  if(small_bin == 1){
    Legendre_coefficients_wtheta   = return_Legendres_small_bins(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_small_bins(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_small_bins(0, 1, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_small_bins(4, 1, ell_max, bins);
  }
  else if(small_bin == 2){ // This is flat sky!!!!
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins_flat_sky(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins_flat_sky(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins_flat_sky(0, 1, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins_flat_sky(4, 1, ell_max, bins);
  }
  else{
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins(2, 0, ell_max, bins);
    //Legendre_coefficients_xi_plus  = return_Legendres_finite_bins(0, 1, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins(0, 0, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins(4, 1, ell_max, bins);
  }
  
  vector<vector<double> > *pointer_to_Legendres_i;
  vector<vector<double> > *pointer_to_Legendres_j;
  
  
  //int n_data = (number_of_2pt_functions + number_of_shear_fields*(number_of_shear_fields+1)/2)*N_bin; // xi_minus makes everything more complicated...
  int n_data = (number_of_shear_fields*(number_of_shear_fields+1) + number_of_shear_fields*(number_of_fields - number_of_shear_fields) + (number_of_fields - number_of_shear_fields))*N_bin;
  int i_without_counting_2_cosmic_shear_2pt_functions = 0;
  int j_without_counting_2_cosmic_shear_2pt_functions = 0;
  int angular_index_i;
  int angular_index_j;
  int field_index_1;
  int field_index_2;
  int field_index_3;
  int field_index_4;
  
  /**** The array elements field_indeces_1[i] and field_indeces_2[i] store between which fields the i'th correlation function is formed. ****/
  index = 0;
  vector<int> field_indeces_1(n_data/N_bin, 0);
  vector<int> field_indeces_2(n_data/N_bin, 0);
  
  // Xi_plus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // Xi_minus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // gamma_t:
  for(int j = number_of_shear_fields; j < number_of_fields; j++){
    for(int i = 0; i < number_of_shear_fields; i++){
      field_indeces_1[index] = j;
      field_indeces_2[index] = i;
      index++;
    }
  }
  
  // wtheta:
  for(int i = number_of_shear_fields; i < number_of_fields; i++){
    field_indeces_1[index] = i;
    field_indeces_2[index] = i;
    index++;
  }
  /**** End of generating the arrays field_indeces_1 and field_indeces_2. ****/
  
  for(int i = 0; i < n_data; i++){
    angular_index_i = i%N_bin;
    field_index_1 = field_indeces_1[i/N_bin];
    field_index_2 = field_indeces_2[i/N_bin];
    
    if(i < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_plus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_minus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_gt;
    }
    else{
      pointer_to_Legendres_i = &Legendre_coefficients_wtheta;
    }
    
    data_vector[i] = xi_curved_sky(angular_index_i, &Cells[field_index_1][field_index_2], pointer_to_Legendres_i);    
  }
  
}







extern "C" void return_pair_counts(double* pair_counts, double* mask_wtheta_values, const char *config_file_cstr, const char *Cell_files_cstr, double theta_min, double theta_max, int N_bin, double f_sky, int l_max, int small_bin){
  
  // The following inputs are read from config_file_cstr :
  int number_of_fields;
  int number_of_shear_fields = 0;
  int number_of_2pt_functions;
  vector<string> Cell_files(0);
  vector<TYPE_OF_FIELD> field_type(0);
  vector<double> number_density(0,0); // should be given in arcmin^-2
  vector<double> shape_noise(0,0); // set to 0 for galagy density field
  vector<double> log_normal_shift(0,0);
  
  /**** Reading input from configuration file. ****/
  fstream input_stream;
  string  input_string;
  int     input_int;
  double  input_double;
  input_stream.open(config_file_cstr);
  getline(input_stream, input_string); // first line is just head of table
  
  while(input_stream.good()){
    input_stream >> input_int;
    if(input_stream.good()){
      field_type.push_back(static_cast<TYPE_OF_FIELD>(input_int)); // 0 == galaxy density ; 1 == shear field
      if(input_int == 1) number_of_shear_fields++;
      
      input_stream >> input_double; number_density.push_back(input_double/constants::arcmin/constants::arcmin);
      
      input_stream >> input_double; shape_noise.push_back(input_double);
      
      input_stream >> input_double; log_normal_shift.push_back(input_double);
    }
  }
  input_stream.close();
  number_of_fields = field_type.size();
  
  input_stream.open(Cell_files_cstr);
  int i = 0;
  while(input_stream.good()){
    input_stream >> input_string;
    if(input_stream.good()){
      Cell_files.push_back(input_string);
    }
  }
  input_stream.close();
  number_of_2pt_functions = Cell_files.size();
  assert(number_of_2pt_functions == number_of_fields*(number_of_fields+1)/2);
  /**** Done reading input from configuration file. ****/
  
  
  /**** Reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  int ell_max = 0;
  int index = 0;
  vector<vector<vector<double> > > Cells(number_of_fields, vector<vector<double> >(number_of_fields, vector<double>(0, 0.0)));
  for(int i = 0; i < number_of_fields; i++){
    for(int j = i; j < number_of_fields; j++){
      Cells[i][j] = return_Cells_from_file(Cell_files[index]);
      Cells[j][i] = Cells[i][j];
      ell_max = max(ell_max, int(Cells[i][j].size()-1));
      index++;
    }
  }
  
  ell_max = min(l_max, ell_max);
  
  /**** Done reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  
  vector<double> bins = set_angular_bins_in_radians(theta_min, theta_max, N_bin);
  
  
  vector<vector<double> > Legendre_coefficients_wtheta;
  vector<vector<double> > Legendre_coefficients_gt;
  vector<vector<double> > Legendre_coefficients_xi_plus;
  vector<vector<double> > Legendre_coefficients_xi_minus;
  
  if(small_bin == 1){
    Legendre_coefficients_wtheta   = return_Legendres_small_bins(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_small_bins(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_small_bins(0, 1, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_small_bins(4, 1, ell_max, bins);
  }
  else if(small_bin == 2){ // This is flat sky!!!!
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins_flat_sky(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins_flat_sky(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins_flat_sky(0, 1, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins_flat_sky(4, 1, ell_max, bins);
  }
  else{
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins(0, 0, ell_max, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins(2, 0, ell_max, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins(0, 1, ell_max, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins(4, 1, ell_max, bins);
  }
  
  
  vector<vector<double> > *pointer_to_Legendres_i;
  vector<vector<double> > *pointer_to_Legendres_j;
  
  
  int n_data = (number_of_shear_fields*(number_of_shear_fields+1) + number_of_shear_fields*(number_of_fields - number_of_shear_fields) + (number_of_fields - number_of_shear_fields))*N_bin;
  int i_without_counting_2_cosmic_shear_2pt_functions = 0;
  int j_without_counting_2_cosmic_shear_2pt_functions = 0;
  int angular_index_i;
  int angular_index_j;
  int field_index_1;
  int field_index_2;
  int field_index_3;
  int field_index_4;
  
  /**** The array elements field_indeces_1[i] and field_indeces_2[i] store between which fields the i'th correlation function is formed. ****/
  index = 0;
  vector<int> field_indeces_1(n_data/N_bin, 0);
  vector<int> field_indeces_2(n_data/N_bin, 0);
  
  // Xi_plus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // Xi_minus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // gamma_t:
  for(int j = number_of_shear_fields; j < number_of_fields; j++){
    for(int i = 0; i < number_of_shear_fields; i++){
      field_indeces_1[index] = j;
      field_indeces_2[index] = i;
      index++;
    }
  }
  
  // wtheta:
  for(int i = number_of_shear_fields; i < number_of_fields; i++){
    field_indeces_1[index] = i;
    field_indeces_2[index] = i;
    index++;
  }
  /**** End of generating the arrays field_indeces_1 and field_indeces_2. ****/
  
  
  
  //vector<double> mask_Cl = return_Cells_from_file("mask_cell_assuming_non-fractional_mask.dat");
  vector<double> mask_Cl = return_Cells_from_file("../shape_noise_tests/Cl_footprint_DES_y3_03_31_30.txt");
  vector<double> mask_wtheta(N_bin, 0.0);
  vector<vector<double> > Legendre_coefficients_mask;
  if(small_bin == 0){
    Legendre_coefficients_mask = return_Legendres_finite_bins(0, 0, mask_Cl.size()-1, bins);
  }
  if(small_bin == 1){
    Legendre_coefficients_mask = return_Legendres_small_bins(0, 0, mask_Cl.size()-1, bins);
  }
  else if(small_bin == 2){ // This is flat sky!!!!
    Legendre_coefficients_mask = return_Legendres_finite_bins_flat_sky(0, 0, mask_Cl.size()-1, bins);
  }
  vector<vector<double> > Legendre_coefficients_mask_small_bin   = return_Legendres_small_bins(0, 0, mask_Cl.size()-1, bins);
  for(int i = 0; i < N_bin; i++){
    mask_wtheta[i] = xi_curved_sky(i, &mask_Cl, &Legendre_coefficients_mask)-1.0;
  }
  
  vector<double> data_vector_dummy(n_data, 0.0);
  vector<double> Gaussian_dummy(n_data, 0.0);
  
  for(int i = 0; i < n_data; i++){
    angular_index_i = i%N_bin;
    field_index_1 = field_indeces_1[i/N_bin];
    field_index_2 = field_indeces_2[i/N_bin];
    
    if(i < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_plus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_minus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_gt;
    }
    else{
      pointer_to_Legendres_i = &Legendre_coefficients_wtheta;
    }
    
    mask_wtheta_values[i] = mask_wtheta[i%N_bin];
    pair_counts[i] = (1.0+mask_wtheta[i%N_bin])*number_density[field_index_1]*number_density[field_index_2]*f_sky*4.0*constants::pi*2.0*constants::pi*(cos(bins[i%N_bin])-cos(bins[i%N_bin+1]));
    if(i>=800){
      pair_counts[i] /= 2.0;//*(1.0+xi_curved_sky(angular_index_i, &Cells[field_index_1][field_index_2], pointer_to_Legendres_i));
    }
  }
  
  
  
}

/*
 * return_covariance:
 * Compute multi-probe 2-point covariance matrix.
 *
 * Params:
 * - double* covariance : pointer to output, see below
 * - const char *config_file_cstr : file containing various input quantities (see inline comments below)
 * - const char *Cell_files_cstr : file listing the file_names containing the various auto- and cross-power spectra
 * - double theta_min : minimal angular scale in arcmin
 * - double theta_max : maximal angular scale in arcmin
 * - int N_bins : number of log-spaced angular bins
 *
 * returns (in the form of pointers):
 * - double array covariance : array storing values of the covariance matrix
 *
 */
extern "C" void return_skewness_and_variance(double* variance_values, double* skewness_values, const char *config_file_cstr, const char *Cell_files_cstr, double theta_min, double theta_max, int N_bin, double f_sky, int small_bin, int LSZ_correction_input, int Npair_modus_input){
    
  LSZ_correction_on = LSZ_correction_input;
  Npair_modus = Npair_modus_input;
  // The following inputs are read from config_file_cstr :
  int number_of_fields;
  int number_of_shear_fields = 0;
  int number_of_2pt_functions;
  vector<string> Cell_files(0);
  vector<TYPE_OF_FIELD> field_type(0);
  vector<double> number_density(0,0); // should be given in arcmin^-2
  vector<double> shape_noise(0,0); // set to 0 for galagy density field
  vector<double> log_normal_shift(0,0);
  
  /**** Reading input from configuration file. ****/
  fstream input_stream;
  string  input_string;
  int     input_int;
  double  input_double;
  input_stream.open(config_file_cstr);
  getline(input_stream, input_string); // first line is just head of table
  
  while(input_stream.good()){
    input_stream >> input_int;
    if(input_stream.good()){
      field_type.push_back(static_cast<TYPE_OF_FIELD>(input_int)); // 0 == galaxy density ; 1 == shear field
      if(input_int == 1) number_of_shear_fields++;
      
      input_stream >> input_double; number_density.push_back(input_double/constants::arcmin/constants::arcmin);
      
      input_stream >> input_double; shape_noise.push_back(input_double);
      
      input_stream >> input_double; log_normal_shift.push_back(input_double);
    }
  }
  input_stream.close();
  number_of_fields = field_type.size();
  
  input_stream.open(Cell_files_cstr);
  int i = 0;
  while(input_stream.good()){
    input_stream >> input_string;
    if(input_stream.good()){
      Cell_files.push_back(input_string);
      //cout << input_string << '\n';
    }
  }
  input_stream.close();
  number_of_2pt_functions = Cell_files.size();
  assert(number_of_2pt_functions == number_of_fields*(number_of_fields+1)/2);
  /**** Done reading input from configuration file. ****/
  
  /**** Reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  int ell_max = 0;
  int index = 0;
  vector<vector<vector<double> > > Cells(number_of_fields, vector<vector<double> >(number_of_fields, vector<double>(0, 0.0)));
  for(int i = 0; i < number_of_fields; i++){
    for(int j = i; j < number_of_fields; j++){
      Cells[i][j] = return_Cells_from_file(Cell_files[index]);
      Cells[j][i] = Cells[i][j];
      ell_max = max(ell_max, int(Cells[i][j].size()-1));
      index++;
    }
  }
  
  //cout << "ell_max = " << ell_max << '\n';
  //ell_max = min(ell_max, 8192);
  //ell_max = min(ell_max, 19000);
  
  
  /**** Done reading Cell_files and determining the maximum ell-mode appearing in the input Cells. ****/
  
  vector<double> bins = set_angular_bins_in_radians(theta_min, theta_max, N_bin);
  
  
  vector<vector<double> > Legendre_coefficients_wtheta;
  vector<vector<double> > Legendre_coefficients_gt;
  vector<vector<double> > Legendre_coefficients_xi_plus;
  vector<vector<double> > Legendre_coefficients_xi_minus;
  
  int ell_max_Legendre = 10000;
  ell_max_Legendre = ell_max;
  
  if(small_bin == 1){
    Legendre_coefficients_wtheta   = return_Legendres_small_bins(0, 0, ell_max_Legendre, bins);
    Legendre_coefficients_gt       = return_Legendres_small_bins(2, 0, ell_max_Legendre, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_small_bins(0, 1, ell_max_Legendre, bins);
    Legendre_coefficients_xi_minus = return_Legendres_small_bins(4, 1, ell_max_Legendre, bins);
  }
  else if(small_bin == 2){ // This is flat sky!!!!
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins_flat_sky(0, 0, ell_max_Legendre, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins_flat_sky(2, 0, ell_max_Legendre, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins_flat_sky(0, 1, ell_max_Legendre, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins_flat_sky(4, 1, ell_max_Legendre, bins);
  }
  else{
    Legendre_coefficients_wtheta   = return_Legendres_finite_bins(0, 0, ell_max_Legendre, bins);
    Legendre_coefficients_gt       = return_Legendres_finite_bins(2, 0, ell_max_Legendre, bins);
    Legendre_coefficients_xi_plus  = return_Legendres_finite_bins(0, 1, ell_max_Legendre, bins);
    Legendre_coefficients_xi_minus = return_Legendres_finite_bins(4, 1, ell_max_Legendre, bins);
  }
  
  vector<vector<double> > *pointer_to_Legendres_i;
  vector<vector<double> > *pointer_to_Legendres_j;
  
  
  //int n_data = (number_of_2pt_functions + number_of_shear_fields*(number_of_shear_fields+1)/2)*N_bin; // xi_minus makes everything more complicated...
  int n_data = (number_of_shear_fields*(number_of_shear_fields+1) + number_of_shear_fields*(number_of_fields - number_of_shear_fields) + (number_of_fields - number_of_shear_fields))*N_bin;
  int i_without_counting_2_cosmic_shear_2pt_functions = 0;
  int j_without_counting_2_cosmic_shear_2pt_functions = 0;
  int angular_index_i;
  int angular_index_j;
  int field_index_1;
  int field_index_2;
  int field_index_3;
  int field_index_4;
  
  /**** The array elements field_indeces_1[i] and field_indeces_2[i] store between which fields the i'th correlation function is formed. ****/
  index = 0;
  vector<int> field_indeces_1(n_data/N_bin, 0);
  vector<int> field_indeces_2(n_data/N_bin, 0);
  
  // Xi_plus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // Xi_minus:
  for(int i = 0; i < number_of_shear_fields; i++){
    for(int j = i; j < number_of_shear_fields; j++){
      field_indeces_1[index] = i;
      field_indeces_2[index] = j;
      index++;
    }
  }
  
  // gamma_t:
  for(int j = number_of_shear_fields; j < number_of_fields; j++){
    for(int i = 0; i < number_of_shear_fields; i++){
      field_indeces_1[index] = j;
      field_indeces_2[index] = i;
      index++;
    }
  }
  
  // wtheta:
  for(int i = number_of_shear_fields; i < number_of_fields; i++){
    field_indeces_1[index] = i;
    field_indeces_2[index] = i;
    index++;
  }
  /**** End of generating the arrays field_indeces_1 and field_indeces_2. ****/
    
  vector<double> noise_terms(number_of_fields, 0.0);
  double n_random = 1.36/constants::arcmin/constants::arcmin;
  for(int i = 0; i < number_of_fields; i++){
    switch(field_type[i]){
      case GALAXY_DENSITY:
        if(LSZ_correction_on)
          noise_terms[i] = sqrt(1.0/number_density[i]+1.0/n_random);
        else
          noise_terms[i] = sqrt(1.0/number_density[i]);
        break;
      case SHEAR_FIELD:
        noise_terms[i] = 1.0*shape_noise[i]/sqrt(number_density[i]);
        break;
    }
    
  }
  
  vector<double> data_vector(n_data, 0.0);
  for(int i = 0; i < n_data; i++){
    angular_index_i = i%N_bin;
    field_index_1 = field_indeces_1[i/N_bin];
    field_index_2 = field_indeces_2[i/N_bin];
    
    if(i < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_plus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_minus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_gt;
    }
    else{
      pointer_to_Legendres_i = &Legendre_coefficients_wtheta;
    }
    
    data_vector[i] = xi_curved_sky(angular_index_i, &Cells[field_index_1][field_index_2], pointer_to_Legendres_i);    
  }
  
  
  vector<vector<double> > lognormal_terms(number_of_fields, vector<double>(number_of_fields, 0.0));
  double survey_radius = sqrt(4.0*f_sky);
  for(int i = 0; i < number_of_fields; i++){
    for(int j = 0; j < number_of_fields; j++){
      lognormal_terms[i][j] = return_variance(survey_radius, Cells[i][j])/log_normal_shift[i]/log_normal_shift[j];
    }
  }
  
  
  index = 0;
  
  //vector<double> mask_Cl = return_Cells_from_file("mask_cell_assuming_non-fractional_mask.dat");
  vector<double> mask_Cl = return_Cells_from_file("mask_cell_assuming_circular_mask_and_my_normalisation_and_Tims_area.dat");
  //vector<double> mask_Cl = return_Cells_from_file("Cl_footprint_DES_y1.dat");
  vector<double> mask_wtheta(N_bin, 0.0);
  vector<vector<double> > Legendre_coefficients_mask;
  if(small_bin == 0){
    Legendre_coefficients_mask = return_Legendres_finite_bins(0, 0, mask_Cl.size()-1, bins);
  }
  else if(small_bin == 1){ // This is flat sky!!!!
    Legendre_coefficients_mask = return_Legendres_finite_bins_flat_sky(0, 0, mask_Cl.size()-1, bins);
  }
  else if(small_bin == 2){ // This is flat sky!!!!
    Legendre_coefficients_mask = return_Legendres_finite_bins_flat_sky(0, 0, mask_Cl.size()-1, bins);
  }
  vector<vector<double> > Legendre_coefficients_mask_small_bin   = return_Legendres_small_bins(0, 0, mask_Cl.size()-1, bins);
  for(int i = 0; i < N_bin; i++){
    mask_wtheta[i] = xi_curved_sky(i, &mask_Cl, &Legendre_coefficients_mask);
    //cout << i << "  " << mask_wtheta[i] << '\n';
  }


  
  for(int i = 0; i < n_data; i++){
    //cout << i << '\n';
    angular_index_i = i%N_bin;
    field_index_1 = field_indeces_1[i/N_bin];
    field_index_2 = field_indeces_2[i/N_bin];
    
    if(i < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_plus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_xi_minus;
    }
    else if(i < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
      pointer_to_Legendres_i = &Legendre_coefficients_gt;
    }
    else{
      pointer_to_Legendres_i = &Legendre_coefficients_wtheta;
    }
    
    
    for(int j = 0; j < n_data; j++){
      angular_index_j = j%N_bin;
      field_index_3 = field_indeces_1[j/N_bin];
      field_index_4 = field_indeces_2[j/N_bin];
    
      if(j < number_of_shear_fields*(number_of_shear_fields+1)/2*N_bin){
        pointer_to_Legendres_j = &Legendre_coefficients_xi_plus;
      }
      else if(j < number_of_shear_fields*(number_of_shear_fields+1)*N_bin){
        pointer_to_Legendres_j = &Legendre_coefficients_xi_minus;
      }
      else if(j < number_of_shear_fields*(number_of_shear_fields+1)*N_bin + number_of_shear_fields*(number_of_fields - number_of_shear_fields)*N_bin){
        pointer_to_Legendres_j = &Legendre_coefficients_gt;
      }
      else{
        pointer_to_Legendres_j = &Legendre_coefficients_wtheta;
      }
      
      if(i == j){
        variance_and_skewness(&variance_values[i], &skewness_values[i], angular_index_i, f_sky, noise_terms[field_index_1], noise_terms[field_index_2], noise_terms[field_index_3], noise_terms[field_index_4], field_index_1, field_index_2, field_index_3, field_index_4, field_type[field_index_1], field_type[field_index_2], field_type[field_index_3], field_type[field_index_4], &bins, &Cells[field_index_1][field_index_3], &Cells[field_index_2][field_index_4], &Cells[field_index_1][field_index_4], &Cells[field_index_2][field_index_3], pointer_to_Legendres_i, &data_vector, &lognormal_terms, &mask_wtheta);
      }
            
      index++;
    }
  }
  
  
}
