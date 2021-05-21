#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include "ConfigFile.h"

using namespace std;

double const pi=3.14159265358979323846264338327950288419716939937510582097494459230e0;

typedef vector<complex<double> > vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// ou A est une matrice tridiagonale
template <class T> void triangular_solve(vector<T> const& diag,
                                         vector<T> const& lower,
                                         vector<T> const& upper,
                                         vector<T> const& rhs,
                                         vector<T>& solution)
{
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;

  // forward elimination
  for(unsigned int i(1); i<diag.size(); ++i)
  {
    T pivot = lower[i-1] / new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution.resize(diag.size());

  // solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  // backward substitution
  for(int i = int(diag.size()) - 2; i >= 0; --i)
  {
    solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
  }
}


// Potentiel V(x) : // TODO
double V(double const& x, double const& omega2, double const& Delta)
{
  return 0.5*omega2*min(pow(x-Delta,2),pow(x+Delta,2)); //masse pas déclarée à ce stade
}

// Declaration des diagnostics de la particule d'apres sa fonction d'onde psi :
//  - prob calcule la probabilite de trouver la particule entre les points de maillage nL et nR,
//  - E calcule son energie moyenne,
//  - xmoy calcule sa position moyenne,
//  - x2moy calcule sa position au carre moyenne,
//  - pmoy calcule sa quantite de mouvement moyenne,
//  - p2moy calcule sa quantite de mouvement au carre moyenne.
double prob(vec_cmplx const& psi, int nL, int nR, double dx);
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx);
double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double pmoy(vec_cmplx const& psi, double const& dx);
double p2moy(vec_cmplx const& psi, double const& dx);

// Fonction pour normaliser une fonction d'onde :
vec_cmplx normalize(vec_cmplx const& psi, double const& dx);

// Les definitions de ces fonctions sont en dessous du main.


int main(int argc,char **argv)
{
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres physiques :
  double tfin           = configFile.get<double>("tfin");
  double xL             = configFile.get<double>("xL");
  double xR             = configFile.get<double>("xR");
  double xda            = configFile.get<double>("xda");
  double xdb            = configFile.get<double>("xdb");
  double hbar           = configFile.get<double>("hbar");;
  double m              = configFile.get<double>("mass");
  double omega          = configFile.get<double>("omega");
  double Delta          = configFile.get<double>("Delta");
  double x0             = configFile.get<double>("x0");
  double k0             = 2. * pi * double(configFile.get<int>("n")) / (xR-xL);
  double sigma0         = configFile.get<double>("sigma_norm") * (xR-xL);
  double t_detect       = configFile.get<double>("t_detect");

  double omega2 = m*omega*omega;

  // Parametres numeriques :
  double dt      = configFile.get<double>("dt");
  int Ninters    = configFile.get<int>("Ninters");
  int Npoints    = Ninters + 1;
  double dx      = (xR-xL) / Ninters;

  // Maillage :
  vector<double> x(Npoints);
  for(int i(0); i<Npoints; ++i)
    x[i] = xL + i*dx;

  // Initialisation de la fonction d'onde :
  vec_cmplx psi(Npoints);
  // TODO: initialiser le paquet d'onde, equation (4.116) du cours
  for(int i(0); i<Npoints; ++i)
    psi[i] = /*1/sqrt(sigma0*sqrt(pi)) **/ exp(complex_i*k0*x[i])*exp( -pow(x[i]-x0,2)/(2*pow(sigma0,2)) ); // MODIFY
  // Modifications des valeurs aux bords :
  psi[0] = 0;
  psi.back() = 0;
  // Normalisation :
  psi = normalize(psi, dx);

  // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
  vec_cmplx dH(Npoints), aH(Ninters), cH(Ninters); // matrice Hamiltonienne
  vec_cmplx dA(Npoints), aA(Ninters), cA(Ninters); // matrice du membre de gauche de l'equation (4.99)
  vec_cmplx dB(Npoints), aB(Ninters), cB(Ninters); // matrice du membre de droite de l'equation (4.99)

  complex<double> a(0,hbar*dt/(4*m*pow(dx,2))); // Coefficient complexe a, Eq.(4.100)

  // TODO: calculer les elements des matrices A, B et H.
  // Ces matrices sont stockees sous forme tridiagonale, d:diagonale, c et a: diagonales superieures et inferieures
  for(int i(0); i<Npoints; ++i) // Boucle sur les points de maillage
  { //diagonale
    complex<double> b(0,dt*V(x[i],omega2,Delta)/(2*hbar));
    dH[i] = pow(hbar,2)/(m*pow(dx,2)) + V(x[i],omega2,Delta);
    dA[i] = 1.0+2.0*a+b;
    dB[i] = 1.0-2.0*a-b;
  }
  for(int i(0); i<Ninters; ++i) // Boucle sur les intervalles
  { //sous-diagonale et sur-diagonale
    aH[i] = -pow(hbar,2)/(2*m*pow(dx,2));
    cH[i] = -pow(hbar,2)/(2*m*pow(dx,2));
    aA[i] = -a;
    cA[i] = -a;
    aB[i] =  a;
    cB[i] =  a;
  }

  // Conditions aux limites: psi nulle aux deux bords
  // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
  dA[0] = 1; dA.back() = 1;
  dB[0] = 0; dB.back() = 0;

  cA[0] = 0; /*cA.back() = 0;*/
  /*aA[0] = 0;*/ aA.back() = 0;

  cB[0] = 0; /*cB.back() = 0;*/
  /*aB[0] = 0*/; aB.back() = 0;

  // Fichiers de sortie :
  ofstream fichier_potentiel((configFile.get<string>("output_potential")).c_str());
  fichier_potentiel.precision(15);
  for(int i(0); i<Npoints; ++i)
    fichier_potentiel << x[i] << " " << V(x[i], omega2, Delta) << endl;
  fichier_potentiel.close();

  ofstream fichier_psi((configFile.get<string>("output_squared_wave")).c_str());
  fichier_psi.precision(15);

  ofstream fichier_observables((configFile.get<string>("output_observables")).c_str());
  fichier_observables.precision(15);

  // ecrire position en x
  fichier_psi << 0.0 << " ";
  for(int i(0); i<Npoints-1; ++i){
    fichier_psi << x[i] << " ";
  }
  fichier_psi << x[Npoints-1] << endl;

  // Boucle temporelle :
  valarray<double> print_array=valarray<double>(0.e0,Npoints+1);
  double t,window;
  for(t=0.; t+dt/*/2.*/<=tfin; t+=dt)
  {
    // Detection de la particule  
    if(round(t/dt) == round(t_detect/dt))
    {

      for(int i(0); i<abs(Ninters*xL/(xL-xR)); ++i){
        psi[i] = complex<double> (0.,0.);
      }
      for(int i(abs(Ninters*xL/(xL-xR))); i<abs(Ninters*(xda-xL)/(xL-xR)); ++i){
        window = pow(sin(0.5*pi*x[i]/xda),2.e0);
        psi[i] = polar(window*abs(psi[i]),arg(psi[i]));
      }
      for(int i(abs(Ninters*(xdb-xL)/(xL-xR))); i<Ninters; ++i){
        window = pow(0.5*pi*x[i]/(xR-xdb),2.0);
        psi[i] = polar(window*abs(psi[i]),arg(psi[i]));
      }
      psi = normalize(psi, dx); // normalise psi pour que la proba totale soit 1
    }

    // Ecriture de |psi|^2 :
    print_array[0] = t;
    for(int i(1); i<Npoints+1; ++i){
      print_array[i] = norm(psi[i-1]); // la fonction C++ norm prend le module au carre
    }
    for(int i(0); i<Npoints; ++i){
      fichier_psi << print_array[i] << " ";
    }
    fichier_psi << print_array[Npoints] << endl;

    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi,0,abs(Ninters*xL/(xL-xR)),dx) << " "          // probabilite que la particule soit en x < 0
                        << prob(psi,abs(Ninters*xL/(xL-xR)),Ninters,dx) << " "    // probabilite que la particule soit en x > 0
                        << prob(psi,0,Ninters,dx) << " " 		   	                  // probabilite totale
                        << E(psi,dH,aH,cH,dx) << " "                       	      // Energie
                        << xmoy(psi,x,dx) << " "                           	      // Position moyenne
                        << x2moy(psi,x,dx) << " "                          	      // Position^2 moyenne
                        << pmoy(psi,dx) << " "                             	      // Quantite de mouvement moyenne
                        << p2moy(psi,dx) << " "                           	      // (Quantite de mouvement)^2 moyenne
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx))*\
			   sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << " "       // Heisenberg index
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx)) << " " // incertitude en x
                        << sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << endl;     // incertitude en p

    // Calcul du membre de droite :
    vec_cmplx psi_tmp(Npoints,0.);

    // Multiplication psi_tmp = B * psi :
    for(int i(0); i<Npoints; ++i)
      psi_tmp[i] = dB[i] * psi[i];
    for(int i(0); i<Ninters; ++i)
    {
      psi_tmp[i] += cB[i] * psi[i+1];
      psi_tmp[i+1] += aB[i] * psi[i];
    }

    // Resolution de A * psi = psi_tmp :
    triangular_solve(dA, aA, cA, psi_tmp, psi);
    //if(t+dt > tfin) dt = tfin - t;
  } // Fin de la boucle temporelle

    // ecrire |psi|^2
    print_array[0] = t;
    for(int i(0); i<Npoints; ++i){
      print_array[i+1] = norm(psi[i]);
    }
    for(int i(0); i<Npoints; ++i){
      fichier_psi << print_array[i] << " ";
    }
    fichier_psi << print_array[Npoints] << endl;

    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi,0,abs(Ninters*xL/(xL-xR)),dx) << " "              // probabilite que la particule soit en x < 0
                        << prob(psi,abs(Ninters*xL/(xL-xR)),Ninters,dx) << " "        // probabilite que la particule soit en x > 0
                        << prob(psi,0,Ninters,dx) << " " 		   	      // probabilite totale
                        << E(psi,dH,aH,cH,dx) << " "                       	      // Energie
                        << xmoy(psi,x,dx) << " "                           	      // Position moyenne
                        << x2moy(psi,x,dx) << " "                          	      // Position^2 moyenne
                        << pmoy(psi,dx) << " "                             	      // Quantite de mouvement moyenne
                        << p2moy(psi,dx) << " "                           	      // (Quantite de mouvement)^2 moyenne
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx))*\
			   sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << " "       // Heisenberg index
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx)) << " " // incertitude en x
                        << sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << endl;     // incertitude en p

  fichier_observables.close();
  fichier_psi.close();

}

double prob(vec_cmplx const& psi, int nL, int nR, double dx)
{
  //TODO: calculer la probabilite de trouver la particule entre les points nL et nR
  //Utiliser la formule des trapezes pour l'integration
  double resultat(0.);
  for(int i(nL); i < nR-1 ; ++i)  resultat += norm(psi[i+1]) + norm(psi[i]);
  return resultat * dx * 0.5;
}

//TODO: Calculer les valeurs moyennes des observables E, x, p, x^2, p^2
//Utiliser la formule des trapezes pour l'integration
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx)
{
  vec_cmplx psi_tmp(psi.size()); // vecteur pour stocker H*psi
  double resultat(0.); // initialiser

  // H(psi): produit de la matrice H et du  vecteur psi
  for(unsigned int i(1); i<diagH.size()-1; ++i)
    psi_tmp[i] = lowerH[i-1]*psi[i-1] + diagH[i]*psi[i] + upperH[i]*psi[i+1];
  psi_tmp[0] = diagH[0]*psi[0] + upperH[0]*psi[1];
  psi_tmp.back() = lowerH.back()*psi[diagH.size()-2] + diagH.back()*psi.back();

  // Integrale de psi* H(psi) dx
  for(size_t i(0); i < psi.size()-1 ; ++i)
    resultat += real(conj(psi[i+1])*psi_tmp[i+1] + conj(psi[i])*psi_tmp[i]); //on prend la partie réelle pour éviter double = complex
  return 0.5*dx*resultat;
}

double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  double resultat(0.);
  for(size_t i(0); i < psi.size()-1 ; ++i)
    resultat += real(conj(psi[i+1])*x[i+1]*psi[i+1] + conj(psi[i])*x[i]*psi[i]);
  return 0.5*dx*resultat;
}

double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  double resultat(0.);
  for(size_t i(0); i < psi.size()-1 ; ++i)
    resultat += real(conj(psi[i+1])*pow(x[i+1],2)*psi[i+1] + conj(psi[i])*pow(x[i],2)*psi[i]);
  return 0.5*dx*resultat;
}

double pmoy(vec_cmplx const& psi, double const& dx) //hbar pas inclu dans les arguments et pas défini !!
{
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i
  //unsigned int N(psi.size());
  // Utiliser la definition de p = -i hbar d/dx
  // Utiliser les differences finies centrees pour d/dx
  // Utiliser la formule des trapezes pour l'integration sur x
  // Ignorer la contribution du premier et du dernier point de maillage
  vec_cmplx lambda(psi.size());
  for(unsigned int i(1); i<lambda.size()-1; ++i)
    lambda[i] = 0.5/dx * (psi[i+1]-psi[i-1]);
  lambda[0] = (psi[1] - psi[0])/dx;
  lambda.back() = (psi.back() - psi[psi.size()-2])/dx;

  double resultat(0.);
  for(size_t i(0); i < psi.size()-2 ; ++i)
    resultat += real( -complex_i*(conj(psi[i+1])*lambda[i+1] + conj(psi[i])*lambda[i]) ) ;
  return 0.5*dx*resultat;
}

double p2moy(vec_cmplx const& psi, double const& dx) //hbar pas défini !!
{
  double resultat(0.);
  // Utiliser la definition de p^2 = hbar^2 d^2/dx2
  // Utiliser les differences finies centrees pour d^2/dx^2
  // Utiliser la formule des trapezes pour l'integration sur x
  // Ignorer la contribution du premier et du dernier point de maillage

  vec_cmplx lambda(psi.size());
  for(unsigned int i(1); i<lambda.size()-1; ++i)
    lambda[i] = pow(dx,-2) * (psi[i+1]-complex<double>(2,0)*psi[i] + psi[i-1]);
  lambda[0] = 0;
  lambda.back() = 0;

  for(size_t i(0); i < psi.size()-1 ; ++i)
    resultat += real((conj(psi[i+1])*lambda[i+1] + conj(psi[i])*lambda[i])) ;
  return -0.5*dx*resultat;
}

vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{
  vec_cmplx psi_norm(psi.size());
  double norm = sqrt(prob(psi,0,psi.size()-1,dx));
  for(unsigned int i(0); i<psi.size(); ++i)
    psi_norm[i] = psi[i]/norm;
  return psi_norm;
}
