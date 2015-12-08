#include <iostream>
#include <cmath>

using namespace std;

// vektorielle DGL
void vecDiv(double* F, double t, double* vec, double eta){
  F[0] = vec[1];
  F[1] = (eta - vec[0]*vec[0]) * vec[0];
}

// vektorielles Runke-Kutta-Verfahren mit
//      0   |   0    0    0
//      0.5 |   0.5  0    0
//      1   |  -1    2    0
//      ----------------------
//          |   0.16 0.67 0.16
// als Butcher-Tableau
void vecRuKu(void(*vecDiv)(double*, double, double*, double),
	     double dis, double t, double* vec, double eta){
  double k[2] = {0.,0.};
  double tmp[2]; // Kopie der urspruenglichen Werte (nicht ueberschr.)
  double tmpk[2] = {0., 0.}; // temporaere Variable fuer k
  for(int i = 0; i < 2; i++){
    tmp[i] = vec[i];
  }
  
  vecDiv(k, t, tmp, eta); // k1
  for(int i = 0; i < 2; i++){
    tmp[i] = tmp[i] + 0.5 * dis * k[i];
    tmpk[i] = k[i];
    vec[i] += dis/6.0 * k[i];
  }
  
  vecDiv(k, t, tmp, eta); // k2
  for(int i = 0; i < 2; i++){
    tmp[i] = tmp[i] + 2 * dis * k[i] - tmpk[i] * dis;
    vec[i] += 2 * dis/3.0 * k[i];
  }
  
  vecDiv(k, t, tmp, eta); // k3
  for(int i = 0; i < 2; i++){
    vec[i] += dis/6.0 * k[i];
  }
}

int main(void){
  const double eta = 10.0; // ETA-Parameter
  double t = 0.0; // Start: 0.0, Ende: 100.0
  
  int N; int ind = 0;
  //cout << "Bitte Schrittanzahl eingeben:" << endl;
  cin >> N;

  const double dt = 100.0/N; // Schrittlaenge
  
  //cout << dt << endl;
  
  double Yps[2];//[N]; // Yps = {[x][y][z]}
  Yps[0] = 5.; Yps[1] = sqrt(eta)*Yps[0]; // initiale Werte x(0), y(0) und z(0)
  
  cout << "# t \t x \t y \t z" << endl;
  cout << t << "\t" << Yps[0] << "\t" << Yps[1] << endl; // Ausgabe Y(0)
  
  while(t < 100.0){
    vecRuKu(vecDiv, dt, t, Yps, eta); // Runge-Kutta
    cout << t << "\t" << Yps[0] << "\t" << Yps[1] << endl; // Ausgabe Yn, Y'n
    t += dt;
    ind++;
  }

  return 0;
}