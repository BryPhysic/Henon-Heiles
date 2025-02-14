#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>


double E;       // 
double T_MAX;   // 
double DT;      // 


double potential(double x, double y) {
    
    return 0.5 * (x*x + y*y) + 1.0*(x*x*y - (1.0/3.0) * y*y*y);
}

// ----------------------------------------------------------------------
// 
// ----------------------------------------------------------------------
void gradV(double x, double y, double &dVdx, double &dVdy) {
    dVdx = x + 2.0 * x * y;
    dVdy = y + x*x - y*y;
}

double hamiltonian(double x, double y, double px, double py) {
    double kinetic = 0.5 * (px*px + py*py);
    return kinetic + potential(x, y);
}

// ----------------------------------------------------------------------
// 
// ----------------------------------------------------------------------
void leapfrog_step(double &x, double &y, double &px, double &py, double dt) {
    double dVdx, dVdy;
    gradV(x, y, dVdx, dVdy);

    double half_dt = 0.5 * dt;

    // 1) Medio paso en el momento
    px -= half_dt * dVdx;
    py -= half_dt * dVdy;

    // 2) Paso completo en coordenadas
    x += dt * px;
    y += dt * py;

    // 3) Otro medio paso en el momento
    gradV(x, y, dVdx, dVdy);
    px -= half_dt * dVdx;
    py -= half_dt * dVdy;
}

// ----------------------------------------------------------------------
// 
// E = 0.5 px^2 + V(0,y0)  =>  px^2 = 2[E - V(0,y0)]
// ----------------------------------------------------------------------
double initial_px(double y0) {
    double V0 = potential(0.0, y0);
    double arg = 2.0 * (E - V0);
    if (arg < 0.0) {
        return -1.0; // No hay solución real
    }
    return std::sqrt(arg);
}

// ======================================================================
// 
// ======================================================================
int main() {
    // Lectura de parámetros
    std::cout << "E (energia)   : ";
    std::cin >> E;
    std::cout << "T_MAX (tiempo): ";
    std::cin >> T_MAX;
    std::cout << "DT (paso)     : ";
    std::cin >> DT;

    // 
    std::ostringstream fname;
    fname << "HHE_" << T_MAX << "_" << DT << "_" << E << ".dat";
    std::ofstream fout(fname.str().c_str());
    if (!fout.is_open()) {
        std::cerr << "No se pudo abrir " << fname.str() << std::endl;
        return 1;
    }
    fout << std::fixed << std::setprecision(6);

    // 
    std::vector<double> y_iniciales;
    for (double y0 = -0.5; y0 <= 0.5; y0 += 0.01) {
        y_iniciales.push_back(y0);
    }

    //
    long long Nsteps = static_cast<long long>(T_MAX / DT);

    // 
    int y0_index = 0;
    for (double y0 : y_iniciales) {
        double px0 = initial_px(y0);
        if (px0 < 0.0) {
            continue; // 
        }

        
        double x  = 0.0;
        double y  = y0;
        double px = px0;
        double py = 0.0;
        std::cout << "Integrando con y0=" << y0 
                  << ", px0=" << px0 
                  << ", E=" << E << "\n";
       
        for (long long step = 0; step < Nsteps; step++) {
            double x_antes  = x;
            double y_antes  = y;
            double px_antes = px;
            double py_antes = py;

            
            leapfrog_step(x, y, px, py, DT);

            double x_despues  = x;
            double y_despues  = y;
            double px_despues = px;
            double py_despues = py;

 
            if (x_antes < 0.0 && x_despues > 0.0 && px_despues > 0.0) {
                double alpha_t = -x_antes / (x_despues - x_antes);
                double y_int   = y_antes  + alpha_t * (y_despues - y_antes);
                double py_int  = py_antes + alpha_t * (py_despues - py_antes);

                fout << y_int << " " << py_int << " " << y0_index <<  "\n";
            }
        }
        y0_index++;
    }

    fout.close();
    std::cout << "Datos guardados en " << fname.str() << "\n";

    return 0;
}
