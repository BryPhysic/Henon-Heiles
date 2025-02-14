#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>

// Parámetros globales
double E;       // Energía total
double T_MAX;   // Tiempo total de integración
double DT;      // Paso de integración


double potential(double x, double y) {
    
    return 0.5 * (x*x + y*y) + 1.0*(x*x*y - (1.0/3.0) * y*y*y);
}

// ----------------------------------------------------------------------
// Gradiente del potencial: dV/dx, dV/dy
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
// Integrador Leapfrog (Stormer-Verlet)
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
// Calcula px(0) para energía E y posición inicial (x=0, y=y0, py=0):
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
// Programa Principal
// ======================================================================
int main() {
    // Lectura de parámetros
    std::cout << "E (energia)   : ";
    std::cin >> E;
    std::cout << "T_MAX (tiempo): ";
    std::cin >> T_MAX;
    std::cout << "DT (paso)     : ";
    std::cin >> DT;

    // Crear archivo de salida
    std::ostringstream fname;
    fname << "HHE_" << T_MAX << "_" << DT << "_" << E << ".dat";
    std::ofstream fout(fname.str().c_str());
    if (!fout.is_open()) {
        std::cerr << "No se pudo abrir " << fname.str() << std::endl;
        return 1;
    }
    fout << std::fixed << std::setprecision(6);

    // Rango de valores iniciales de y0
    std::vector<double> y_iniciales;
    for (double y0 = -0.5; y0 <= 0.5; y0 += 0.01) {
        y_iniciales.push_back(y0);
    }

    // Número de pasos de integración
    long long Nsteps = static_cast<long long>(T_MAX / DT);

    // Integrar para cada condición inicial en y0
    int y0_index = 0;
    for (double y0 : y_iniciales) {
        double px0 = initial_px(y0);
        if (px0 < 0.0) {
            continue; // No hay solución real
        }

        // Condiciones iniciales
        double x  = 0.0;
        double y  = y0;
        double px = px0;
        double py = 0.0;
        std::cout << "Integrando con y0=" << y0 
                  << ", px0=" << px0 
                  << ", E=" << E << "\n";
        // Integración
        for (long long step = 0; step < Nsteps; step++) {
            double x_antes  = x;
            double y_antes  = y;
            double px_antes = px;
            double py_antes = py;

            // Integrar con Leapfrog
            leapfrog_step(x, y, px, py, DT);

            double x_despues  = x;
            double y_despues  = y;
            double px_despues = px;
            double py_despues = py;

            // Detección de cruce en x=0 con px>0
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
