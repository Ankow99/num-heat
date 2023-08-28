
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to set initial conditions
void setInitialConditions(double *T, int N, double T0, double Tm, double TL, double L) {
    double dx = L / (N - 1);
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        T[i] = Tm * sin(M_PI * x / L) + (((TL - T0) / L) * x) + T0;
    }
}

// Thomas algorithm for solving tridiagonal systems
void thomasAlgorithm(double *a, double *b, double *c, double *d, double *T, int N) {
    // Forward elimination
    for (int i = 1; i < N; i++) {
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    // Back substitution
    T[N - 1] = d[N - 1] / b[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        T[i] = (d[i] - c[i] * T[i + 1]) / b[i];
    }
}

int main() {
    // Parameters
    double L = 1.0;
    double alpha = 0.1;
    double Tm = 300.0;
    double T0 = 300.0;
    double TL = 1000.0;
    int numTimeSteps = 2000;
    int numSpatialSteps = 100;

    // Calculate dx
    double dx = L / (numSpatialSteps - 1);

    // Calculate dt based on the stability condition
    double dt = 2.0 * dx * dx / alpha;

    // Calculate r
    double r = alpha * dt / (dx * dx);

    // Initialize gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    fprintf(gnuplotPipe, "set title 'Heat Equation Solution Comparison'\n");
    fprintf(gnuplotPipe, "set xlabel 'Position (x)'\n");
    fprintf(gnuplotPipe, "set ylabel 'Temperature (T)'\n");

    // Allocate memory 
    double *T = (double *)malloc(numSpatialSteps * sizeof(double));
    double *a = (double *)malloc(numSpatialSteps * sizeof(double));
    double *b = (double *)malloc(numSpatialSteps * sizeof(double));
    double *c = (double *)malloc(numSpatialSteps * sizeof(double));
    double *d = (double *)malloc(numSpatialSteps * sizeof(double));
    double *errors = (double *)malloc(numTimeSteps * sizeof(double));

    if (T == NULL || a == NULL || b == NULL || c == NULL || d == NULL || errors == NULL) {
        perror("Memory allocation failed");
        return 1; // Exit with an error
    }

    // Set initial conditions
    setInitialConditions(T, numSpatialSteps, T0, Tm, TL, L);

    // Time-stepping loop
    for (int j = 1; j <= numTimeSteps; j++) {

        // Set coefficients for the tridiagonal system
        for (int i = 0; i < numSpatialSteps; i++) {
            a[i] = -1;
            b[i] = 2 * ((1 / r) + 1);
            c[i] = -1;
            d[i] = T[i-1] + 2 * ((1 / r) - 1) * T[i] + T[i+1];

            // Apply boundary conditions
            if (i == 0) {
                d[i] = T0;
                b[i] = 1;
                c[i] = 0;
            } else if (i == numSpatialSteps - 1) {
                d[i] = TL;
                a[i] = 0;
                b[i] = 1;
            }
        }

        // Solve tridiagonal system using Thomas algorithm
        thomasAlgorithm(a, b, c, d, T, numSpatialSteps);

        // Plot numerical and analytical solutions
        fprintf(gnuplotPipe, "plot '-' with lines title 'Numerical Time Step %d', '-' with lines title 'Analytical'\n", j);

        // Plot numerical solution
        for (int i = 0; i < numSpatialSteps; i++) {
            fprintf(gnuplotPipe, "%f %f\n", i * dx, T[i]);
        }
        fprintf(gnuplotPipe, "e\n");

        // Calculate error
        double error = 0.0;
        
        // Plot analytical solution
        for (int i = 0; i < numSpatialSteps; i++) {
            double x = i * dx;
            double t = j * dt;
            double analytical = Tm * sin(M_PI * x / L) * exp(-alpha * M_PI * M_PI * t / (L * L)) + (((TL - T0) / L) * x) + T0;
            fprintf(gnuplotPipe, "%f %f\n", x, analytical);
            error += fabs(analytical - T[i]);
        }
        fprintf(gnuplotPipe, "e\n");

        error /= numSpatialSteps;
        errors[j - 1] = error;
        
        fflush(gnuplotPipe);
    }

    // Plot the error over time steps
    FILE *errorGnuplotPipe = popen("gnuplot -persist", "w");
    fprintf(errorGnuplotPipe, "set title 'Error Comparison'\n");
    fprintf(errorGnuplotPipe, "set xlabel 'Time Step'\n");
    fprintf(errorGnuplotPipe, "set ylabel 'Average Error'\n");
    fprintf(errorGnuplotPipe, "plot '-' with lines title 'Error'\n");
    for (int j = 0; j < numTimeSteps; j++) {
        fprintf(errorGnuplotPipe, "%d %f\n", j + 1, errors[j]);
    }
    fprintf(errorGnuplotPipe, "e\n");

    fflush(errorGnuplotPipe);

    // Close gnuplot
    pclose(gnuplotPipe);
    pclose(errorGnuplotPipe);

    // Free memory
    free(T);
    free(a);
    free(b);
    free(c);
    free(d);
    free(errors);

    return 0;
}
