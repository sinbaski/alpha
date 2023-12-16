#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <gsl/gsl_math.h>

// #define N_MIX 3
// #define N_VAR 4
// #define MP_IDX 0
// #define MG_IDX 1
// #define MC_IDX 2
// #define I_IDX 2
// #define T_IDX 3
using namespace std;

const double T0 = 2300;

enum StateVariable {
    Plastic, Graphene, CarbonBlack, n_ingredients,
    Current=2, Temperature, n_statvar
};

const double resistivity[] = {
    pow(1.0e+13, 0.5) * pow(1.0e+16, 0.5), // HDPE, LDPE - High/Low density Polyethylene
    pow(1.0e-7, 0.5) * pow(1.0e-2, 0.5), // graphene nano-platelets:
    0.1 // carbon black
};

const double density[] = {
    961.0, //HDPE: 970 kg/m^3
    0.065 * 1.0e-3 * 1.0e+6, // graphene: 0.065 g/cc
    1.79 * 1.0e-3 * 1.0e+6 // carbon black: 1.79 g/cc
};

struct measurements {
    double radius;
    double height;
    double V0;
    double capacity;
    double Cmass;
    double ambient_temperature;
};


double calc_resistivity(const array<double, n_ingredients> &masses)
{
    array<double, n_ingredients> V;
    double V_tot = 0;
    double rho = 1;
    
    #pragma omp parallel for
    for (size_t i = 0; i < n_ingredients; i++) {
        V[i] = masses[i]/density[i];
        V_tot += V[i];
    }
    #pragma omp parallel for
    for (double &x: V) {
        x /= V_tot;
    }
    for (size_t i = 0; i < n_ingredients; i++) {
        rho *= pow(resistivity[i], V[i]);
    }
    return rho;
}

double calc_graphene_ratio(double T)
{
    static const double bbar = 0.8571428571428571;
    static const double b1 = 0.55;
    static const double T1 = 3100;
    static const double alpha = -log(bbar - b1)/(T1 - T0);
    return bbar - exp(-alpha * (T - T0));
}

void calc_time_deriv(const struct measurements &msmt, const array<double, n_statvar> &statvar,
                     double r, array<double, n_statvar> &deriv)
{
    deriv[Plastic] = -r * (statvar[Temperature] - T0) * statvar[Plastic];
    deriv[Graphene] = -deriv[Plastic] * calc_graphene_ratio(statvar[Temperature]);
    static double A = M_PI * gsl_pow_2(msmt.radius);
    double R = calc_resistivity(statvar) * msmt.height / A;
    
}

void simulate_pulse(const struct measurements &msmt, array<double, n_statvar> &statvar)
{
    double A = M_PI * gsl_pow_2(msmt.radius);
    array<double, n_ingredients> masses = {
        statvar[Plastic], statvar[Graphene], msmt.Cmass
    };
    double rho = calc_resistivity(masses);
    statvar[Current] = msmt.V0/(rho * msmt.height / A);
    
    // const double dt = 1.0e-4;
    while (statvar[Current] > 0) {
        if (statvar[Temperature] < T0) {
        } else {
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc < 4) {
        printf("Usage: jh.exe m_P m_C\n");
        return 0;
    }
    
    double mp0 = strtod(argv[1], NULL);
    double carbon_black_ratio = strtod(argv[2], NULL);

    struct measurements msmt = {
        .radius = 4.3e-2,
        .height = 10.5e-2,
        .V0 = 220,
        .capacity = 60e-3,
        .Cmass = carbon_black_ratio * mp0,
        .ambient_temperature = 273.15 + 25
    };
    // The masses of plastic, graphene and carbon black
    array<double, n_ingredients> masses = {
        mp0, 0, mp0 * carbon_black_ratio
    };
    // m_P, I, T
    array<double, n_statvar> state_var;
    state_var[Plastic] = masses[Plastic];
    state_var[Graphene] = masses[Graphene];
    state_var[CarbonBlack] = masses[CarbonBlack];
    state_var[Temperature] = msmt.ambient_temperature;
    
    // unsigned int n_iter = 10000;
    
    
    // printf("Plastic:\t%.4f\n", masses[0]);
    // printf("Graphene:\t%.4f\n", masses[1]);
    // printf("Carbon black:\t%.4f\n", masses[2]);
    return 0;
}
