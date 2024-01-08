#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <cmath>
#include <array>
#include <vector>
#include <memory>
#include <numeric>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_num.h>
#include "DebyeModel.hpp"

using namespace std;

shared_ptr<FILE> outfile;

const double T0 = 2300;
const double T1 = 3100;
const double destruction_factor = -2*log(0.01)/(T1 - T0);

enum StateVariable {
    Plastic, Graphene, CarbonBlack, n_ingredients,
    ElectricCharge=2, Temperature, n_statvar
};

const double resistivity[] = {
    // HDPE, LDPE - High/Low density Polyethylene
    pow(1.0e+13, 0.5) * pow(1.0e+16, 0.5),
    // graphene nano-platelets:
    pow(1.0e-7, 0.5) * pow(1.0e-2, 0.5),
    0.1 // carbon black
};

const double density[] = {
    961.0, //HDPE: 970 kg/m^3
    2.3 * 1.0e-3 * 1.0e+6, // graphene: 2.3 g/cc
    1.79 * 1.0e-3 * 1.0e+6 // carbon black: 1.79 g/cc
};

//kg/mol
const double molar_mass[] = {
    6.27e+3, //HDPE, average
    12.011e-3, //Carbon molar mass
    12.011e-3 //Carbon Black
};

const double Debye_temperature[] = {
    NAN, // HDPE, average
    1813, // Graphene
    1813, // Carbon black
};

const double specific_heat[] = {
    1.82e+3, //HDPE, average
    NAN, //Carbon molar mass
    NAN //Carbon Black
};


struct measurements {
    double radius;
    double height;
    double V0;
    double capacity;
    double Cmass;
    double ambient_temperature;
};


void calc_volume_ratios(const double *masses, double *V)
{
    double V_tot = 0;    
    for (size_t i = 0; i < n_ingredients; i++) {
        V[i] = masses[i]/density[i];
        V_tot += V[i];
    }
    #pragma omp parallel for
    for (size_t i = 0; i < n_ingredients; i++) {
        V[i] /= V_tot;
    }
}
    
double calc_conductivity(const double *masses, const double *V_ratios)
{
    double sigma = 0;
    
    for (size_t i = 0; i < n_ingredients; i++) {
        sigma += V_ratios[i] / resistivity[i];
    }
    return sigma;
}

double calc_graphene_ratio(double T)
{
    static const double bbar = 0.8571428571428571;
    static const double b1 = bbar * 0.99;
    // static const double alpha = -log(bbar - b1)/(T1 - T0);
    static const double alpha = -log(1 - b1/bbar)/(T1 - T0);
    return bbar * (1 - exp(-alpha * (T - T0)));
}

double calc_heat_capacity(int mtrl, double mass, double T)
{
    double Cv;
    double N;
    if (mass == 0) return 0;
    
    switch (mtrl) {
    case Plastic:
	Cv = specific_heat[Plastic] * mass;
	break;
    case Graphene:
    case CarbonBlack:
	N = GSL_CONST_NUM_AVOGADRO * mass/molar_mass[mtrl];
	Cv = Debye_heat_capacity(T, N, Debye_temperature[mtrl]);
	break;
    default:
	Cv = NAN;
    }
    return Cv;
}

void calc_time_deriv(const struct measurements &msmt,
		     const array<double, n_statvar> &statvar,
                     double r, array<double, n_statvar> &deriv)
{
    static double A = M_PI * gsl_pow_2(msmt.radius);
    double masses[n_ingredients], V[n_ingredients];
    memcpy(masses, statvar.data(), (n_ingredients - 1) * sizeof(double));
    masses[CarbonBlack] = msmt.Cmass;
    calc_volume_ratios(masses, V);

    deriv[Plastic] = -r * (statvar[Temperature] - T0) * statvar[Plastic];
    deriv[Graphene] = -deriv[Plastic] * calc_graphene_ratio(statvar[Temperature]);
    
    double R;
    double mixture_conductivity = calc_conductivity(masses, V);
    R = msmt.height / A / mixture_conductivity;

    double I = statvar[ElectricCharge]/msmt.capacity/R;
    deriv[ElectricCharge] = -I;

    double Cv = 0;
    for (size_t i = 0; i < n_ingredients; i++) {
	Cv += calc_heat_capacity(i, masses[i], statvar[Temperature]);
    }
    deriv[Temperature] = gsl_pow_2(I) * R / Cv;
}

void print_statvar(double t, const array<double, n_statvar> &statvar)
{
    if (! outfile) {
	outfile.reset(fopen("./trajectory.csv", "w"), fclose);
    }
    fprintf(outfile.get(), "%.4f, ", t);
    for (const auto &x: statvar) {
	fprintf(outfile.get(), "%.4f, ", x);
    }
    fprintf(outfile.get(), "\n");
    fflush(outfile.get());
}

void simulate_pulse(const struct measurements &msmt, array<double, n_statvar> &statvar, double &t)
{
    double ratio = msmt.height / (M_PI * gsl_pow_2(msmt.radius));
    double dt = 1.0e-3;
    double Q = statvar[ElectricCharge];
    
    while (statvar[ElectricCharge]/Q > 0.01) {
	print_statvar(t, statvar);
	
        if (statvar[Temperature] < T0) {
	    double masses[n_ingredients], V[n_ingredients];
	    memcpy(masses, statvar.data(), sizeof(double) * (n_ingredients - 1));
	    masses[CarbonBlack] = msmt.Cmass;
	    calc_volume_ratios(masses, V);
	    double R = ratio / calc_conductivity(masses, V);
	    double I = statvar[ElectricCharge]/msmt.capacity/R;
	    double dE = gsl_pow_2(I) * R * dt;
	    double Cv = 0;
	    for (size_t i = 0; i < n_ingredients; i++) {
		Cv += calc_heat_capacity(i, masses[i], statvar[Temperature]);
	    }
	    statvar[Temperature] += dE/Cv;
	    statvar[ElectricCharge] -= I * dt;
        } else {
	    array<double, n_statvar> deriv;
	    calc_time_deriv(msmt, statvar, destruction_factor, deriv);
	    auto statvar_old = statvar;
	    bool negflag = false;
	    for (size_t i = 0; i < n_statvar; i++) {
		statvar[i] += deriv[i] * dt;
		negflag = negflag || statvar[i] < 0;
	    }
	    if (statvar[Plastic] < 0 && statvar[ElectricCharge] > 0) {
		dt = statvar_old[Plastic]/fabs(deriv[Plastic]);
	    } else if (statvar[Plastic] > 0 && statvar[ElectricCharge] < 0) {
		dt = statvar_old[ElectricCharge]/fabs(deriv[ElectricCharge]);
	    } else if (statvar[Plastic] < 0 && statvar[ElectricCharge] < 0) {
		double dt1 = statvar_old[Plastic]/fabs(deriv[Plastic]);
		dt = statvar_old[ElectricCharge]/fabs(deriv[ElectricCharge]);
		dt = min(dt, dt1);
	    }
	    if (negflag) {
		for (size_t i = 0; i < n_statvar; i++) {
		    statvar[i] = statvar_old[i] + deriv[i] * dt;
		}
	    }
        }
	t += dt;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 4) {
        printf("Usage: jh.exe m_P r_C height\n");
        return 0;
    }
    double mp0 = strtod(argv[1], NULL);
    double carbon_black_ratio = strtod(argv[2], NULL);
    double mc = carbon_black_ratio * mp0;
    double init_vol = mp0/density[Plastic] + mc/density[CarbonBlack];
    
    struct measurements msmt = {
        .radius = NAN,
        .height = strtod(argv[3], NULL),
        .V0 = 1000,
        .capacity = 1,
        .Cmass = mc,
        .ambient_temperature = 273.15 + 25
    };
    msmt.radius = sqrt(init_vol/msmt.height/M_PI);
    // The masses of plastic, graphene and carbon black
    double masses[] = {
        mp0, 0, mp0 * carbon_black_ratio
    };
    
    array<double, n_statvar> state_var;
    state_var[Plastic] = masses[Plastic];
    state_var[Graphene] = masses[Graphene];
    state_var[Temperature] = msmt.ambient_temperature;
    state_var[ElectricCharge] = msmt.V0 * msmt.capacity;

    double t = 0;
    do {
	simulate_pulse(msmt, state_var, t);
	msmt.V0 = 220;
	state_var[ElectricCharge] = msmt.capacity * msmt.V0;
	t += 1.0;
    } while (state_var[Plastic] / masses[Plastic] > 0.02);

    /* while (state_var[Plastic] / masses[Plastic] > 0.05) { */
    /* 	simulate_pulse(msmt, state_var); */
    /* } */
    return 0;
}
