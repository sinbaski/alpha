#!/usr/bin/python

import numpy as np
from numpy import real as Re, imag as Im, sin, cos, exp, sqrt, nan, inf, pi
from scipy.fft import fft, ifft, fftfreq
from scipy.special import j0, jvp, gammainc
from scipy import integrate
import sqlite3
import matplotlib.pyplot as plt

kelvin_offset = 273.15;

def cylinderTemperature(params, r, t):
    r0 = params["r0"];
    h = params['h'];
    L = params['L'];
    k = params['kappa'];
    tau = params['tau'];
    L = params['L'];
    T = params['T'];
    
    b = sqrt(25/k[0]/t) if t > 0 else 512;
    radial = integrate.quad(
        lambda l: 2*l * j0(r0 * l) * j0(r * l) / b**2 / (j0(r0 * b)**2 + jvp(0, r0 * b)**2) * exp(-l**2 * k[0] * t),
        0, b, limit = 100
    )[0];

    dz = h/400;
    Z = np.linspace(0, h, int(h/dz));
    F = fft(gammainc(0.5, (Z + L)**2/k[1]/tau) * sqrt(pi));
    Omega = 2*pi * fftfreq(len(F), dz);
    vertical = Re(
        ifft(exp(-Omega**2 * k[0] * t) * F)
    );
    u = -(T[1] - T[0])/sqrt(pi) * radial * vertical + T[1];
    return u;

def barTemperature(params, x, t):
    k = params['kappa'];
    T = params['T'];
    A = T[1] - T[0];
    u = -A * gammainc(0.5, x**2/k/t) + T[1];
    return u;

if __name__ == "__main__":
    material = 'W';
    conn = sqlite3.connect("../materials.db");
    cs = conn.cursor();
    cs.execute(F"""
    select thermal_conductivity/specific_heat/density from Materials
    where formula = '{material}';
    """);
    k = cs.fetchall()[0][0];
    electrode_params = {
        'T': [kelvin_offset + 25, 3100],
        'kappa': k
    };
    L = 10.0e-2;
    tau = 1;
    T = barTemperature(electrode_params, L, tau);
    print(T - kelvin_offset);


    fml = 'Al2O3';
    r = 8.0e-2;
    cs.execute(F"""
    select formula, thermal_conductivity/specific_heat/density
    from Materials
    where formula = '{fml}';
    """);
    kappa = {data[0]: data[1] for data in cs.fetchall()};
    cs.close();
    conn.close();

    shell_params = {
        'r0': 5.0e-3/2,
        'h': 1.0e-2,
        'L': L,
        'kappa': [kappa[fml], k],
        'tau': tau,
        'T': [kelvin_offset + 25, 3100]
    };
    dr = r/100;

    R = np.arange(shell_params['r0'], r, dr);
    T = [cylinderTemperature(shell_params, x, 1.0)[0] for x in R];
    plt.plot(R, T);
    plt.grid();
    plt.show();
    
    # cs.execute(F"""
    # select formula, thermal_conductivity/specific_heat/density
    # from Materials
    # where formula != 'W'
    # and max_temperature > {T - kelvin_offset};
    # """);
    # kappa = {data[0]: data[1] for data in cs.fetchall()};
    # for fml in kappa:
    #     shell_params = {
    #         'r0': 5.0e-3/2,
    #         'h': 1.0e-2,
    #         'L': L,
    #         'kappa': [kappa[fml], k],
    #         'tau': tau,
    #         'T': [kelvin_offset + 25, 3100]
    #     };
    #     T = cylinderTemperature(shell_params, r, 0, 0.01);
    #     print(fml, T - kelvin_offset);
    
    
