/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "Models.h"


namespace MaterialLib
{

double SalineWaterDensityIAPWS_IF97(double Press, double TempK, double Conc)
{
    /*int c_idx;*/
    double rho_0;
    double GammaPi, Pressurevar, Tau, pressure_average, temperature_average;
    double Tstar, Pstar,GazConst;
    double L[35],J[35],n[35];
    int i;
    double salinity;

    pressure_average = Press;
    temperature_average = TempK;
    salinity = Conc;
    Tstar = 1386;
    Pstar = 16.53e6;                      // MPa
    GazConst = 0.461526e3;                //

    n[0] = 0.0;
    n[1] = 0.14632971213167;
    n[2] = -0.84548187169114;
    n[3] = -0.37563603672040e1;
    n[4] = 0.33855169168385e1;
    n[5] = -0.95791963387872;
    n[6] = 0.15772038513228;
    n[7] = -0.16616417199501e-1;
    n[8] = 0.81214629983568e-3;
    n[9] = 0.28319080123804e-3;
    n[10] = -0.60706301565874e-3;
    n[11] = -0.18990068218419e-1;
    n[12] = -0.32529748770505e-1;
    n[13] = -0.21841717175414e-1;
    n[14] = -0.52838357969930e-4;
    n[15] = -0.47184321073267e-3;
    n[16] = -0.30001780793026e-3;
    n[17] = 0.47661393906987e-4;
    n[18] = -0.44141845330846e-5;
    n[19] = -0.72694996297594e-15;
    n[20] = -0.31679644845054e-4;
    n[21] = -0.28270797985312e-5;
    n[22] = -0.85205128120103e-9;
    n[23] = -0.22425281908000e-5;
    n[24] = -0.65171222895601e-6;
    n[25] = -0.14341729937924e-12;
    n[26] = -0.40516996860117e-6;
    n[27] = -0.12734301741641e-8;
    n[28] = -0.17424871230634e-9;
    n[29] = -0.68762131295531e-18;
    n[30] = 0.14478307828521e-19;
    n[31] = 0.26335781662795e-22;
    n[32] = -0.11947622640071e-22;
    n[33] = 0.18228094581404e-23;
    n[34] = -0.93537087292458e-25;

    L[0] = 0.;
    L[1] = 0.;
    L[2] = 0.;
    L[3] = 0.;
    L[4] = 0.;
    L[5] = 0.;
    L[6] = 0.;
    L[7] = 0.;
    L[8] = 0.;
    L[9] = 1.;
    L[10] = 1.;
    L[11] = 1.;
    L[12] = 1.;
    L[13] = 1.;
    L[14] = 1.;
    L[15] = 2.;
    L[16] = 2.;
    L[17] = 2.;
    L[18] = 2.;
    L[19] = 2.;
    L[20] = 3.;
    L[21] = 3.;
    L[22] = 3.;
    L[23] = 4.;
    L[24] = 4.;
    L[25] = 4.;
    L[26] = 5.;
    L[27] = 8.;
    L[28] = 8.;
    L[29] = 21.;
    L[30] = 23.;
    L[31] = 29.;
    L[32] = 30.;
    L[33] = 31.;
    L[34] = 32.;

    J[0] = -2.;
    J[1] = -2.;
    J[2] = -1.;
    J[3] = 0.;
    J[4] = 1.;
    J[5] = 2.;
    J[6] = 3.;
    J[7] = 4.;
    J[8] = 5.;
    J[9] = -9.;
    J[10] = -7.;
    J[11] = -1.;
    J[12] = 0.;
    J[13] = 1.;
    J[14] = 3.;
    J[15] = -3.;
    J[16] = 0.;
    J[17] = 1.;
    J[18] = 3.;
    J[19] = 17.;
    J[20] = -4.;
    J[21] = 0.;
    J[22] = 6.;
    J[23] = -5.;
    J[24] = -2.;
    J[25] = 10.;
    J[26] = -8.;
    J[27] = -11.;
    J[28] = -6.;
    J[29] = -29.;
    J[30] = -31.;
    J[31] = -38.;
    J[32] = -39.;
    J[33] = -40.;
    J[34] = -41.;

    Pressurevar = pressure_average / Pstar;
    Tau = Tstar / temperature_average;

    GammaPi = 0.;
    for (i = 1; i < 35; i++)
        GammaPi = GammaPi - (n[i]) * (L[i]) *
                  (pow((7.1 - Pressurevar),(L[i] - 1))) * (pow((Tau - 1.222),J[i]));

    rho_0 = pressure_average / (GazConst * temperature_average * Pressurevar * GammaPi);
    return rho_0 + salinity;
}

double LiquidViscosity_HP(double Press,double TempK,double C_0)
{
    double A1,A2,A3,A4,A5,A6,A7,A8;       /*constants*/
    double TempC,TempF, Pbar,Salinity;    /* Temperature [K], Temperature [F], Pressure [bar]*/
    double my_Zero,PsatBar, PsatKPa;      /*my pure water, my saline water, Saturation pressure [bar], Saturation pressure [KPa]*/
                                          /*intermediate values*/
    double sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8, exponent;
    //Link to function from ALR
    Salinity = C_0 / 1000.;
    /***constant values*******/
    A1 = -7.419242;
    A2 = -0.29721;
    A3 = -0.1155286;
    A4 = -0.008685635;
    A5 = 0.001094098;
    A6 = 0.00439993;
    A7 = 0.002520658;
    A8 = 0.0005218684;

    /*Unit conversions*/
    TempC = TempK - 273.15;
    TempF = TempC * 1.8 + 32.0;
    Pbar = Press / 100000.0;
    /*end of units conversion*/

    /*Calculation of the saturation pressure*/
//   sum1=pow((0.65-0.01*TempK),0)*A1;
    sum1 = 1.0 * A1;
//   sum2=pow((0.65-0.01*TempK),1)*A2;
    sum2 = (0.65 - 0.01 * TempK) * A2;
    sum3 = (0.65 - 0.01 * TempK) * (0.65 - 0.01 * TempK) * A3;
    sum4 = std::pow((0.65 - 0.01 * TempK),3) * A4;
    sum5 = std::pow((0.65 - 0.01 * TempK),4) * A5;
    sum6 = std::pow((0.65 - 0.01 * TempK),5) * A6;
    sum7 = std::pow((0.65 - 0.01 * TempK),6) * A7;
    sum8 = std::pow((0.65 - 0.01 * TempK),7) * A8;

    /*intermediate value*/
    exponent = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8;
    exponent = exponent * (374.136 - TempC) / TempK; /*intermediate value*/

    PsatKPa = exp(exponent) * 22088;      /*saturation pressure in kPa*/
    PsatBar = PsatKPa / (1000 * 100000);  /*Saturation pressure in bar*/

    /*Viscosity of pure water in Pa-S*/
    my_Zero = 243.18e-7 *
              (pow(10.,
                   (247.8 / (TempK - 140)))) * (1 + (Pbar - PsatBar) * 1.0467e-6 * (TempK - 305));

    /*Viscosity of saline water in Pa-S*/
    double viscosity = my_Zero *
                (1 - 0.00187 * (sqrt(Salinity)) + 0.000218 *
                 (std::pow(sqrt(Salinity),
                                   5)) +
                 (sqrt(TempF) - 0.0135 *
              TempF) * (0.00276 * Salinity - 0.000344 * (std::pow(sqrt(Salinity),3))));
    return viscosity;
}

/**
 * Dynamic viscosity of high-concentration salt water based on Lever&Jackson(1985),
 * Hassanizadeh(1988), and Mercer&Pinder(1974)
 *
 * Reference: FEFLOW Reference manual pg 31, (1-100)
 *
 * @param c
 * @param T
 * @return
 */
double LiquidViscosity_LJH_MP1(double c,double T, double rho, double my_0)
{
    double omega = c / rho;
    double sigma = (T - 150.0) / 100.0;

    double f = (1 + 0.7063 * sigma - 0.04832 * sigma * sigma * sigma) / (1. + 1.85 * omega - 4.1 * omega * omega + 44.5 * omega * omega * omega);
    double vis = my_0 / f;
    return vis;
}


double LiquidViscosity_LJH_MP2(double c,double T, double rho, double rho0, double my_C0, double my_T0, double my_0)
{
    double omega0 = -1, sigma0, f1_0, f2_0;
    omega0 = my_C0 / rho0;
    sigma0 = (my_T0 - 150.) / 100.;
    f1_0 = 1. + 1.85 * omega0 - 4.1 * omega0 * omega0 + 44.5 * omega0 * omega0 * omega0;
    f2_0 = 1./(1 + 0.7063 * sigma0 - 0.04832 * sigma0 * sigma0 * sigma0);

    double omega, sigma;
    double f1, f2, mu;

    omega = c / rho;
    sigma = (T - 150.) / 100.;

    f1 = f1_0 / (1. + 1.85 * omega - 4.1 * omega * omega + 44.5 * omega * omega * omega);
    f2 = (1 + 0.7063 * sigma - 0.04832 * sigma * sigma * sigma) * f2_0;
    mu = my_0 / (f1*f2);

    return mu;
}

double SalineWaterSpecificHeat_IAPWS_IF97(double Press, double TempK, double Conc)
{
    double Pressurevar, Tau, pressure_average, temperature_average, Tstar, Pstar,GazConst;
    double GammaPi, GammaPiTau, GammaPiPi, GammaTauTau;
    double L[35],J[35],n[35];
    int i;
    double Cp;

    pressure_average = Press;
    temperature_average = TempK;

    Tstar = 1386;

    Pstar = 16.53e6;                      // MPa
    GazConst = 0.461526e3;                //

    n[0] = 0.0;
    n[1] = 0.14632971213167;
    n[2] = -0.84548187169114;
    n[3] = -0.37563603672040e1;
    n[4] = 0.33855169168385e1;
    n[5] = -0.95791963387872;
    n[6] = 0.15772038513228;
    n[7] = -0.16616417199501e-1;
    n[8] = 0.81214629983568e-3;
    n[9] = 0.28319080123804e-3;
    n[10] = -0.60706301565874e-3;
    n[11] = -0.18990068218419e-1;
    n[12] = -0.32529748770505e-1;
    n[13] = -0.21841717175414e-1;
    n[14] = -0.52838357969930e-4;
    n[15] = -0.47184321073267e-3;
    n[16] = -0.30001780793026e-3;
    n[17] = 0.47661393906987e-4;
    n[18] = -0.44141845330846e-5;
    n[19] = -0.72694996297594e-15;
    n[20] = -0.31679644845054e-4;
    n[21] = -0.28270797985312e-5;
    n[22] = -0.85205128120103e-9;
    n[23] = -0.22425281908000e-5;
    n[24] = -0.65171222895601e-6;
    n[25] = -0.14341729937924e-12;
    n[26] = -0.40516996860117e-6;
    n[27] = -0.12734301741641e-8;
    n[28] = -0.17424871230634e-9;
    n[29] = -0.68762131295531e-18;
    n[30] = 0.14478307828521e-19;
    n[31] = 0.26335781662795e-22;
    n[32] = -0.11947622640071e-22;
    n[33] = 0.18228094581404e-23;
    n[34] = -0.93537087292458e-25;

    L[0] = 0.;
    L[1] = 0.;
    L[2] = 0.;
    L[3] = 0.;
    L[4] = 0.;
    L[5] = 0.;
    L[6] = 0.;
    L[7] = 0.;
    L[8] = 0.;
    L[9] = 1.;
    L[10] = 1.;
    L[11] = 1.;
    L[12] = 1.;
    L[13] = 1.;
    L[14] = 1.;
    L[15] = 2.;
    L[16] = 2.;
    L[17] = 2.;
    L[18] = 2.;
    L[19] = 2.;
    L[20] = 3.;
    L[21] = 3.;
    L[22] = 3.;
    L[23] = 4.;
    L[24] = 4.;
    L[25] = 4.;
    L[26] = 5.;
    L[27] = 8.;
    L[28] = 8.;
    L[29] = 21.;
    L[30] = 23.;
    L[31] = 29.;
    L[32] = 30.;
    L[33] = 31.;
    L[34] = 32.;

    J[0] = -2.;
    J[1] = -2.;
    J[2] = -1.;
    J[3] = 0.;
    J[4] = 1.;
    J[5] = 2.;
    J[6] = 3.;
    J[7] = 4.;
    J[8] = 5.;
    J[9] = -9.;
    J[10] = -7.;
    J[11] = -1.;
    J[12] = 0.;
    J[13] = 1.;
    J[14] = 3.;
    J[15] = -3.;
    J[16] = 0.;
    J[17] = 1.;
    J[18] = 3.;
    J[19] = 17.;
    J[20] = -4.;
    J[21] = 0.;
    J[22] = 6.;
    J[23] = -5.;
    J[24] = -2.;
    J[25] = 10.;
    J[26] = -8.;
    J[27] = -11.;
    J[28] = -6.;
    J[29] = -29.;
    J[30] = -31.;
    J[31] = -38.;
    J[32] = -39.;
    J[33] = -40.;
    J[34] = -41.;

    Pressurevar = pressure_average / Pstar;
    Tau = Tstar / temperature_average;

    /*BEGIN:Calculation of GammaPi*/
    GammaPi = 0;

    for (i = 1; i < 35; i++)
        GammaPi = GammaPi - (n[i]) * (L[i]) *
                  (pow((7.1 - Pressurevar),(L[i] - 1.))) * (pow((Tau - 1.222),J[i]));
    /*END: Calculation of GammaPi*/

    /*BEGIN:Calculation of GammaPiTau*/
    GammaPiTau = 0;
    for (i = 1; i < 35; i++)
        GammaPiTau = GammaPiTau - (n[i]) * (L[i]) *
                     (pow((7.1 - Pressurevar),
                          (L[i] - 1.))) * (J[i]) * (pow((Tau - 1.222),(J[i] - 1.)));
    /*END: Calculation of GammaPiTau*/

    /*BEGIN:Calculation of GammaTauTau*/
    GammaTauTau = 0;
    for (i = 1; i < 35; i++)
        GammaTauTau = GammaTauTau + (n[i]) *
                      (pow((7.1 - Pressurevar),
                           (L[i]))) * (J[i]) * (J[i] - 1.) * (pow((Tau - 1.222),(J[i] - 2)));
    /*END: Calculation of GammaTauTau*/

    /*BEGIN:Calculation of GammaPiPi*/
    GammaPiPi = 0;
    for (i = 1; i < 35; i++)
        GammaPiPi = GammaPiPi + (n[i]) * (L[i]) *
                    (L[i] -
                     1) * (pow((7.1 - Pressurevar),(L[i] - 2.))) * (pow((Tau - 1.222),(J[i])));
    /*END: Calculation of GammaPiPi*/

    /*************************Partial derivatives calculation*****************************************/
    /*************************************************************************************************/
    /*************************************************************************************************/

    /*BEGIN: Fluid isobaric heat capacity*/
    Cp = -(Tau * Tau) * (GammaTauTau) * (GazConst);
    /*END: Fluid isobaric heat capacity*/

    /*BEGIN: Fluid isochoric heat capacity*/
    /* Cv is not used currently 9.2003*/
    //WW Cv = (- (pow(Tau,2))* (GammaTauTau) + pow(GammaPi - Tau * (GammaPiTau),2) / GammaPiPi) * GazConst;
    /*BEGIN: Fluid isochoric heat capacity*/

    return Cp;
}

double SalineWaterThermalConductivity_IAPWS_IF97(double Press, double TempK, double Conc)
{
    int i, j;
    double TauTC, Nabla, Delta, Nabla0, Nabla1, Nabla2;
    double heat_conductivity, Rho, temperature_average, pressure_average, viscosity;
    double Rhostar, TstarTC, Lambdastar; //WW, Pstar1;
    double nZero[4];
    double n[5][6];
    double A1,A2,A3,A4,A5,A6,A7,A8;       /*constants*/
    double TempC, TempF, Pbar, Salinity;  /* Temperature [K], Temperature [F], Pressure [bar]*/
    double my_Zero,PsatBar, PsatKPa;      /*my pure water, my saline water, Saturation pressure [bar], Saturation pressure [KPa]*/
                                          /*intermediate values*/
    double sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8, exponent;

    /*************************************************************************************************/
    /*************************************************************************************************/
    /*************************Partial derivatives calculation*****************************************/
    double GammaPi, GammaPiTau, GammaPiPi, Pii, Tau,  GazConst;
    double LGamma[35];
    double JGamma[35];
    double nGamma[35];
    double Tstar, Pstar;
    double TstarTilda, PstarTilda, RhostarTilda;
    double First_derivative, Second_derivative;

    pressure_average = Press;
    temperature_average = TempK;
    Salinity = Conc;
    Tstar = 1386;
    Pstar = 16.53e6;                      // MPa
    GazConst = 0.461526e3;                //!!!!Given by equation (2.1)
    TstarTilda = 647.226;
    PstarTilda = 22.115e6;
    RhostarTilda = 317.763;
    Lambdastar = 0.4945;
    Rhostar = 317.763;
    TstarTC = 647.226;
    //WW Pstar1 = 22.115e-6;

    /*BEGIN: reduced dimensions*/
    TauTC = TstarTC / temperature_average;
    //Delta = Rho / Rhostar;
    //WW PiiTC = pressure_average / Pstar1;
    /*END: reduced dimensions*/

    nGamma[1] = 0.14632971213167;
    nGamma[2] = -0.84548187169114;
    nGamma[3] = -0.37563603672040e1;
    nGamma[4] = 0.33855169168385e1;
    nGamma[5] = -0.95791963387872;
    nGamma[6] = 0.15772038513228;
    nGamma[7] = -0.16616417199501e-1;
    nGamma[8] = 0.81214629983568e-3;
    nGamma[9] = 0.28319080123804e-3;
    nGamma[10] = -0.60706301565874e-3;
    nGamma[11] = -0.18990068218419e-1;
    nGamma[12] = -0.32529748770505e-1;
    nGamma[13] = -0.21841717175414e-1;
    nGamma[14] = -0.52838357969930e-4;
    nGamma[15] = -0.47184321073267e-3;
    nGamma[16] = -0.30001780793026e-3;
    nGamma[17] = 0.47661393906987e-4;
    nGamma[18] = -0.44141845330846e-5;
    nGamma[19] = -0.72694996297594e-15;
    nGamma[20] = -0.31679644845054e-4;
    nGamma[21] = -0.28270797985312e-5;
    nGamma[22] = -0.85205128120103e-9;
    nGamma[23] = -0.22425281908000e-5;
    nGamma[24] = -0.65171222895601e-6;
    nGamma[25] = -0.14341729937924e-12;
    nGamma[26] = -0.40516996860117e-6;
    nGamma[27] = -0.12734301741641e-8;
    nGamma[28] = -0.17424871230634e-9;
    nGamma[29] = -0.68762131295531e-18;
    nGamma[30] = 0.14478307828521e-19;
    nGamma[31] = 0.26335781662795e-22;
    nGamma[32] = -0.11947622640071e-22;
    nGamma[33] = 0.18228094581404e-23;
    nGamma[34] = -0.93537087292458e-25;

    LGamma[1] = 0.;
    LGamma[2] = 0.;
    LGamma[3] = 0.;
    LGamma[4] = 0.;
    LGamma[5] = 0.;
    LGamma[6] = 0.;
    LGamma[7] = 0.;
    LGamma[8] = 0.;
    LGamma[9] = 1.;
    LGamma[10] = 1.;
    LGamma[11] = 1.;
    LGamma[12] = 1.;
    LGamma[13] = 1.;
    LGamma[14] = 1.;
    LGamma[15] = 2.;
    LGamma[16] = 2.;
    LGamma[17] = 2.;
    LGamma[18] = 2.;
    LGamma[19] = 2.;
    LGamma[20] = 3.;
    LGamma[21] = 3.;
    LGamma[22] = 3.;
    LGamma[23] = 4.;
    LGamma[24] = 4.;
    LGamma[25] = 4.;
    LGamma[26] = 5.;
    LGamma[27] = 8.;
    LGamma[28] = 8.;
    LGamma[29] = 21.;
    LGamma[30] = 23.;
    LGamma[31] = 29.;
    LGamma[32] = 30.;
    LGamma[33] = 31.;
    LGamma[34] = 32.;

    JGamma[1] = -2.;
    JGamma[2] = -1.;
    JGamma[3] = 0.;
    JGamma[4] = 1.;
    JGamma[5] = 2.;
    JGamma[6] = 3.;
    JGamma[7] = 4.;
    JGamma[8] = 5.;
    JGamma[9] = -9.;
    JGamma[10] = -7.;
    JGamma[11] = -1.;
    JGamma[12] = 0.;
    JGamma[13] = 1.;
    JGamma[14] = 3.;
    JGamma[15] = -3.;
    JGamma[16] = 0.;
    JGamma[17] = 1.;
    JGamma[18] = 3.;
    JGamma[19] = 17.;
    JGamma[20] = -4.;
    JGamma[21] = 0.;
    JGamma[22] = 6.;
    JGamma[23] = -5.;
    JGamma[24] = -2.;
    JGamma[25] = 10.;
    JGamma[26] = -8.;
    JGamma[27] = -11.;
    JGamma[28] = -6.;
    JGamma[29] = -29.;
    JGamma[30] = -31.;
    JGamma[31] = -38.;
    JGamma[32] = -39.;
    JGamma[33] = -40.;
    JGamma[34] = -41.;

    Pii = pressure_average / Pstar;
    Tau = Tstar / temperature_average;

    /*BEGIN:Calculation of GammaPi*/
    GammaPi = 0;

    for (i = 1; i < 35; i++)
        GammaPi = GammaPi - (nGamma[i]) * (LGamma[i]) *
                  (pow((7.1 - Pii),(LGamma[i] - 1))) * (pow((Tau - 1.222),JGamma[i]));
    /*END: Calculation of GammaPi*/

    /*BEGIN:Calculation of GammaPiTau*/
    GammaPiTau = 0;
    for (i = 1; i < 35; i++)
        GammaPiTau = GammaPiTau - (nGamma[i]) * (LGamma[i]) *
                     (pow((7.1 - Pii),
                          (LGamma[i] -
                           1))) * (JGamma[i]) * (pow((Tau - 1.222),(JGamma[i] - 1)));
    /*END: Calculation of GammaPiTau*/

    /*BEGIN:Calculation of GammaPiPi*/
    GammaPiPi = 0;
    for (i = 1; i <= 34; i++)
        GammaPiPi = GammaPiPi + (nGamma[i]) * (LGamma[i]) *
                    (LGamma[i] -
                     1) *
                    (pow((7.1 - Pii),(LGamma[i] - 2))) * (pow((Tau - 1.222),(JGamma[i])));
    /*END: Calculation of GammaPiPi*/

    /*BEGIN:Calculation of derivative*/
    First_derivative =
            ((TstarTilda) * (Pstar) *
             ((GammaPiTau) * (Tstar) - (GammaPi) *
          (temperature_average))) /
            (PstarTilda * temperature_average * temperature_average * GammaPiPi),
    Second_derivative =
            ((-1) * (PstarTilda) *
             (GammaPiPi) ) /
            ( (RhostarTilda) * (temperature_average) * (GazConst) * ((GammaPi * GammaPi)));
    /*End:Calculation of derivative*/

    /*BEGIN: Calculation of density*/
    Rho = pressure_average / (GazConst * (temperature_average) * (Pii) * (GammaPi));
    /*END: Calculation of density*/

    /*************************Partial derivatives calculation*****************************************/
    /*************************************************************************************************/
    /*************************************************************************************************/

    /*BEGIN: Constant values*/

    Lambdastar = 0.4945;
    Rhostar = 317.763;
    TstarTC = 647.226;
    //WW Pstar1 = 22.115e6;

    nZero[0] = 0.1e1;
    nZero[1] = 0.6978267e1;
    nZero[2] = 0.2599096e1;
    nZero[3] = -0.998254;

    n[0][0] = 0.13293046e1;
    n[0][1] = -0.40452437;
    n[0][2] = 0.24409490;
    n[0][3] = 0.18660751e-1;
    n[0][4] = -0.12961068;
    n[0][5] = 0.44809953e-1;

    n[1][0] = 0.17018363e1;
    n[1][1] = -0.22156845e1;
    n[1][2] = 0.16511057e1;
    n[1][3] = -0.76736002;
    n[1][4] = 0.37283344;
    n[1][5] = -0.11203160;

    n[2][0] = 0.52246158e1;
    n[2][1] = -0.10124111e2;
    n[2][2] = 0.49874687e1;
    n[2][3] = -0.27297694;
    n[2][4] = -0.43083393;
    n[2][5] = 0.13333849;

    n[3][0] = 0.87127675e1;
    n[3][1] = -0.95000611e1;
    n[3][2] = 0.43786606e1;
    n[3][3] = -0.91783782;
    n[3][4] = 0;
    n[3][5] = 0;

    n[4][0] = -0.18525999e1;
    n[4][1] = 0.93404690;
    n[4][2] = 0;
    n[4][3] = 0;
    n[4][4] = 0;
    n[4][5] = 0;
    /*END: Constant values*/

    /*BEGIN: reduced dimensions*/
    TauTC = TstarTC / temperature_average;
    Delta = Rho / Rhostar;
    //WW PiiTC = pressure_average / Pstar1;
    /*END: reduced dimensions*/

    /*BEGIN: Nabla0*/
    Nabla0 = 0;
    for (i = 0; i <= 3; i++)
        Nabla0 = Nabla0 + (nZero[i]) * (std::pow(TauTC,i));

    Nabla0 = Nabla0 * (sqrt(TauTC));
    Nabla0 = 1 / Nabla0;
    /*END: Nabla0*/

    /*BEGIN: Nabla1*/
    Nabla1 = 0;
    for (i = 0; i <= 4; i++)
    {
        const double t(std::pow(TauTC - 1, i)); // TF
        for (j = 0; j <= 5; j++)
            Nabla1 = Nabla1 + (n[i][j]) * (t * std::pow(Delta - 1, j));
    }

    Nabla1 = Delta * (Nabla1);
    Nabla1 = exp(Nabla1);
    /*END: Nabla1*/

    /*Calculate Viscosity*/
    /*Link to function from ALR*/

    TempK = temperature_average;
    Press = pressure_average;
    Salinity = Conc / 1000.;

    /***constant values*******/
    A1 = -7.419242;
    A2 = -0.29721;
    A3 = -0.1155286;
    A4 = -0.008685635;
    A5 = 0.001094098;
    A6 = 0.00439993;
    A7 = 0.002520658;
    A8 = 0.0005218684;

    /*Unit conversions*/
    TempC = TempK - 273.15;
    TempF = TempC * 1.8 + 32.0;
    if (TempF < 0.0)
        TempF = 0.0;
    Pbar = Press / 100000.0;
    /*end of units conversion*/

    /*Calculation of the saturation pressure*/
// TF   sum1=pow((0.65-0.01*TempK),0)*A1;
    sum1 = A1;
// TF  sum2=pow((0.65-0.01*TempK),1)*A2;
    sum2 = (0.65 - 0.01 * TempK) * A2;
    sum3 = (0.65 - 0.01 * TempK) * (0.65 - 0.01 * TempK) * A3;
    sum4 = std::pow((0.65 - 0.01 * TempK),3) * A4;
    sum5 = std::pow((0.65 - 0.01 * TempK),4) * A5;
    sum6 = std::pow((0.65 - 0.01 * TempK),5) * A6;
    sum7 = std::pow((0.65 - 0.01 * TempK),6) * A7;
    sum8 = std::pow((0.65 - 0.01 * TempK),7) * A8;

    /*intermediate value*/
    exponent = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8;
    exponent = exponent * (374.136 - TempC) / TempK; /*intermediate value*/

    PsatKPa = exp(exponent) * 22088;      /*saturation pressure in kPa*/
    PsatBar = PsatKPa / (1000 * 100000);  /*Saturation pressure in bar*/

    /*Viscosity of pure water in Pa-S*/
    my_Zero = 243.18e-7 *
              (pow(10.,
                   (247.8 / (TempK - 140)))) * (1 + (Pbar - PsatBar) * 1.0467e-6 * (TempK - 305));

    /*Viscosity of saline water in Pa-S*/
    viscosity = my_Zero *
                (1 - 0.00187 * (sqrt(Salinity)) + 0.000218 *
                 (std::pow(sqrt(Salinity),
                                   5)) +
                 (sqrt(TempF) - 0.0135 *
              TempF) * (0.00276 * Salinity - 0.000344 * (std::pow(sqrt(Salinity),3))));

    /* End of viscosity function*/

    /*BEGIN: Nabla2*/
    Nabla2 = 0.0013848 /
             ((viscosity) /
              55.071e-6) *
             (1.0 /
              (TauTC *
           Delta) *
              (TauTC *
           Delta)) *
             (First_derivative *
              First_derivative) *
             (pow((Delta * (Second_derivative)),0.4678)) * (sqrt(Delta)) * exp(
            -18.66 * ((1 / TauTC - 1) * (1 / TauTC - 1)) - (std::pow(Delta - 1,4)));
    /*END: Nabla2*/

    /*BEGIN: Nabla => heat_conductivity*/
    Nabla = Nabla0 * (Nabla1) + Nabla2;
    heat_conductivity = Nabla * (Lambdastar);
    /*END: Nabla => Lambda*/

    return heat_conductivity;
}

} //end
