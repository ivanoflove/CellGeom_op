/*
 * * Copyright 1987-2023 ANSYS, Inc. All Rights Reserved.
 */

#include "sofc.h"

/*
 * Change these to 1 to use user-defined properties
 * Relevant functions need to be modified.
 */
#define USER_DEFINED_ZONE_CONDUCTIVITIES 1
#define USER_DEFINED_CONTACT_RESIST 0
#define USER_DEFINED_LEAKAGE_CURRENT 0

/****************************************************************/
/****************************************************************/
/*  CONSTITUTIVE LAWS - NERNST, ACTIVATION, CONDUCTIVITY, ETC.  */
/****************************************************************/
/****************************************************************/

real Solve_Butler_Volmer_NR (real i, real i0, real A, real B, real *detadi);

real Nernst (real sp_a[],  real sp_c[],
             real P_a, real P_c, real T_a, real T_c, real *DNDY)
{

  real E, G, G0, T, Q;
  real X_a[NO_OF_FC_SPECIES], X_c[NO_OF_FC_SPECIES];
  int i, mm;

  Domain *domain = Get_Domain(1);
  Material *mat = mixture_material(domain);

  for (mm = 0; mm < MIXTURE_NSPECIES(mat); mm++)
  {
    sp_a[mm] = MAX(sp_a[mm], 1.0e-30);
    sp_c[mm] = MAX(sp_c[mm], 1.0e-30);
  }

  Mole_Fraction(mat, sp_a, X_a);
  Mole_Fraction(mat, sp_c, X_c);

  for (mm = 0; mm < MIXTURE_NSPECIES(mat); mm++)
  {
    X_a[mm] = MAX(X_a[mm], 1.0e-30);
    X_c[mm] = MAX(X_c[mm], 1.0e-30);
  }

  T = 0.5 * (T_a + T_c);

  P_a = fabs(P_a / 101325.0);
  P_c = fabs(P_c / 101325.0);

  if ((Eq_CO == -1) || (Eq_CO2 == -1) || CO_Disabled) /* H2-H20 version */
  {
    G0 = (-1.192684E5 - 15.950548 * T - 9.909308E-4 * T * T
          - 9.198735E-8 * T * T * T + 5.786936 * log(T) * T)
         * 1.0e3 * 2.0;   /* J/K-kmole H2 */

    if (Electrolysis_Mode)
     Q = (P_a * X_a[Eq_H2]) * sqrt(P_c * X_c[Eq_O2]) / (P_a * X_a[Eq_H2O]);
    else 
     Q = (P_a * X_a[Eq_H2]) * sqrt(P_c * X_c[Eq_O2]) / (P_a * X_a[Eq_H2O]);
    

    G = G0 - Rgas * T * log(Q); /* J/K-kmole H2 */

    E = G / (-2.0 * Faraday);
  }
  else /* 3 H2 + 2 O2 + 1 CO --> 3 H2O + CO2 */
  {
    G0 = 4.0 * (-1.244683E5 - 11.554384 * T - 0.002622 * T * T
                + 1.828235E-7 * T * T * T + 6.056462 * log(T) * T)
         * 1.0e3 * 2.0;   /* J/K-kmole H2 */

    if (Electrolysis_Mode)
     Q = 1.0;
    else 
     Q = (X_a[Eq_H2] * X_a[Eq_H2] * X_a[Eq_H2]) * (X_a[Eq_CO]) *
        X_c[Eq_O2] * X_c[Eq_O2] * P_c * P_c
        /
        (X_a[Eq_H2O] * X_a[Eq_H2O] * X_a[Eq_H2O] *
         X_a[Eq_CO2]);

    G = G0 - Rgas * T * log(Q); /* J/K-kmole H2 */

    E = G / (-8.0 * Faraday);
  }
  if (E < 0.0) E = 0.0;

  for (i = 0; i < MIXTURE_NSPECIES(mat); i++)
  {
    DNDY[i] = 0.0;
  }
   if(!Electrolysis_Mode) 
  {
    real DE_DG, DG_DQ, DQ_DXH2, DQ_DXO2, DQ_DXH2O;
    real DXH2_DYH2, DXO2_DYO2, DXH2O_DYH2O;
    real DM_DY, Moles_a, Moles_c;

    DE_DG = -1.0 / (2.0 * Faraday);
    DG_DQ = -1.0 * Rgas * T / Q;

    DQ_DXH2 = sqrt(P_c * X_c[Eq_O2]) / X_a[Eq_H2O];
    DQ_DXO2 = 0.5 * X_a[Eq_H2] * P_c / (X_a[Eq_H2O] * sqrt(P_c * X_c[Eq_O2]));
    DQ_DXH2O = 0.0;


    Moles_a = 0.0;
    for (mm = 0; mm < MIXTURE_NSPECIES(mat); mm++)
    {
      Moles_a += sp_a[mm] / MW[mm];
    }

    Moles_c = 0.0;
    for (mm = 0; mm < MIXTURE_NSPECIES(mat); mm++)
    {
      Moles_c += sp_c[mm] / MW[mm];
    }

    DXH2_DYH2 = (Moles_a * MW[Eq_H2] - sp_a[Eq_H2]) /
                (Moles_a * Moles_a * MW[Eq_H2] * MW[Eq_H2]);

    DXO2_DYO2 = (Moles_c * MW[Eq_O2] - sp_c[Eq_O2]) /
                (Moles_c * Moles_c * MW[Eq_O2] * MW[Eq_O2]);

    DXH2O_DYH2O = (Moles_a * MW[Eq_H2O] - sp_a[Eq_H2O]) /
                  (Moles_a * Moles_a * MW[Eq_H2O] * MW[Eq_H2O]);

    DNDY[Eq_H2O] = DE_DG * DG_DQ * DQ_DXH2O * DXH2O_DYH2O;
    DNDY[Eq_H2] = DE_DG * DG_DQ * DQ_DXH2 * DXH2_DYH2;
    DNDY[Eq_O2] = DE_DG * DG_DQ * DQ_DXO2 * DXO2_DYO2;
  }

  if ((E < 10.0) && (E > -5.0))
    return E;
  else
  {
    Message0("Nernst out of bounds-Value: E:%e, Q: %e, P_a: %e, P_c: %e\n",
             E, Q, P_a, P_c);
    return 1.0e-10;
  }
}

/*  Activation Overpotential  */

real Activation (real Y[], real P, real T, real i, real i_0,
                 real alpha_a, real alpha_b, int species, real *dA_di,
                 cell_t c0, Thread *t0)
{
  real A, alpha = 0.5, n = 2.0;
  real X[NO_OF_FC_SPECIES], i_0_p = 1.0, sign;
  Domain *domain = Get_Domain(1);
  Material *mat = mixture_material(domain);
  int mm;

  *dA_di = 0.0;

  if (Temp_Dependent_I_0 && (species == Eq_O2))  /*for cathode, original*/
    i_0 = A_I_0 / exp ( B_I_0 / (Rgas * T));
  if (an_Temp_Dependent_I_0 && (species == Eq_H2))   /*for anode, added in 23.1*/
    i_0 = an_A_I_0 / exp ( an_B_I_0 / (Rgas * T));

  sign = 1.0;
  if (i < 0.0) sign = -1.0;

  Mole_Fraction(mat, Y, X);

  for (mm = 0; mm < MIXTURE_NSPECIES(mat); mm++)
    X[mm] = MAX(X[mm], 1.0e-30);


  if (species == Eq_H2)
  {
    n = 2.0;
    if (Electrolysis_Mode)
      i_0_p = pow(X[Eq_H2] / MoleF_ref_H2, exponent_h2) *
            pow(X[Eq_H2O] / MoleF_ref_H2O, exponent_h2o);
    else 
      i_0_p = pow(X[Eq_H2] / MoleF_ref_H2, exponent_h2) *
            pow(X[Eq_H2O] / MoleF_ref_H2O, exponent_h2o);
  }
  if (species == Eq_O2)
  {
    n = 4.0;
    if (Electrolysis_Mode)
      i_0_p = pow(X[Eq_O2] / MoleF_ref_O2, exponent_o2); 
    else 
      i_0_p = pow(X[Eq_O2] / MoleF_ref_O2, exponent_o2);
  }

  i_0_p *= fabs(i_0);

  if (i_0_p < 0.0) Internal_Error("BAD VALUE IN ACTIVATION!");

  if (alpha_a == alpha_b)
  {


    if (i_0_p < 0.0) Internal_Error("BAD VALUE IN ACTIVATION!");


    {
      real asinh_value, arg, log_arg;
      arg = (0.5 * fabs(i) / i_0_p);
      asinh_value = log(fabs(arg + sqrt(1.0 + arg * arg)));
      A = Rgas * T * asinh_value / (n * Faraday * alpha_a);
      *dA_di = Rgas * T /
               (Faraday * n * sqrt(1.0 + arg * arg) * 2.0 * i_0_p * alpha_a);
    }





    if (fabs(i) < 1e-15) *dA_di = 0.0;

    if ((A < 3.0) && (A > -3.0))
    {
      return A * sign;
    }
    else
      return 0.0;
  }
  else /* Use Newton Method */
  {
    real inv_P = n * Faraday / (Rgas * T);

    A = Solve_Butler_Volmer_NR(i, i_0_p, alpha_a * inv_P, alpha_b * inv_P, dA_di);
  }

  if ((A < 3.0) && (A > -3.0))
    return A;
  else
    return 0.0;

}


/*  Resistivity of Electrolyte  */

real Resistivity (real T)
{
  real resist;

  if (!mem_con_on)
  {
    resist = membrane_resistivity;
  }
  else
  {
    /*real temp = T;
    if (temp < 1073) temp = 1073;
    if (temp > 1373) temp = 1373;
    resist = ((0.3685 + 0.002838 * exp(10300 / temp)) / 100.0);*/
    real temp = T;
    resist = (1 / (3.6e7 / temp * exp(-8e4 / (Rgas * temp))));
    /*resist = (1 / (3.34e4 * exp(-10300 / temp)));*/
  }

  return resist;
}

/* Leakage Current */
/*
 * ft_an, f_an, ft_ca, f_ca are anode and cathode threads and faces, respectively
 */
real Leakage_Current_Density (face_t f_an, Thread *ft_an, face_t f_ca, Thread *ft_ca)
{
  Thread *t_an = ft_an->t0, *t_ca = ft_ca->t0;
  cell_t c_an = F_C0(f_an, ft_an), c_ca = F_C0(f_ca, ft_ca);
  real T_an, T_ca, T_ave;
  real leak_i = 0.0;

  if (!USER_DEFINED_LEAKAGE_CURRENT)
  {
    /* The default value is Leakage_Current_Density_const, taken from the user interface.*/
    leak_i = Leakage_Current_Density_const;
  }
  else
  {
    /* The default value may be overwritten as a function of Temperature or other parameters. */
    /* user needs to add code here... below is an example*/
    T_an = F_T(f_an, ft_an);
    T_ca = F_T(f_ca, ft_ca);
    T_ave = 0.5 * (T_an + T_ca);
    leak_i = Leakage_Current_Density_const * pow(T_ave, 0.1);
  }

  return leak_i;
}


real Solve_Butler_Volmer_NR (real i, real i0, real A, real B, real *detadi)
{
  real eta = 0.0;
  real res = 1.0;
  real dfdeta0, f1, f0;
  real termA, termB;
  int count = 0;

  termA = i0 * exp(A * eta);
  termB = i0 * exp(-B * eta);
  f0 = termA - termB - i;
  dfdeta0 = A * termA + B * termB;


  do
  {
    eta -= f0 / dfdeta0;

    termA = i0 * exp(A * eta);
    termB = i0 * exp(-B * eta);
    f0 = termA - termB - i;
    dfdeta0 = A * termA + B * termB;

    res = fabs(f0);
    count++;
  }
  while ((res > 1.e-10) && (count < 200));

  *detadi = 1.0 / dfdeta0;
  return eta;
}

real CONDUCTIVITY_CELL(cell_t c, Thread *t)
{
  real c_cond = 0.0;

  if (!USER_DEFINED_ZONE_CONDUCTIVITIES)   /* Returns contant values as defined in Efield Panel */
    c_cond = CONDUCTIVITY_CELL_THREAD(t); 
  else                                     /* User function to compute condutivity at this cell */
   {  
      int zone_id = THREAD_ID(t);
      switch(zone_id)
      {
        case 1: /*anode*/
        {
          real temp = C_T(c, t);
          /*c_cond = 9.7e4 / temp * exp(-2100 / temp);*/
          c_cond = 9.5e7 / temp * exp(-1150 / temp);
          break;
        }
        case 2: /*cathode*/
        {
          real temp = C_T(c, t);
          c_cond = 4.2e7 / temp * exp(-1200 / temp);
          /*c_cond = 9.5e7 / temp * exp(-1150 / temp);*/
          break;
        }
        default:
        {
          c_cond = 7.768e5;
        }
      }


   }

  return c_cond;
}

real CONTACT_RESIST_FACE (face_t f, Thread *t)
{
  real contact_r = 0.0;

  if (!USER_DEFINED_CONTACT_RESIST)  /* Default returns contant values as defined in Efiled Panel (GUI)  */
    contact_r = CONTACT_RESIST_FACE_THREAD(t);
  else                               /* User function to compute contance resistance at this face  */
  {
     real temp = F_T(f, t); /* EXAMPLE - DON'T USE */
     contact_r = 10.0 * temp;
  }

  return contact_r;
}

real h2_co_split_func(cell_t c_an, Thread *t_an)
{

  /* By default, h2_co_split at the next-to-"anode electrolyte" cell
   * is computed as X_h2/(X_h2+X_co);
   * Users can overwrite this value in this udf */

  Domain *domain = Get_Domain(1);
  Material *mat = mixture_material(domain);

  real h2_co_split = 1.0;
  real Y_an[NO_OF_FC_SPECIES], mole_frac[NO_OF_FC_SPECIES];
  int mm;

  for (mm = 0; mm < MIXTURE_NSPECIES(mat); mm++)
  {
    Y_an[mm] = MAX(1.0e-30, C_YI(c_an, t_an, mm));
  }
  Mole_Fraction(mat, Y_an, mole_frac);
  for (mm = 0; mm < MIXTURE_NSPECIES(mat); mm++)
  {
    mole_frac[mm] = MAX(1.0e-30, mole_frac[mm]);
  }

  h2_co_split = mole_frac[Eq_H2] / (mole_frac[Eq_H2] + mole_frac[Eq_CO]);

  /* users can specify their own "h2_co_split" below */

  return h2_co_split;
}


/*=====================================================================*/
/*
 * Interface to FLUENT. Do not modify.
 *
 */

DEFINE_DIFFUSIVITY(diffusivity, c, t, i)
{
  return sofc_diffusivity(c, t, i);
}

DEFINE_DIFFUSIVITY(E_Conductivity, c, t, i)
{
  return CONDUCTIVITY_CELL(c, t);
}

DEFINE_SOURCE(source, c, t, dS, eqn)
{
  return sofc_source(c, t, dS, eqn);
}

DEFINE_ADJUST(adjust_function, domain)
{
  sofc_adjust_function(domain);
  return;
}

DEFINE_INIT(sofc_init, domain)
;
DEFINE_EXECUTE_ON_LOADING(sofc_autorun, sofclib)
;
DEFINE_EXECUTE_AFTER_DATA(sofc_counter_reset, sofclib)
;
DEFINE_ON_DEMAND(list_sofc_udm)
{
  int n;

  Message0("\n\nSOFC UDMs are defined as follows:\n");
  for (n = 0; n < NO_OF_UDMS; n++)
    Message0("UDM %2d: %32s [%s]\n", n, user_memory_vars[n].name, user_memory_vars[n].units);
  Message0("\n");
}
