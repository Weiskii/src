# MILP HEMS MODEL (ELEC5206 SUSTAINABLE ENERGY SYSTEMS)
# Donald Azuatalam (donald.azuatalam@sydney.edu.au)
# Dr. Gregor Verbic and Dr. Archie Chapman

# SETS
set D; /* Set of half-hourly time steps for 2 days */

# PARAMETERS
# Tariff data
param c_g{d in D} >= 0; /* Time of use tariff */
param c_flat = 0.31317; /* Flat tariff */
param c_pv = 0.08; /* Feed-in-tariff */

# Demand Data
param Pd{d in D} >= 0;  /* Electrical demand (daily load profile) in kW */
param Ppv{d in D} >= 0;  /* Solar PV output (daily PV output profile) in kW */
param PgM = 15;  /* Maximum capacity of grid connection in kW */

# Battery and Inverter Data
param ebM >= 0;  /* Battery maximum storage limit [kWh] */
param ebm >= 0;  /* Battery minimum storage limit [kWh]*/
param eb1 >= 0;  /* Start-of-day battery state of charge (SOC) */
param ebN = 0.2*ebM;  /* End-of-day (day 2 in a 2-day rolling horizon) battery state of charge (SOC) -  20% Max SOC*/
param PbM >= 0;  /* Battery maximum charging rate [kW] */  
param Pbm = PbM; /* Battery maximum discharge rate [kW] */
param etaBc >= 0; /* Battery charging efficiency */
param etaBd = etaBc; /* Battery discharging efficiency */
param etaI = 1; /* Inverter efficiency. This is because the inverter efficiency is already accounted for in the PV data */ 
param dt = 24/48; /* Half hourly time steps */
param N = 96; /* Total number of time-slots for 2 days */

# VARIABLES
var Pgplus{d in D} >= 0, <= PgM;  /* Power flowing from grid to customer in kW */
var Pgminus{d in D} >= 0, <= PgM;  /* Power flowing from customer to grid in kW */
var Pbplus{d in D} >= 0, <= PbM;  /* Battery charge power in kW */
var Pbminus{d in D} >= 0, <= Pbm;  /* Battery discharge power in kW */
var eb{d in D} >= ebm, <= ebM;  /* Battery state of charge */
var dg{d in D} binary;   /* Direction of grid power flow (0: demand->grid, 1: grid->demand) */
var sb{d in D} binary;  /* Battery charging status (0: discharge, 1: charge) */

# OBJECTIVE FUNCTION
minimize cost: sum{d in D} (dt*c_g[d]*Pgplus[d] - dt*c_pv*Pgminus[d]); /* Objective fucntion with Time of Use Tariff */
#minimize cost: sum{d in D} (dt*c_flat*Pgplus[d] - dt*c_pv*Pgminus[d]); /* Objective function with Flat Tariff */
  
subject to
# CONSTRAINTS
# Equality Constraints
Equality_constraint1 {d in D}: Pgplus[d] - Pgminus[d] = etaI*(etaBc*Pbplus[d] - (1/etaBd)*Pbminus[d]) - etaI*Ppv[d] + Pd[d]; /* Power balance constraint */
Equality_constraint2 {d in D}: eb[d] = if d = 1 then eb1 else if d = N then ebN else eb[d-1] + dt*etaBc*Pbplus[d-1] - dt*(1/etaBd)*Pbminus[d-1]; /* Battery operation */

# Inequality Constraints
Inequality_constraint1 {d in D}: Pgplus[d] <= PgM*dg[d];  /* Power limiting strategy  */
Inequality_constraint2 {d in D}: Pgminus[d] <= PgM*(1 - dg[d]);  /* Import/Export constraint - can't import and export at the same time */
Inequality_constraint3 {d in D}: Pbplus[d] <= PbM*sb[d];   /* Battery charging constraint */  
Inequality_constraint4 {d in D}: Pbminus[d] <= Pbm*(1 - sb[d]); /* Battery charge/discharge constraint - can't charge and discharge at the same time */