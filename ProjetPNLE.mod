/*********************************************
 * OPL 22.1.1.0 Model
 * Author: ineshafassamaiza
 * Creation Date: Jan 15, 2025 at 4:38:13 PM
 *********************************************/
// Sets
int nbClients = ...;
int nbVehicules = ...;

range K1 = 1..3;             // Long-term vehicles
range K2 = 4..nbVehicules;   // Short-term vehicles
range C = 1..nbClients;      // Clients
range N = 0..nbClients;      // Nodes (0 = depot)

// Parameters
float d[N][N]=...;      // Distance matrix
float v[N]=...;         // Demand of each client
float cap[K1]=...;      // Capacity of long-term vehicles
float speed[K1]=...;    // Speed of long-term vehicles
float fc[1..5]=...;     // Fixed costs for all vehicles (K1 + K2)
float Tsoft[K1]=...;    // Soft time limit for long-term vehicles
float Thard[K1]=...;    // Hard time limit for long-term vehicles
float Dsoft[K1]=...;    // Soft distance limit for long-term vehicles
float Dhard[K2]=...;    // Hard distance limit for short-term vehicles
float Pt[K1]=...;       // Time penalty for long-term vehicles
float Pd[K1]=...;       // Distance penalty for long-term vehicles

// Decision Variables
dvar boolean x[N][N][K1];           // Long-term vehicle arcs
dvar boolean z[C][K2];              // Short-term vehicle assignment
dvar float+ et[K1];                 // Excess time for long-term vehicles
dvar float+ ed[K1];                 // Excess distance for long-term vehicles         
dvar float+ u[N][K1];				// Updated capacity

// Objective Function
minimize 
  sum(k in K1, j in N) (fc[k]*x[0][j][k] + Pt[k] * et[k] + Pd[k] * ed[k]) + 
  sum(k in K2, j in C) (z[j][k] * fc[k]);

// Constraints
subject to {
  // Each client is visited exactly once
  forall(j in N: j != 0)
    sum(i in N: i != j, k in K1) x[i][j][k] + sum(k in K2) z[j][k] == 1;

  // Flow conservation for long-term vehicles
  forall(k in K1, j in N)
    sum(i in N: i != j) x[i][j][k] == sum(i in N) x[j][i][k];

  // Capacity constraints for long-term vehicles
  forall(k in K1)
    sum(i in N, j in N: j != 0 && i != j) v[i] * x[i][j][k] <= cap[k];

  // Time constraints for long-term vehicles
  forall(k in K1)
    et[k] <=  Thard[k]-Tsoft[k];

  forall(k in K1) 
  	sum(i in N, j in N: i != j) (d[i][j]/speed[k]) * x[i][j][k] <= et[k] + Tsoft[k] ;
  	
 // Distance constraints for long-term vehicles
  forall(k in K1) 
  	sum(i in N, j in N: i != j) d[i][j] * x[i][j][k] <= ed[k] + Dsoft[k] ;
  	
// Distance constraints for long-term vehicles
  forall(k in K1, i in N, j in N: i != j) {
    u[i][k]-u[j][k] + cap[k]*x[i][j][k] <= cap[k]-v[i];
    v[i] <= u[i][k] <= cap[k];
  }
}


// Résolution du modèle
execute AfficherSolution{
  writeln("Solution trouvée !");
  writeln("Valeur de la fonction objective : ", cplex.getObjValue());
}