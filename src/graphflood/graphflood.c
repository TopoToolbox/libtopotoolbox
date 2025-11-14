#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "gf_utils.h"
#include "graphflood/define_types.h"
#include "pq_maxheap.h"
#include "topotoolbox.h"

/*
 * GRAPHFLOOD FULL IMPLEMENTATION
 *
 * This module implements the complete GraphFlood algorithm for transient flood
 * simulation. Unlike the metrics function which calculates steady-state flow,
 * this implementation solves the time-dependent shallow water equations using
 * an iterative approach.
 *
 * The algorithm supports both:
 * - Single Flow Direction (SFD): Flow follows steepest descent path
 * - Multiple Flow Direction (MFD): Flow is distributed proportionally to all
 * downslope neighbors
 *
 * Key Features:
 * - Transient simulation with explicit time stepping
 * - Priority flooding to handle depressions
 * - Mass conservation through continuity equation
 * - Manning's friction for discharge calculations
 * - Flexible boundary conditions
 */

// ============================================================================
// UTILITY FUNCTIONS: Min/Max operations for floating point values
// ============================================================================

// Helper functions to mimic C++ std::min and std::max
inline GF_FLOAT max_float(GF_FLOAT a, GF_FLOAT b) { return (a > b) ? a : b; }
inline GF_FLOAT min_float(GF_FLOAT a, GF_FLOAT b) { return (a < b) ? a : b; }

// ============================================================================
// SINGLE FLOW DIRECTION IMPLEMENTATION
// ============================================================================

/*
 * Internal function implementing GraphFlood with Single Flow Direction (SFD)
 *
 * In SFD mode:
 * - Each cell directs all its flow to exactly one downstream neighbor
 * - Flow direction is determined by steepest hydraulic gradient
 * - More computationally efficient but less realistic for divergent flow areas
 * - Better suited for channelized flow systems
 *
 * Algorithm:
 * 1. Build flow graph using priority flooding
 * 2. Accumulate flow following single steepest paths
 * 3. Calculate discharge using Manning's equation
 * 4. Update water depths using continuity equation
 * 5. Repeat for specified number of iterations
 */
void _graphflood_full_sfd(
    GF_FLOAT* Z,               // Digital elevation model [input/const]
    GF_FLOAT* hw,              // Water depth array [input/output]
    uint8_t* BCs,              // Boundary conditions [input/const]
    GF_FLOAT* Precipitations,  // Precipitation rates [input/const]
    GF_FLOAT* manning,         // Manning's roughness [input/const]
    GF_UINT* dim,              // Grid dimensions [nx, ny] [input/const]
    GF_FLOAT dt,               // Time step size [input]
    GF_FLOAT dx,               // Grid spacing [input]
    bool SFD,                  // Flow direction flag [input]
    bool D8,                   // Connectivity flag [input]
    GF_UINT N_iterations,      // Number of iterations [input]
    GF_FLOAT step)             // Flooding step size [input]
{
  // --------------------------------------------------------------------------
  // MEMORY ALLOCATION: Create working arrays for SFD algorithm
  // --------------------------------------------------------------------------

  // Hydraulic surface elevation (ground + water)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i)
    Zw[i] = Z[i] + hw[i];  // Water surface = elevation + depth

  // Flow graph data structures for single flow direction
  GF_UINT* Sreceivers =
      (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));  // Single receiver per cell
  GF_FLOAT* distToReceivers =
      (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));  // Distance to receiver
  GF_UINT* Sdonors = (GF_UINT*)malloc(
      sizeof(GF_UINT) * nxy(dim) * (D8 ? 8 : 4));  // Multiple donors per cell
  uint8_t* NSdonors = (uint8_t*)malloc(sizeof(uint8_t) *
                                       nxy(dim));  // Number of donors per cell
  GF_UINT* Stack =
      (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));  // Processing order stack
  GF_FLOAT* Qwin = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) *
                                     nxy(dim));  // Input discharge per cell

  // Cell area for volume calculations
  GF_FLOAT cell_area = dx * dx;

  // --------------------------------------------------------------------------
  // MAIN ITERATION LOOP: Time stepping for transient simulation
  // --------------------------------------------------------------------------

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // ------------------------------------------------------------------------
    // STEP 1: BUILD FLOW GRAPH with priority flooding
    // ------------------------------------------------------------------------
    /*
     * Priority flooding ensures:
     * - All depressions in hydraulic surface are filled
     * - Flow paths are monotonically decreasing (no local minima)
     * - Single receiver is identified for each cell
     * - Topological ordering is established
     *
     * This step is crucial for numerical stability and physical realism
     */
    compute_sfgraph_priority_flood(Zw, Sreceivers, distToReceivers, Sdonors,
                                   NSdonors, Stack, BCs, dim, dx, D8, step);

    // ------------------------------------------------------------------------
    // STEP 2: FLOW ACCUMULATION following single flow paths
    // ------------------------------------------------------------------------
    /*
     * Accumulate discharge from precipitation and upstream contributions
     * Following single steepest descent paths established by flow graph
     */
    compute_weighted_drainage_area_single_flow(Qwin, Precipitations, Sreceivers,
                                               Stack, dim, dx);

    // ------------------------------------------------------------------------
    // STEP 3: UPDATE WATER DEPTHS using continuity equation
    // ------------------------------------------------------------------------
    /*
     * Process cells in reverse topological order (upstream to downstream)
     * This ensures hydraulic gradients remain explicit during updates
     */
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      // Get current cell (reverse stack order for proper upstream-downstream
      // processing)
      GF_UINT node = Stack[nxy(dim) - i - 1];
      GF_UINT rec = Sreceivers[node];

      // Skip cells that flow to themselves (pits/outlets)
      if (rec == node) continue;

      // Skip boundary cells that can discharge out of domain
      if (can_out(node, BCs)) continue;

      // Skip dry cells with no input (computational efficiency)
      if (Zw[node] == Z[node] && Qwin[node] == 0) continue;

      // ----------------------------------------------------------------------
      // HYDRAULIC GRADIENT CALCULATION
      // ----------------------------------------------------------------------

      // Calculate water surface slope to receiver
      // Use minimum slope to ensure numerical stability
      GF_FLOAT tSw =
          min_float(Zw[node] - Zw[rec], (GF_FLOAT)1e-6) / distToReceivers[node];

      // ----------------------------------------------------------------------
      // MANNING'S DISCHARGE CALCULATION
      // ----------------------------------------------------------------------

      /*
       * Manning's equation for open channel flow:
       * Q = (1/n) × A × R^(2/3) × S^(1/2)
       *
       * For 2D overland flow approximation:
       * Q ≈ (width/n) × depth^(5/3) × slope^(1/2)
       */
      GF_FLOAT tQwout = 0.0;
      if (Zw[node] > Z[node]) {  // Only if water is present
        GF_FLOAT depth = Zw[node] - Z[node];
        tQwout = (GF_FLOAT)(distToReceivers[node] / manning[node] *
                            pow(depth, 5.0 / 3.0) * sqrt(tSw));
      }

      // ----------------------------------------------------------------------
      // CONTINUITY EQUATION: Update water surface elevation
      // ----------------------------------------------------------------------

      /*
       * Explicit finite difference for continuity equation:
       * ∂h/∂t = (Q_in - Q_out) / A_cell
       *
       * h^(n+1) = h^n + dt × (Q_in - Q_out) / A_cell
       *
       * Ensures water surface never goes below ground level
       */
      Zw[node] =
          max_float(Z[node], Zw[node] + dt * (Qwin[node] - tQwout) / cell_area);
    }
  }

  // --------------------------------------------------------------------------
  // FINALIZATION: Convert hydraulic surface back to water depth
  // --------------------------------------------------------------------------

  // Extract water depths from hydraulic surface
  for (GF_UINT i = 0; i < nxy(dim); ++i)
    hw[i] = max_float(0.0, Zw[i] - Z[i]);  // Ensure non-negative depths

  // Free allocated memory
  free(Zw);
  free(Qwin);
  free(Sreceivers);
  free(distToReceivers);
  free(Sdonors);
  free(NSdonors);
  free(Stack);
}

// ============================================================================
// MULTIPLE FLOW DIRECTION IMPLEMENTATION
// ============================================================================

/*
 * Internal function implementing GraphFlood with Multiple Flow Direction (MFD)
 *
 * In MFD mode:
 * - Each cell can direct flow to multiple downstream neighbors
 * - Flow is distributed proportionally based on hydraulic gradients
 * - More physically realistic for divergent flow areas (alluvial fans,
 * floodplains)
 * - Computationally more intensive due to flow splitting calculations
 *
 * Algorithm:
 * 1. Apply priority flooding to hydraulic surface
 * 2. For each cell, calculate gradients to all valid neighbors
 * 3. Distribute flow proportionally based on gradient weights
 * 4. Calculate discharge using Manning's equation along steepest path
 * 5. Update water depths using continuity equation
 * 6. Repeat for specified iterations
 */
void _graphflood_full_mfd(
    GF_FLOAT* Z,               // Digital elevation model [input/const]
    GF_FLOAT* hw,              // Water depth array [input/output]
    uint8_t* BCs,              // Boundary conditions [input/const]
    GF_FLOAT* Precipitations,  // Precipitation rates [input/const]
    GF_FLOAT* manning,         // Manning's roughness [input/const]
    GF_UINT* dim,              // Grid dimensions [nx, ny] [input/const]
    GF_FLOAT dt,               // Time step size [input]
    GF_FLOAT dx,               // Grid spacing [input]
    bool SFD,                  // Flow direction flag [input]
    bool D8,                   // Connectivity flag [input]
    GF_UINT N_iterations,      // Number of iterations [input]
    GF_FLOAT step)             // Flooding step size [input]
{
  // --------------------------------------------------------------------------
  // NEIGHBOR CONNECTIVITY SETUP
  // --------------------------------------------------------------------------

  // Initialize neighbor offset arrays for grid traversal
  GF_INT offset[8];
  (D8 == false)
      ? generate_offset_D4_flat(offset, dim)  // 4-connectivity (cardinal)
      : generate_offset_D8_flat(offset,
                                dim);  // 8-connectivity (cardinal + diagonal)

  // Initialize distance arrays for each neighbor direction
  GF_FLOAT offdx[8];
  (D8 == false)
      ? generate_offsetdx_D4(offdx, dx)   // Cardinal distances only
      : generate_offsetdx_D8(offdx, dx);  // Cardinal + diagonal distances

  // Diagonal distance and cell area for calculations
  GF_FLOAT dxy = (GF_FLOAT)sqrt(2) * dx;  // Diagonal distance
  GF_FLOAT cell_area = dx * dx;           // Cell area for volume calculations

  // --------------------------------------------------------------------------
  // MEMORY ALLOCATION: Create working arrays for MFD algorithm
  // --------------------------------------------------------------------------

  // Hydraulic surface elevation (ground + water)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i) Zw[i] = Z[i] + hw[i];

  GF_FLOAT* Qwin =
      (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));  // Input discharge
  GF_FLOAT* Qwout =
      (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));  // Output discharge
  GF_UINT* Stack =
      (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));  // Processing order

  // Initialize arrays
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    Qwin[i] = 0.0;   // Start with zero input discharge
    Qwout[i] = 0.0;  // Start with zero output discharge
    Stack[i] = i;    // Initial stack order
  }

  // Flow weight array for multiple flow direction calculations
  GF_FLOAT weights[8];

  // --------------------------------------------------------------------------
  // MAIN ITERATION LOOP: Time stepping for transient simulation
  // --------------------------------------------------------------------------

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // ------------------------------------------------------------------------
    // STEP 1: PRIORITY FLOODING and topological ordering
    // ------------------------------------------------------------------------
    /*
     * Priority flooding fills depressions and establishes flow paths
     * Topological ordering ensures proper upstream-downstream processing
     */
    compute_priority_flood_plus_topological_ordering(Zw, Stack, BCs, dim, D8,
                                                     step);

    // Reset discharge arrays for this iteration
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      Qwin[i] = 0.0;
      Qwout[i] = 0.0;
    }

    // ------------------------------------------------------------------------
    // STEP 2: FLOW ACCUMULATION with multiple flow direction
    // ------------------------------------------------------------------------
    /*
     * Process cells in reverse topological order (upstream to downstream)
     * This maintains explicit hydraulic gradients during flow distribution
     */
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      // Get current cell index (reverse order processing)
      GF_UINT node = Stack[nxy(dim) - i - 1];

      // Skip invalid cells
      if (is_nodata(node, BCs)) continue;

      // Skip boundary cells that discharge out of domain
      if (can_out(node, BCs)) continue;

      // ----------------------------------------------------------------------
      // LOCAL PRECIPITATION INPUT
      // ----------------------------------------------------------------------

      // Add precipitation contribution (rate × cell area = volume/time)
      Qwin[node] += Precipitations[node] * dx * dx;

      // ----------------------------------------------------------------------
      // GRADIENT ANALYSIS: Calculate slopes to all valid neighbors
      // ----------------------------------------------------------------------

      GF_FLOAT sumslope = 0.0;  // Total weighted slope for flow partitioning
      GF_FLOAT maxslope = 0.0;  // Steepest slope for Manning's calculation
      GF_FLOAT dxmaxdir = dx;   // Distance to steepest neighbor

      // Examine all potential neighbors
      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        // Check if neighbor is valid and within bounds
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) {
          weights[n] = 0;
          continue;
        }

        GF_UINT nnode = node + offset[n];  // Neighbor cell index

        // Skip uphill neighbors or flow-restricted cells
        if (Zw[nnode] >= Zw[node] ||             // Neighbor is uphill
            can_receive(nnode, BCs) == false ||  // Neighbor cannot receive flow
            can_give(node, BCs) == false) {  // Current cell cannot give flow
          weights[n] = 0;
          continue;
        }

        // Calculate hydraulic gradient (slope) to this neighbor
        // Ensure minimum slope to prevent numerical issues
        GF_FLOAT tSw =
            max_float((GF_FLOAT)1e-8, (Zw[node] - Zw[nnode]) / offdx[n]);

        // Calculate flow weight proportional to gradient and effective flow
        // width Use appropriate distance based on flow direction (cardinal vs
        // diagonal)
        weights[n] = tSw * ((dx == offdx[n] || D8 == false) ? dx : dxy);

        // Accumulate total slope for proportional distribution
        sumslope += weights[n];

        // Track steepest gradient for discharge calculations
        if (tSw > maxslope) {
          maxslope = tSw;
          dxmaxdir = offdx[n];
        }
      }

      // ----------------------------------------------------------------------
      // FLOW DISTRIBUTION: Partition flow among downstream neighbors
      // ----------------------------------------------------------------------

      /*
       * Multiple Flow Direction (MFD) algorithm:
       * Each downstream neighbor receives flow proportional to its gradient
       * weight Flow fraction = (neighbor_weight / total_weights)
       *
       * This approach naturally handles:
       * - Flow convergence (multiple inputs to one cell)
       * - Flow divergence (one cell distributing to multiple neighbors)
       * - Flow routing through complex topography
       */
      if (sumslope > 0) {
        for (GF_UINT n = 0; n < N_neighbour(D8); ++n) {
          if (weights[n] == 0) continue;  // Skip neighbors with no flow

          // Distribute flow proportionally to gradient weight
          Qwin[node + offset[n]] += weights[n] / sumslope * Qwin[node];
        }
      }

      // ----------------------------------------------------------------------
      // MANNING'S DISCHARGE CALCULATION
      // ----------------------------------------------------------------------

      /*
       * Calculate discharge using Manning's equation along steepest gradient
       * This represents the maximum conveyance capacity of the cell
       */
      if (Zw[node] > Z[node]) {  // Only if water is present
        GF_FLOAT depth = Zw[node] - Z[node];
        Qwout[node] = (GF_FLOAT)(dxmaxdir / manning[node] *
                                 pow(depth, 5.0 / 3.0) * sqrt(maxslope));
      }
    }

    // ------------------------------------------------------------------------
    // STEP 3: UPDATE WATER DEPTHS using continuity equation
    // ------------------------------------------------------------------------

    for (GF_UINT node = 0; node < nxy(dim); ++node) {
      // ----------------------------------------------------------------------
      // CONTINUITY EQUATION: Update hydraulic surface
      // ----------------------------------------------------------------------

      /*
       * Explicit time stepping for continuity equation:
       *
       * ∂(h)/∂t + ∇·(q) = P
       *
       * Where:
       * - h = water depth
       * - q = discharge per unit width
       * - P = precipitation rate
       *
       * Discretized as:
       * h^(n+1) = h^n + dt × (Q_in - Q_out) / A_cell
       *
       * Constraint: Water surface cannot go below ground level
       */
      Zw[node] = max_float(
          Z[node], Zw[node] + dt * (Qwin[node] - Qwout[node]) / cell_area);
    }
  }

  // --------------------------------------------------------------------------
  // FINALIZATION: Convert hydraulic surface back to water depths
  // --------------------------------------------------------------------------

  // Extract water depths, ensuring non-negative values
  for (GF_UINT i = 0; i < nxy(dim); ++i) hw[i] = max_float(0.0, Zw[i] - Z[i]);

  // Free allocated memory
  free(Zw);
  free(Qwin);
  free(Qwout);
  free(Stack);
}

// ============================================================================
// PUBLIC API FUNCTION: Main GraphFlood interface
// ============================================================================

/*
 * GRAPHFLOOD_FULL: Complete transient flood simulation
 *
 * This is the main public interface for running complete flood simulations.
 * It supports both steady-state and transient modeling with flexible
 * flow direction schemes and boundary conditions.
 *
 * PARAMETERS:
 * -----------
 * Z            : Digital elevation model [m] (input, constant)
 * hw           : Water depth at each cell [m] (input/output, modified)
 * BCs          : Boundary condition flags (input, constant)
 * Precipitations: Precipitation rate per cell [m/s] (input, constant)
 * manning      : Manning's roughness coefficient [s/m^(1/3)] (input, constant)
 * dim          : Grid dimensions [nx, ny] (input, constant)
 * dt           : Time step size [s] (input, constant)
 * dx           : Grid cell spacing [m] (input, constant)
 * SFD          : Single Flow Direction flag (true=SFD, false=MFD) (input)
 * D8           : Connectivity scheme (true=8-neighbor, false=4-neighbor)
 * (input) N_iterations : Number of time steps to simulate (input) step : Step
 * size for priority flooding algorithm (input)
 *
 * ALGORITHM SELECTION:
 * -------------------
 * The function automatically selects the appropriate algorithm based on the SFD
 * flag:
 *
 * SFD=true:  Single Flow Direction
 *            - Faster computation
 *            - Good for channelized flow
 *            - Each cell flows to exactly one neighbor
 *
 * SFD=false: Multiple Flow Direction
 *            - More realistic for overland flow
 *            - Better for complex topography
 *            - Flow can be distributed among multiple neighbors
 *
 * NUMERICAL CONSIDERATIONS:
 * ------------------------
 * - Time step (dt) should satisfy CFL condition for numerical stability
 * - Grid resolution (dx) affects computational cost and accuracy
 * - Number of iterations determines simulation duration
 * - Priority flooding ensures numerical stability in depressions
 *
 * APPLICATIONS:
 * ------------
 * - Urban flood modeling
 * - Dam break simulations
 * - Overland flow routing
 * - Watershed hydrology
 * - Emergency response planning
 */
TOPOTOOLBOX_API
void graphflood_full(GF_FLOAT* Z,   // Digital elevation model [input]
                     GF_FLOAT* hw,  // Water depths [input/output]
                     uint8_t* BCs,  // Boundary conditions [input]
                     GF_FLOAT* Precipitations,  // Precipitation rates [input]
                     GF_FLOAT* manning,         // Manning's roughness [input]
                     GF_UINT* dim,              // Grid dimensions [input]
                     GF_FLOAT dt,               // Time step size [input]
                     GF_FLOAT dx,               // Grid spacing [input]
                     bool SFD,                  // Flow direction scheme [input]
                     bool D8,                   // Connectivity scheme [input]
                     GF_UINT N_iterations,      // Number of iterations [input]
                     GF_FLOAT step)             // Flooding step size [input]
{
  // Route to appropriate algorithm based on flow direction scheme
  if (SFD) {
    // Single Flow Direction: computationally efficient
    _graphflood_full_sfd(Z, hw, BCs, Precipitations, manning, dim, dt, dx, SFD,
                         D8, N_iterations, step);
  } else {
    // Multiple Flow Direction: more physically realistic
    _graphflood_full_mfd(Z, hw, BCs, Precipitations, manning, dim, dt, dx, SFD,
                         D8, N_iterations, step);
  }
}

/*
 * GRAPHFLOOD METRICS FUNCTION
 *
 * This function implements a hydrological flow routing algorithm using priority
 * flooding and topological ordering to calculate various discharge metrics for
 * flood modeling.
 *
 * The algorithm performs:
 * 1. Priority flooding to determine hydraulic surfaces and flow directions
 * 2. Topological ordering to ensure upstream-to-downstream processing
 * 3. Flow accumulation based on precipitation inputs
 * 4. Discharge calculations using Manning's equation
 *
 * Key metrics calculated:
 * - Qi: Stationary discharge based on drainage area accumulation
 * - Qo: Calculated discharge from the model (should equal Qi at convergence)
 * - q: Discharge per unit width
 * - u: Flow velocity
 * - Sw: Water surface slope
 */

TOPOTOOLBOX_API
void graphflood_metrics(
    GF_FLOAT* Z,               // Digital elevation model (DEM) [input]
    GF_FLOAT* hw,              // Water depth at each cell [input]
    uint8_t* BCs,              // Boundary conditions array [input]
    GF_FLOAT* Precipitations,  // Precipitation rate per cell [input]
    GF_FLOAT* manning,         // Manning's roughness coefficient [input]
    GF_UINT* dim,              // Grid dimensions [nx, ny] [input]
    GF_FLOAT dx,               // Grid cell size (spacing) [input]
    bool D8,        // Flow direction scheme: true=D8, false=D4 [input]
    GF_FLOAT step,  // Step size for flooding algorithm [input]
    GF_FLOAT* Qi,   // Input discharge (accumulated flow) [output]
    GF_FLOAT* Qo,   // Output discharge (Manning's equation) [output]
    GF_FLOAT* qo,   // Discharge per unit width [output]
    GF_FLOAT* u,    // Flow velocity [output]
    GF_FLOAT* Sw)   // Water surface slope [output]
{
  // ============================================================================
  // INITIALIZATION: Set up neighbor connectivity and distance arrays
  // ============================================================================

  // Initialize offset array for neighboring cell operations
  // D4: 4 cardinal directions (N,S,E,W), D8: 8 directions (including diagonals)
  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)   // 4-connectivity
                : generate_offset_D8_flat(offset, dim);  // 8-connectivity

  // Initialize distance to each neighbor
  // Cardinal directions = dx, diagonal directions = dx*sqrt(2)
  GF_FLOAT offdx[8];
  (D8 == false)
      ? generate_offsetdx_D4(offdx, dx)   // Only cardinal distances
      : generate_offsetdx_D8(offdx, dx);  // Cardinal + diagonal distances

  // Diagonal distance for flow calculations
  GF_FLOAT dxy = (GF_FLOAT)sqrt(2) * dx;

  // ============================================================================
  // MEMORY ALLOCATION: Create working arrays
  // ============================================================================

  // Create hydraulic surface array (elevation + water depth)
  // Zw represents the water surface elevation at each cell
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i)
    Zw[i] = Z[i] + hw[i];  // Water surface = ground elevation + water depth

  // Stack for topologically ordered processing of cells
  // Will contain cell indices in processing order (upstream to downstream)
  GF_UINT* Stack = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));

  // Initialize input discharge array to zero and populate stack with cell
  // indices
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    Qi[i] = 0.0;   // Start with no accumulated flow
    Stack[i] = i;  // Initially, stack contains all cells in order
  }

  // Array to store flow weights for multi-direction flow distribution
  GF_FLOAT weights[8];

  // ============================================================================
  // PRIORITY FLOODING: Determine flow directions and processing order
  // ============================================================================

  /*
   * Priority flooding algorithm:
   * 1. Fills depressions in the hydraulic surface (Zw)
   * 2. Ensures monotonic flow paths (no local minima except outlets)
   * 3. Creates topological ordering in Stack for upstream-to-downstream
   * processing
   *
   * This is critical because it ensures that when we process flow accumulation,
   * all upstream cells are processed before downstream cells.
   */
  compute_priority_flood_plus_topological_ordering(Zw, Stack, BCs, dim, D8,
                                                   step);

  // Re-initialize input discharge array after flooding
  // (flooding may have modified the processing order)
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    Qi[i] = 0.0;
  }

  // ============================================================================
  // FLOW ACCUMULATION: Process cells in topological order
  // ============================================================================

  /*
   * Process nodes from downstream to upstream (reverse stack order)
   * This ensures that receivers are never processed before donors,
   * maintaining explicit hydraulic slopes throughout the calculation
   */
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    // Get current cell index (processing in reverse topological order)
    GF_UINT node = Stack[nxy(dim) - i - 1];

    // Skip cells with no data
    if (is_nodata(node, BCs)) continue;

    // Skip boundary cells that can flow out of the domain
    // These cells don't accumulate flow from upstream
    if (can_out(node, BCs)) continue;

    // ------------------------------------------------------------------------
    // PRECIPITATION INPUT: Add local water input
    // ------------------------------------------------------------------------

    // Add precipitation contribution to input discharge
    // Precipitation rate × cell area gives volumetric flow rate
    Qi[node] += Precipitations[node] * dx * dx;

    // ------------------------------------------------------------------------
    // FLOW DIRECTION ANALYSIS: Calculate gradients to all neighbors
    // ------------------------------------------------------------------------

    GF_FLOAT sumslope = 0.0;  // Sum of all downslope gradients (for weighting)
    GF_FLOAT maxslope = 0.0;  // Steepest gradient (for Manning's equation)
    GF_FLOAT dxmaxdir = dx;   // Distance to steepest neighbor

    // Examine all possible neighbors
    for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
      // Check if neighbor is within grid bounds and valid
      if (check_bound_neighbour(node, n, dim, BCs, D8) == false) {
        weights[n] = 0;  // No flow to invalid neighbors
        continue;
      }

      GF_UINT nnode = node + offset[n];  // Calculate neighbor index

      // Skip uphill neighbors or neighbors that can't receive/give flow
      if (Zw[nnode] >= Zw[node] ||             // Uphill neighbor
          can_receive(nnode, BCs) == false ||  // Neighbor can't receive flow
          can_give(node, BCs) == false) {      // Current cell can't give flow
        weights[n] = 0;
        continue;
      }

      // Calculate hydraulic gradient (slope) to this neighbor
      // Ensure minimum slope to avoid division by zero
      GF_FLOAT tSw =
          max_float((GF_FLOAT)1e-8, (Zw[node] - Zw[nnode]) / offdx[n]);

      // Calculate flow weight proportional to gradient and effective width
      // For diagonal flow, use diagonal distance; for cardinal, use dx
      weights[n] = tSw * ((dx == offdx[n] || D8 == false) ? dx : dxy);

      // Accumulate total slope for proportional flow distribution
      sumslope += weights[n];

      // Track steepest gradient for velocity calculations
      if (tSw > maxslope) {
        maxslope = tSw;
        dxmaxdir = offdx[n];  // Distance to steepest neighbor
      }
    }

    // ------------------------------------------------------------------------
    // FLOW DISTRIBUTION: Transfer accumulated flow to downstream neighbors
    // ------------------------------------------------------------------------

    /*
     * Distribute the accumulated flow (Qi) to downstream neighbors
     * proportionally based on their hydraulic gradients.
     * This implements multiple flow direction (MFD) routing.
     */
    if (sumslope > 0) {
      for (GF_UINT n = 0; n < N_neighbour(D8); ++n) {
        if (weights[n] == 0) continue;  // Skip neighbors with no flow

        // Transfer flow proportional to gradient weight
        // Each downstream neighbor receives: (its_weight / total_weight) ×
        // total_flow
        Qi[node + offset[n]] += weights[n] / sumslope * Qi[node];
      }
    }

    // ------------------------------------------------------------------------
    // DISCHARGE CALCULATIONS: Apply Manning's equation for flow velocity
    // ------------------------------------------------------------------------

    /*
     * Calculate volumetric discharge using Manning's friction equation:
     * Q = (1/n) × A × R^(2/3) × S^(1/2)
     *
     * Where:
     * - n = Manning's roughness coefficient
     * - A = Cross-sectional area (approximated as width × depth)
     * - R = Hydraulic radius (approximated as depth for wide channels)
     * - S = Hydraulic slope (water surface gradient)
     *
     * Simplified for 2D flow: Q = (width/n) × depth^(5/3) × slope^(1/2)
     */
    GF_FLOAT tQwout = 0.0;

    // Only calculate discharge if there's water depth
    if (Zw[node] > Z[node]) {
      GF_FLOAT depth = Zw[node] - Z[node];  // Water depth

      // Manning's equation: Q = (width/n) × depth^(5/3) × slope^(1/2)
      tQwout = (GF_FLOAT)(dxmaxdir / manning[node] * pow(depth, 5.0 / 3.0) *
                          sqrt(maxslope));
    }

    // ------------------------------------------------------------------------
    // OUTPUT METRICS: Store calculated values
    // ------------------------------------------------------------------------

    Qo[node] = tQwout;             // Volumetric discharge [m³/s]
    qo[node] = tQwout / dxmaxdir;  // Discharge per unit width [m²/s]

    // Calculate flow velocity: v = q/h (discharge per unit width ÷ depth)
    if (Zw[node] - Z[node] > 0) {
      u[node] = qo[node] / (Zw[node] - Z[node]);  // Flow velocity [m/s]
    } else {
      u[node] = 0.0;  // No velocity if no water depth
    }

    Sw[node] = maxslope;  // Water surface slope [-]
  }

  // Free allocated memory
  free(Zw);
  free(Stack);
}

// ============================================================================
// DYNAMIC INDUCED GRAPH IMPLEMENTATION
// ============================================================================

/*
 * GRAPHFLOOD_DYNAMIC_GRAPH: Dynamic induced graph flood simulation
 *
 * This function implements a wavefront-based flood routing approach that
 * processes cells in descending elevation order, starting from specified
 * input discharge locations. Unlike graphflood_full which processes the
 * entire grid, this method dynamically builds a flow graph by propagating
 * from high to low elevations.
 *
 * Key Differences from graphflood_full:
 * - Uses max-heap to process in descending hydraulic elevation order
 * - Starts from cells with input_Qw > 0
 * - Dynamically expands processing to downstream neighbors
 * - Upstream flow accumulation calculated on-the-fly
 *
 * FLOW ACCUMULATION ALGORITHM:
 * ----------------------------
 * The algorithm ensures mass conservation by tracking discharge (Qw) through
 * the priority queue structure:
 *
 * 1. Qw Transport: Each cell in the PQ carries a Qw value representing
 *    discharge being transported to that location
 *
 * 2. Flow Accumulation at each cell:
 *    total_Qw = Qw_transported (from PQ cell)
 *             + extra_Qw[node] (flow added by upstream neighbors already in PQ)
 *             + input_Qw[node] (only on first visit)
 *             + Precipitations[node] * cell_area (only on first visit)
 *
 * 3. Qwin Array: Burns the maximum total_Qw seen at each cell across all visits
 *    This handles cells visited multiple times due to local minima
 *
 * 4. Flow Distribution Rules:
 *    a) If ANY downstream neighbor can_out: give ALL flow to that neighbor
 *    b) Else if ALL downstream neighbors visited: give ALL to steepest
 *    c) Else: split proportionally to NON-VISITED downstream neighbors
 *
 * 5. Flow Transmission:
 *    - If neighbor NOT in PQ: push with Qw value directly
 *    - If neighbor IN PQ: add to extra_Qw buffer (consumed when popped)
 *
 * This approach prevents over-accumulation while handling revisits correctly.
 *
 * PARAMETERS:
 * -----------
 * Z              : Digital elevation model [m] (input, constant)
 * hw             : Water depth at each cell [m] (input/output, modified)
 * BCs            : Boundary condition flags (input, constant)
 * Precipitations : Precipitation rate per cell [m/s] (input, constant)
 * manning        : Manning's roughness coefficient [s/m^(1/3)] (input,
 * constant) input_Qw       : Input discharge locations/values [m³/s] (input,
 * constant) Qwin           : Accumulated discharge array [m³/s] (input/output)
 *                  Reinitialized each iteration, stores max Qw seen at each
 * cell dim            : Grid dimensions [nx, ny] (input, constant) dt : Time
 * step size [s] (input, constant) dx             : Grid cell spacing [m]
 * (input, constant) D8             : Connectivity scheme (true=8-neighbor,
 * false=4-neighbor) N_iterations   : Number of time steps to simulate (input)
 */
TOPOTOOLBOX_API
void graphflood_dynamic_graph(
    GF_FLOAT* Z,               // Digital elevation model [input]
    GF_FLOAT* hw,              // Water depths [input/output]
    uint8_t* BCs,              // Boundary conditions [input]
    GF_FLOAT* Precipitations,  // Precipitation rates [input]
    GF_FLOAT* manning,         // Manning's roughness [input]
    GF_FLOAT* input_Qw,        // Input discharge locations [input]
    GF_FLOAT* Qwin,            // Accumulated input discharge [input/output]
    GF_UINT* dim,              // Grid dimensions [input]
    GF_FLOAT dt,               // Time step size [input]
    GF_FLOAT dx,               // Grid spacing [input]
    bool D8,                   // Connectivity scheme [input]
    GF_UINT N_iterations)      // Number of iterations [input]
{
  // --------------------------------------------------------------------------
  // NEIGHBOR CONNECTIVITY SETUP
  // --------------------------------------------------------------------------
  // Generate offset arrays for neighbor access based on connectivity (D4/D8)
  // offset[n] gives the flat index offset to reach neighbor n
  // offdx[n] gives the distance to neighbor n (dx for cardinal, dx*sqrt(2) for
  // diagonal)

  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)
                : generate_offset_D8_flat(offset, dim);

  GF_FLOAT offdx[8];
  (D8 == false) ? generate_offsetdx_D4(offdx, dx)
                : generate_offsetdx_D8(offdx, dx);

  GF_FLOAT cell_area = dx * dx;  // Area of each grid cell [m²]

  // --------------------------------------------------------------------------
  // MEMORY ALLOCATION
  // --------------------------------------------------------------------------

  GF_UINT tnxy = nxy(dim);  // Total number of cells in the grid

  // Hydraulic surface elevation = ground elevation + water depth [m]
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);
  for (GF_UINT i = 0; i < tnxy; ++i) {
    Zw[i] = Z[i] + hw[i];
  }

  // Flow arrays:
  // - Qwout: Output discharge from each cell calculated via Manning's equation
  // [m³/s]
  // - extra_Qw: Buffer for flow from upstream neighbors already in PQ [m³/s]
  //   This prevents double-counting when cells are pushed multiple times
  GF_FLOAT* Qwout = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);
  GF_FLOAT* extra_Qw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);

  // State tracking arrays:
  // - visited: Has this cell been processed before?
  // - inPQ: Is this cell currently in the priority queue?
  uint8_t* visited = (uint8_t*)malloc(sizeof(uint8_t) * tnxy);
  uint8_t* inPQ = (uint8_t*)malloc(sizeof(uint8_t) * tnxy);

  // --------------------------------------------------------------------------
  // IDENTIFY INPUT CELLS
  // --------------------------------------------------------------------------
  // Pre-compute list of cells with input discharge (input_Qw > 0)
  // These cells serve as entry points for the wavefront propagation

  GF_UINT n_input_cells = 0;
  for (GF_UINT i = 0; i < tnxy; ++i) {
    if (input_Qw[i] > 0.0) {
      n_input_cells++;
    }
  }

  // Store indices of input cells for fast iteration
  GF_UINT* input_indices =
      (GF_UINT*)malloc(sizeof(GF_UINT) * (n_input_cells + 1));
  GF_UINT idx = 0;
  for (GF_UINT i = 0; i < tnxy; ++i) {
    if (input_Qw[i] > 0.0) {
      input_indices[idx++] = i;
    }
  }

  // --------------------------------------------------------------------------
  // MAIN ITERATION LOOP
  // --------------------------------------------------------------------------
  // Each iteration represents one time step (dt) of the simulation
  // The priority queue is rebuilt each iteration from input cells

  MaxHeapPQueue pq;
  maxheap_init(&pq, tnxy);

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // Reset arrays for this iteration
    for (GF_UINT i = 0; i < tnxy; ++i) {
      Qwin[i] = 0.0;       // Accumulated discharge (will track maximum)
      Qwout[i] = 0.0;      // Output discharge via Manning
      extra_Qw[i] = 0.0;   // Flow buffer for cells already in PQ
      visited[i] = false;  // Clear visitation flags
      inPQ[i] = false;     // Clear PQ membership flags
    }

    // ------------------------------------------------------------------------
    // INITIALIZE PRIORITY QUEUE with input cells
    // ------------------------------------------------------------------------
    // Start wavefront from input cells with Qw = 0.0
    // They will pick up input_Qw when first popped

    for (GF_UINT i = 0; i < n_input_cells; ++i) {
      GF_UINT node = input_indices[i];
      maxheap_push_with_qw(&pq, node, Zw[node], 0.0);
      inPQ[node] = true;
    }

    // ------------------------------------------------------------------------
    // PROCESS CELLS IN DESCENDING ELEVATION ORDER
    // ------------------------------------------------------------------------
    // Max-heap ensures cells are processed from high to low elevation
    // This guarantees upstream cells are processed before downstream

    while (maxheap_empty(&pq) == false) {
      // ======================================================================
      // POP CELL FROM PRIORITY QUEUE
      // ======================================================================
      // Extract both the node index and the Qw being transported to this cell
      GF_FLOAT Qw_transported;
      GF_UINT node = maxheap_pop_and_get_key_qw(&pq, &Qw_transported);

      inPQ[node] = false;
      bool was_visited_before = visited[node];
      visited[node] = true;

      // Skip invalid cells (nodata)
      if (is_nodata(node, BCs)) continue;

      // ======================================================================
      // ACCUMULATE TOTAL DISCHARGE AT THIS CELL
      // ======================================================================
      // Combine all sources of discharge:
      // 1. Qw_transported: flow carried in the PQ cell
      // 2. extra_Qw: flow added by upstream neighbors already in PQ
      // 3. input_Qw: external input discharge (only first visit)
      // 4. Precipitation: local water input (only first visit)

      GF_FLOAT total_Qw = Qw_transported;

      // Always consume extra_Qw buffer (reset after consumption)
      total_Qw += extra_Qw[node];
      extra_Qw[node] = 0.0;

      // Add local sources only on first visit (prevents over-accumulation)
      if (!was_visited_before) {
        total_Qw += input_Qw[node] + Precipitations[node] * cell_area;
      }

      // Burn maximum Qw seen into Qwin array
      // Handles cells visited multiple times (e.g., in depressions)
      if (total_Qw > Qwin[node]) {
        Qwin[node] = total_Qw;
      }

      // Skip boundary cells that can flow out (flow exits domain here)
      if (can_out(node, BCs)) continue;

      // ======================================================================
      // DEPRESSION FILLING: Check for lower neighbors and raise if needed
      // ======================================================================
      // If cell has no downstream (lower) neighbors, it's in a depression
      // Raise it slightly above the lowest neighbor to create flow path

      GF_FLOAT min_neighbor_zw = INFINITY;
      bool has_lower_neighbor = false;
      bool has_any_neighbor = false;

      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) continue;

        GF_UINT nnode = node + offset[n];
        if (is_nodata(nnode, BCs)) continue;

        has_any_neighbor = true;

        if (Zw[nnode] < Zw[node]) {
          has_lower_neighbor = true;
        }
        if (Zw[nnode] < min_neighbor_zw) {
          min_neighbor_zw = Zw[nnode];
        }
      }

      // Raise cell slightly above lowest neighbor if trapped in depression
      if (has_any_neighbor && has_lower_neighbor == false) {
        Zw[node] = min_neighbor_zw + 1e-3;
        maxheap_push_with_qw(&pq, node, Zw[node], total_Qw);
        inPQ[node] = true;
        continue;
      }

      // ======================================================================
      // FLOW DISTRIBUTION TO DOWNSTREAM NEIGHBORS
      // ======================================================================
      // Strategy depends on downstream neighbor status:
      // 1. If any neighbor can_out: route ALL flow to boundary
      // 2. Else if all downstream visited: route ALL to steepest (backflow)
      // 3. Else: split proportionally to non-visited downstream neighbors

      // Check for boundary outflow neighbor first (highest priority)
      GF_UINT can_out_neighbor = node;
      bool has_can_out_neighbor = false;

      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) continue;

        GF_UINT nnode = node + offset[n];
        if (is_nodata(nnode, BCs)) continue;

        // Check if downstream (lower elevation)
        if (Zw[nnode] >= Zw[node]) continue;

        if (can_out(nnode, BCs)) {
          has_can_out_neighbor = true;
          can_out_neighbor = nnode;
          break;  // Found boundary - route all flow here
        }
      }

      // ----------------------------------------------------------------------
      // CASE 1: Boundary neighbor found - route ALL flow to domain boundary
      // ----------------------------------------------------------------------
      if (has_can_out_neighbor) {
        // If neighbor already in PQ: add to extra_Qw buffer
        // If not in PQ: push with Qw value directly
        if (inPQ[can_out_neighbor]) {
          extra_Qw[can_out_neighbor] += total_Qw;
        } else {
          maxheap_push_with_qw(&pq, can_out_neighbor, Zw[can_out_neighbor],
                               total_Qw);
          inPQ[can_out_neighbor] = true;
        }
      } else {
        // --------------------------------------------------------------------
        // CASE 2 & 3: No boundary neighbor - analyze downstream neighbors
        // --------------------------------------------------------------------
        // Calculate slopes to all downstream neighbors
        // Track both steepest (for backflow) and sum of non-visited (for
        // splitting)

        GF_FLOAT sum_slopes = 0.0;     // Sum of slopes to non-visited neighbors
        GF_FLOAT maxslope = 0.0;       // Steepest slope (any neighbor)
        GF_FLOAT dxmaxslope = 0.0;     // Distance to steepest neighbor
        GF_UINT steepest_node = node;  // Steepest downstream neighbor
        bool all_downstream_visited =
            true;  // Track if all downstream are visited

        for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
          if (check_bound_neighbour(node, n, dim, BCs, D8) == false) continue;

          GF_UINT nnode = node + offset[n];
          if (is_nodata(nnode, BCs)) continue;

          // Check if downstream (lower elevation)
          if (Zw[nnode] >= Zw[node]) continue;

          GF_FLOAT slope =
              max_float((GF_FLOAT)1e-8, (Zw[node] - Zw[nnode]) / offdx[n]);

          // Track steepest neighbor (for CASE 2: backflow routing)
          if (slope > maxslope) {
            maxslope = slope;
            dxmaxslope = offdx[n];
            steepest_node = nnode;
          }

          // Track non-visited neighbors (for CASE 3: proportional split)
          if (!visited[nnode]) {
            sum_slopes += slope;
            all_downstream_visited = false;
          }
        }

        // ----------------------------------------------------------------------
        // CASE 3: Some downstream neighbors not visited - split flow
        // ----------------------------------------------------------------------
        if (sum_slopes > 0.0) {
          // Distribute flow proportionally based on slopes to non-visited
          // neighbors This is the normal forward routing case
          for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
            if (check_bound_neighbour(node, n, dim, BCs, D8) == false) continue;

            GF_UINT nnode = node + offset[n];
            if (is_nodata(nnode, BCs)) continue;

            // Only consider downstream AND not yet visited neighbors
            if (Zw[nnode] >= Zw[node]) continue;
            if (visited[nnode]) continue;

            // Calculate proportion of flow going to this neighbor
            GF_FLOAT slope =
                max_float((GF_FLOAT)1e-8, (Zw[node] - Zw[nnode]) / offdx[n]);
            GF_FLOAT proportion = slope / sum_slopes;

            // Add to extra_Qw buffer (will be consumed when neighbor is popped)
            extra_Qw[nnode] += proportion * total_Qw;

            // Add neighbor to PQ if not already there (with Qw=0, will pick up
            // extra_Qw)
            if (!inPQ[nnode]) {
              maxheap_push_with_qw(&pq, nnode, Zw[nnode], 0.0);
              inPQ[nnode] = true;
            }
          }
        } else if (all_downstream_visited && steepest_node != node) {
          // --------------------------------------------------------------------
          // CASE 2: All downstream visited - backflow to steepest
          // --------------------------------------------------------------------
          // This handles cells in depressions where downstream was already
          // processed Route all flow to steepest neighbor to allow
          // re-processing
          if (inPQ[steepest_node]) {
            extra_Qw[steepest_node] += total_Qw;
          } else {
            maxheap_push_with_qw(&pq, steepest_node, Zw[steepest_node],
                                 total_Qw);
            inPQ[steepest_node] = true;
          }
        } else {
          // --------------------------------------------------------------------
          // CASE 4: No downstream path - stuck in pit
          // --------------------------------------------------------------------
          // Push cell back with raised elevation to escape depression
          maxheap_push_with_qw(&pq, node, Zw[node] + 1e-3, total_Qw);
          inPQ[node] = true;
          Zw[node] += 1e-3;
          continue;
        }

        // ======================================================================
        // CALCULATE OUTPUT DISCHARGE using Manning's equation
        // ======================================================================
        // Q = (1/n) * A * R^(2/3) * S^(1/2)
        // For wide shallow flow: R ≈ depth
        // Width ≈ dxmaxslope, Area = width * depth

        if (Zw[node] > Z[node] && dxmaxslope > 0. && maxslope > 0. &&
            !was_visited_before) {
          GF_FLOAT depth = max_float(Zw[node] - Z[node], 0.);
          Qwout[node] = (GF_FLOAT)(dxmaxslope / manning[node] *
                                   pow(depth, 5.0 / 3.0) * sqrt(maxslope));
        }
      }
    }  // End while (PQ not empty)

    // ========================================================================
    // UPDATE WATER DEPTHS
    // ========================================================================
    // Apply continuity equation: dh/dt = (Qin - Qout) / Area
    // Update hydraulic elevation: Zw = Zw + dt * (Qwin - Qwout) / cell_area
    // Only update cells that were active (had water or were visited)

    for (GF_UINT node = 0; node < tnxy; ++node) {
      // Only update cells with water or that were visited
      if (Zw[node] > Z[node] || visited[node]) {
        // Apply water balance: increase if Qwin > Qwout, decrease if Qwin <
        // Qwout Ensure Zw never drops below ground elevation Z
        Zw[node] = max_float(
            Z[node], Zw[node] + dt * (Qwin[node] - Qwout[node]) / cell_area);
      }
    }
  }  // End iteration loop

  // ==========================================================================
  // FINALIZATION
  // ==========================================================================
  // Convert hydraulic elevation back to water depth for output

  for (GF_UINT i = 0; i < tnxy; ++i) {
    hw[i] = max_float(0.0, Zw[i] - Z[i]);
  }

  // Clean up allocated memory
  maxheap_free(&pq);
  free(Zw);
  free(Qwout);
  free(extra_Qw);
  free(visited);
  free(inPQ);
  free(input_indices);
}

// ============================================================================
// COMPUTE INPUT DISCHARGE FROM DRAINAGE AREA THRESHOLD
// ============================================================================

/*
 * Computes input discharge array for dynamic graph based on drainage area
 * threshold. Identifies channel heads where drainage area transitions from
 * below to above the threshold and assigns precipitation-weighted discharge.
 *
 * Algorithm:
 * 1. Build single flow graph from hydraulic surface (Z + hw)
 * 2. Calculate both regular and precipitation-weighted drainage areas
 * 3. Identify transition points where:
 *    - Current node has drainage area < threshold
 *    - Receiver node has drainage area >= threshold
 * 4. Assign Qwin at these transition points as entry discharges
 */
TOPOTOOLBOX_API
void compute_input_Qw_from_area_threshold(GF_FLOAT* input_Qw, GF_FLOAT* Z,
                                          GF_FLOAT* hw, uint8_t* BCs,
                                          GF_FLOAT* Precipitations,
                                          GF_FLOAT area_threshold, GF_UINT* dim,
                                          GF_FLOAT dx, bool D8, GF_FLOAT step) {
  // --------------------------------------------------------------------------
  // MEMORY ALLOCATION
  // --------------------------------------------------------------------------

  GF_UINT tnxy = nxy(dim);

  // Hydraulic surface elevation (ground + water)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);
  for (GF_UINT i = 0; i < tnxy; ++i) {
    Zw[i] = Z[i] + hw[i];
  }

  // Flow graph data structures for single flow direction
  GF_UINT* Sreceivers = (GF_UINT*)malloc(sizeof(GF_UINT) * tnxy);
  GF_FLOAT* distToReceivers = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);
  GF_UINT* Sdonors = (GF_UINT*)malloc(sizeof(GF_UINT) * tnxy * (D8 ? 8 : 4));
  uint8_t* NSdonors = (uint8_t*)malloc(sizeof(uint8_t) * tnxy);
  GF_UINT* Stack = (GF_UINT*)malloc(sizeof(GF_UINT) * tnxy);

  // Drainage area arrays
  GF_FLOAT* drainage_area = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);
  GF_FLOAT* Qwin = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);

  // --------------------------------------------------------------------------
  // STEP 1: BUILD SINGLE FLOW GRAPH
  // --------------------------------------------------------------------------

  compute_sfgraph_priority_flood(Zw, Sreceivers, distToReceivers, Sdonors,
                                 NSdonors, Stack, BCs, dim, dx, D8, step);

  // --------------------------------------------------------------------------
  // STEP 2: CALCULATE DRAINAGE AREAS
  // --------------------------------------------------------------------------

  // Regular drainage area (unweighted)
  compute_drainage_area_single_flow(drainage_area, Sreceivers, Stack, dim, dx);

  // Precipitation-weighted drainage area
  compute_weighted_drainage_area_single_flow(Qwin, Precipitations, Sreceivers,
                                             Stack, dim, dx);

  // --------------------------------------------------------------------------
  // STEP 3: IDENTIFY CHANNEL HEADS AND ASSIGN INPUT DISCHARGE
  // --------------------------------------------------------------------------

  // Initialize input_Qw to zero
  for (GF_UINT i = 0; i < tnxy; ++i) {
    input_Qw[i] = 0.;
  }

  // Process nodes in reverse stack order (downstream to upstream)
  for (GF_UINT i = 0; i < tnxy; ++i) {
    GF_UINT node = Stack[tnxy - 1 - i];
    GF_UINT rec = Sreceivers[node];

    // Skip if node is its own receiver (outlet/pit)
    if (node == rec) continue;

    // Check if this is a channel head:
    // - Node has drainage area < threshold
    // - Receiver has drainage area >= threshold
    if (drainage_area[node] < area_threshold &&
        drainage_area[rec] >= area_threshold) {
      // This is an entry point - assign precipitation-weighted discharge
      input_Qw[node] = Qwin[node];
    }
  }

  // --------------------------------------------------------------------------
  // CLEANUP
  // --------------------------------------------------------------------------

  free(Zw);
  free(Sreceivers);
  free(distToReceivers);
  free(Sdonors);
  free(NSdonors);
  free(Stack);
  free(drainage_area);
  free(Qwin);
}
