# Water Balance and Phenology Simulation Pipeline

This project implements a simulation pipeline for calculating **crop phenology** and **soil water balance** using the **FAO-56** methodology. The core script is `01_Balance_agua.R`, which orchestrates the process from data loading to result generation.

## Overview

The pipeline fits into a larger workflow for agricultural risk analysis. It takes historical weather data and soil profiles to simulate crop development (Phenology) and daily water stress/usage (Water Balance) for defined management scenarios (planting dates, initial conditions).

Key components:
-   **Script:** `01_Balance_agua.R`
-   **Core Logic:** `lib/balance_fao56.R` (FAO-56 implementation) and `lib/fenologia.R` (Wheat phenology).
-   **Orchestration:** `lib/workers.R` (Functions dealing with iteration and parallelization preparations).

## Usage

The script is designed to be run via Rscript or interactively in an R environment. It accepts command-line arguments for configuration files.

```bash
Rscript 01_Balance_agua.R [config_file] [params_file] [api_config_file]
```

If no arguments are provided, it defaults to:
1.  `configuracion.yml`
2.  `parametros.yml`
3.  `configuracion_api.yml`

## Configuration

### `configuracion.yml`
Defines the directory structure and simulation settings.
*   `dir`: Paths for libraries, inputs, outputs, shared resources, weather data, and logs.
*   `soil_file`: Name of the DSSAT format soil file (e.g., `SOIL.SOL`).
*   `max.procesos`: Number of parallel processes to use.

### `parametros.yml`
Defines the simulation parameters.
*   `run_id`: Identifier for the current run (used for output folder naming).
*   `crops`: List of crops to simulate (e.g., `["WH"]` for Wheat).
*   `weather`: Start/End years and weather file name.
*   `parametros.fao`: Crop coefficients (Kc values, height) for the FAO model.
*   `initial_conditions`: List of initial soil water content fractions (e.g., `[0.2, 0.5, 1.0]`).
*   `soils`: A matrix defining combinations of soils, weather stations, varieties, and management dates to simulate.

## Detailed Execution Steps

The `01_Balance_agua.R` script executes a series of sequential steps, marked as "PASOS" in the code. Here is a detailed breakdown of what happens in each step:

### Step 1: Initialization (`PASO 1-3`)
*   **Environment Setup:** Clears memory, sets the Time Zone to UTC, and loads necessary R packages (`dplyr`, `caret`, `parallel`, etc.).
*   **Configuration Loading:** Reads the YAML configuration files to determine input/output paths and simulation parameters.
*   **Logging:** Initializes a logging system (using `futile.logger` via a custom `Script` class) to track execution progress and errors. Logs are saved to the `log/` directory.

### Step 2: Data Loading (`PASO 4-5`)
*   **Weather Data:** Reads the full historical weather dataset from the file specified in `parametros.yml`.
*   **Soil Data Prep:**
    *   Iterates through the requested soils.
    *   Calculates the **Maximum Depth** for each soil profile using `buscar_maxima_profundidad`.
    *   Generates water layer vectors for the balance calculation.
    *   Calculates **Initial Conditions** (SLB, SH2O) for each soil based on the requested `initial_conditions` (e.g., 20%, 50%, 100% of field capacity).
*   **Run Control Table:** Creates a master execution table (`datos_corridas`) containing every combination of:
    *   Crop
    *   Soil ID
    *   Weather Station
    *   Realization (if multiple weather realizations exist)
    *   This table controls the parallel execution loops.

### Step 3: Hybrid Series Creation (`PASO 6`)
*   **Parallel Task:** Executes `CrearSeriesHibridas` in parallel.
*   **Purpose:** Prepares the specific weather series needed for each simulation unit. It constructs a continuous daily weather record centered around the `simulation_date` and `planting_date`.
*   **Output:** Generates `series_hibridas.csv` in the output folder.

### Step 4: Phenology Estimation (`PASO 7`)
*   **Parallel Task:** Executes `EstimarFenologiaTrigo` in parallel.
*   **Logic:**
    *   Uses the weather series generated in Step 3.
    *   Loads crop genetic coefficients (`.CUL` and `.ECO` files) for the specified cultivar.
    *   Runs the `determine_phenology_stages_wheat` function to predict dates for key stages:
        *   **Sowing**
        *   **Emergence**
        *   **End Juvenile**
        *   **Anthesis** (Flowering)
        *   **End Grain Filling** (Physiological Maturity)
        *   **Harvest**
*   **Output:** Generates `fenologia.csv` containing the dates (DOY, DAP) for each stage.

### Step 5: Soil Water Balance (`PASO 8`)
*   **Parallel Task:** Executes `BalanceAguaSuelo` in parallel.
*   **Logic:**
    *   Inputs: Weather series, Phenology dates, Soil properties, and Initial conditions.
    *   **Process:**
        1.  Calculates Reference Evapotranspiration (**ETref**) using Helper functions (Hargreaves or Penman-Monteith).
        2.  Maps the calculated phenological stages to FAO-56 growth stages (Initial, Development, Mid, End).
        3.  Runs the `fao56` function (from `lib/balance_fao56.R`) to simulate the daily water balance.
        4.  Updates the soil water content daily based on precipitation (`Rain`), crop evapotranspiration (`ETc`), and percolation (`DP`).
*   **Output:** Generates `balance.agua.csv` with daily variables (Dr, TAW, Ks, ETc, etc.).

### Step 6: Finalization (`PASO 9`)
*   **Cleanup:** Stops the logger and script timer.
*   **Info File:** Copies the run log to the output directory as a `.info` file for record-keeping.

## Detailed File Descriptions

### `01_Balance_agua.R`
The main orchestrator. It sets up parallel processing tasks (`Task` object from `lib/R/Task.R`) for the three main computational steps (Series, Phenology, Balance) to ensure efficiency.

### `lib/workers.R`
Contains the high-level wrapper functions that are called by the main script. These functions:
*   `CrearSeriesHibridas`: Prepares weather data for the specific run.
*   `EstimarFenologiaTrigo`: Iterates over weather stations/years to compute phenology.
*   `BalanceAguaSuelo`: Iterates over computed phenology series to run the FAO-56 model.

### `lib/fenologia.R`
Contains the detailed biological logic for Wheat (Triticum aestivum).
*   `crown_temperatures_wheat`: Simulates crown temperature.
*   `thermal_time_calculation`: Computes Growing Degree Days (GDD).
*   `vernalization`: Simulates vernalization requirements.
*   `determine_phenology_stages_wheat`: The main function determining dates for stages like "End Juvenile", "Anthesis", etc.

### `lib/balance_fao56.R`
Implements the "dual Kc" curve methodology from FAO Irrigation and Drainage Paper No. 56.
*   `fao56`: The main loop simulating daily steps.
*   `advance`: Updates state variables (root depth, plant height, Kc values) and calculates daily stresses and evaporation.

### `lib/helpers.R`
Utility functions for soil profile handling, converting DSSAT soil layers into the format required for the balance model.

## Outputs

All outputs are saved in `outputs/{run_id}/`.

*   **`control.rds`**: An RDS file containing the master table of all simulation runs to be performed (combinations of crop, soil, year, etc.).
*   **`series_hibridas.csv`**: The weather data series used for the simulations.
*   **`fenologia.csv`**: A CSV file containing the date (YYYY-MM-DD), DOY, and DAP (Days After Planting) for each phenological stage for every simulation run.
*   **`balance.agua.csv`**: The comprehensive daily output. Key columns:
    *   `ETc`, `ETref`: Evapotranspiration.
    *   `Rain`: Precipitation.
    *   `Dr`, `TAW`: Root zone depletion and Total Available Water.
    *   `Ks`: Water stress coefficient (0-1).
    *   `Date`: Simulation date.
