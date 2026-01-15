# -----------------------------------------------------------------------------#
# --- PASO 1. Cargar paquetes necesarios ----
rm(list = ls()); gc()
Sys.setenv(TZ = "UTC")
list.of.packages <- c("ADGofTest", "caret", "dplyr", "fitdistrplus", "lmomco", "logspline",
                      "lubridate", "magrittr", "mgcv", "purrr", "SCI", "sirad", "SPEI", 
                      "stats", "stringr", "utils", "WRS2", "yaml", "yardstick", "feather",
                      "R6", "futile.logger", "mgcv", "doSNOW", "foreach", "snow", "parallel",
                      "RPostgres", "data.table", "purrr", "automap", "sf")
for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    stop(paste0("Paquete no encontrado: ", pack))
  }
}
rm(pack); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 2. Leer archivo de configuracion ----
# -----------------------------------------------------------------------------#

normalize_dirnames <- function(dirnames) {
  if (is.atomic(dirnames)) 
    dirnames <- base::sub('/$', '', dirnames)
  if (!is.atomic(dirnames))
    for (nm in names(dirnames)) 
      dirnames[[nm]] <- normalize_dirnames(dirnames[[nm]])
  return (dirnames)
}

# a) YAML de configuracion del cálculo del balance
args <- base::commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  archivo.config <- args[1]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.config <- paste0(getwd(), "/configuracion.yml")
}
if (! file.exists(archivo.config)) {
  stop(paste0("El archivo de configuración ", archivo.config, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.config, "...\n"))
  config <- yaml::yaml.load_file(archivo.config)
  config$dir <- normalize_dirnames(config$dir)
}

# b) YAML de parametros del balance
if (length(args) > 1) {
  archivo.params <- args[2]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.params <- paste0(getwd(), "/parametros.yml")
}
if (! file.exists(archivo.params)) {
  stop(paste0("El archivo de parámetros ", archivo.params, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de parámetros ", archivo.params, "...\n"))
  config$params <- yaml::yaml.load_file(archivo.params)
}

replace_run_identifier <- function(filenames, identifier) {
  if (is.atomic(filenames) && grepl("<\\*idc>", filenames)) 
    filenames <- base::sub('<\\*idc>', identifier, filenames)
  if (!is.atomic(filenames))
    for (nm in names(filenames)) 
      filenames[[nm]] <- replace_run_identifier(filenames[[nm]], identifier)
  return (filenames)
}

# c) YAML de configuraci?n de la API
if (length(args) > 1) {
  archivo.config.api <- args[3]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.config.api <- paste0(getwd(), "/configuracion_api.yml")
}
if (! file.exists(archivo.config.api)) {
  stop(paste0("El archivo de configuraci?n ", archivo.config.api, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuraci?n ", archivo.config.api, "...\n"))
  config$api <- yaml::yaml.load_file(archivo.config.api)
}

rm(archivo.config, archivo.params, archivo.config.api, args); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 3. Cargar librerias propias e iniciar script ----
# -----------------------------------------------------------------------------#

# a) Cargar librerias
#source(glue::glue("{config$dir$lib}/FechaUtils.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/R/Script.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/R/Task.R"), echo = FALSE)
#source(glue::glue("{config$dir$lib}/Helpers.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/helpers.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/workers.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/funciones_auxiliares.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/fenologia.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/funciones_evapotranspiracion.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/balance_fao56.R"), echo = FALSE)

# c.1) Definir nombre del script
script_name <- "BalanceAgua"

# c.2) Crear carpeta de logs (si no existe)
logs_folder <- glue::glue('{here::here()}/{config$dir$logs}')
if (!fs::dir_exists(logs_folder))
  fs::dir_create(logs_folder)

# c.3) Definir nombre del archivo de log 
script_logfile <- glue::glue("{logs_folder}/{script_name}.log")

# c.4) Borrar archivo .log de corridas anteriores
if (fs::file_exists(script_logfile))
  fs::file_delete(script_logfile)

# c.5) Iniciar script
script <- Script$new(run.dir = logs_folder, name = script_name, create.appender = T)
script$start()

# e) Cargar funciones necesarias
source(glue::glue("./lib/crc-api.R"), echo = FALSE)

# f) Variable para almacenar los posibles errores
errors <- c()

# Configuracion SSL
httr::set_config( httr::config(ssl_verifypeer = FALSE))

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# --- PASO 4. Cargar datos climaticos ----
# ---------------------------------------------------------------------------- #

# Cargar datos climaticos completos
datos <- readr::read_csv(glue::glue("./{config$dir$weather}/{config$params$run_id}/{config$params$weather$file}"))

# ---------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------#
# --- PASO 5. Leer/procesar variables de entrada                            ----
# -----------------------------------------------------------------------------#

# Obtener parámetros en parametros.yaml
cultivos <- config$params$crops
años <- base::seq(from = config$params$weather$start_year, to = config$params$weather$end_year)
realizaciones <- base::seq(from = 1, to = config$params$weather$realizations)
suelos <- purrr::pmap_dfr(
  .l = utils::tail(config$params$soils, -1) %>% purrr::transpose(), 
  .f = function(..., nombres_columnas) { 
    tibble::tibble(..., .name_repair = "minimal") %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$soils[[1]]
) 

# Archivo con todos los suelos disponibles! 
soil_file <- glue::glue('{config$dir$shared}/{config$soil_file}')

script$info("Calculando profundidad máxima para cada suelo")
profundidades_maximas <- suelos %>% dplyr::distinct(soil_id) %>%
  dplyr::rowwise() %>% dplyr::mutate(
    max_depth = buscar_maxima_profundidad(soil_id, soil_file),
    water_layers = list(generate_water_layers_vector(max_depth)),
    water_layers_depth = list(generate_water_layers_depth_vector(max_depth))
  ) %>% dplyr::ungroup() 
suelos %<>% dplyr::left_join(profundidades_maximas, by = 'soil_id')

script$info("Calculando condiciones inciales para cada suelo")
condiciones_iniciales <- suelos %>% dplyr::distinct(soil_id) %>%
  tidyr::crossing(initial_water = config$params$initial_conditions) %>%
  dplyr::rowwise() %>% dplyr::mutate(
    condicion_inicial = list(
      calcular_condiciones_iniciales(soil_id, soil_file, initial_water)
    ),
    SLB = list(condicion_inicial$SLB),
    SH2O = list(condicion_inicial$SH2O)
  ) %>% dplyr::ungroup() %>% dplyr::select(-condicion_inicial)

script$info("Controlando datos climáticos")
estaciones_necesarias <- suelos %>% dplyr::distinct(weather_id) %>% dplyr::pull()
estaciones_con_datos  <- unique(datos$station_id)

if (any(!estaciones_necesarias %in% estaciones_con_datos)) {
  warning("No hay datos para algunas de las estaciones. ", 
          "Las estaciones sin datos serán excluidas!!")
  suelos %<>% dplyr::filter(weather_id %in% estaciones_con_datos)
}

script$info("Definiendo tibble de control de la ejecución")
datos_corridas <- tibble::tibble(crop = cultivos) %>%
  tidyr::crossing(realization = realizaciones, suelos) %>%
  dplyr::mutate(run_id = 1:n(),
                first_year = config$params$weather$start_year, 
                last_year = config$params$weather$end_year) %>% 
  dplyr::select(run_id, dplyr::everything())


script$info("Guardando tibble de control de la ejecución")
control_file <- glue::glue('{here::here()}/{config$dir$outputs}/{config$params$run_id}/control.rds')
saveRDS(datos_corridas, file = control_file)

# Parametros del cultivo para el modelo de Balance.
# Fuente: FAO 56
parametros.balance <- purrr::pmap_dfr(
  .l = utils::tail(config$params$parametros.fao, -1) %>% purrr::transpose(), 
  .f = function(..., nombres_columnas) { 
    tibble::tibble(..., .name_repair = "minimal") %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$parametros.fao[[1]]
) 

# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# --- PASO 6. Creacion de series hibirdas ----
# ---------------------------------------------------------------------------- #

# Definir el nombre de la función a ser paralelizada
function_name <- "CrearSeriesHibridas"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{logs_folder}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{logs_folder}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (fs::file_exists(task_logfile))
  fs::file_delete(task_logfile)
if (fs::file_exists(task_outfile))
  fs::file_delete(task_outfile)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Creando Series climáticas híbridas")
# Ejecutar tarea distribuida
series.hibridas.resultados <- task$run(input.values = datos_corridas,
                           climate.data = datos,
                           number.of.processes = config$max.procesos)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.series.hibridas.errors <- task$getErrors()
if (length(task.series.hibridas.errors) > 0) {
  for (error.obj in task.series.hibridas.errors) {
    id_column <- IdentificarIdColumn(ubicaciones_a_procesar %>% dplyr::top_n(1))
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}): {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# Transformar resultados a un objeto de tipo tibble
series.hibridas.resultados.tibble <- series.hibridas.resultados %>% purrr::map_dfr(~.x)

# Guardar resultados en un archivo fácil de compartir
results_filename <- glue::glue("./{config$dir$outputs}/{config$params$run_id}/series_hibridas.csv")
script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
data.table::fwrite(series.hibridas.resultados.tibble, file = results_filename, nThread = config$max.procesos)
# ---------------------------------------------------------------------------- 

# ---------------------------------------------------------------------------- #
# --- PASO 7. Estimacion de fechas criticas (Fenologia) ----
# ---------------------------------------------------------------------------- #

# Definir el nombre de la función a ser paralelizada
function_name <- "EstimarFenologiaTrigo"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{logs_folder}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{logs_folder}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (fs::file_exists(task_logfile))
  fs::file_delete(task_logfile)
if (fs::file_exists(task_outfile))
  fs::file_delete(task_outfile)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Estimando fenología para el cultivo de trigo")
# Ejecutar tarea distribuida
fenologia.resultados <- task$run(input.values = datos_corridas,
                                 climate.data = series.hibridas.resultados.tibble,
                                 script = script,
                                 config = config,
                                 number.of.processes = config$max.procesos)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.fenologia.errors <- task$getErrors()
if (length(task.fenologia.errors) > 0) {
  for (error.obj in task.fenologia.errors) {
    id_column <- IdentificarIdColumn(ubicaciones_a_procesar %>% dplyr::top_n(1))
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}): {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# Transformar resultados a un objeto de tipo tibble
fenologia.resultados.tibble <- fenologia.resultados %>% purrr::map_dfr(~.x)

# Guardar resultados en un archivo fácil de compartir
results_filename <- glue::glue("./{config$dir$outputs}/{config$params$run_id}/fenologia.csv")
script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
data.table::fwrite(fenologia.resultados.tibble, file = results_filename, nThread = config$max.procesos)

# ---------------------------------------------------------------------------- 

# ---------------------------------------------------------------------------- #
# --- PASO 8. Estimacion de Balance de agua en el suelo ----
# ---------------------------------------------------------------------------- #

# Definir el nombre de la función a ser paralelizada
function_name <- "BalanceAguaSuelo"

# Definir nombre de archivos .log y .out de corridas anteriores
task_logfile <- glue::glue("{logs_folder}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{logs_folder}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (fs::file_exists(task_logfile))
  fs::file_delete(task_logfile)
if (fs::file_exists(task_outfile))
  fs::file_delete(task_outfile)

# Crear tarea distribuida y ejecutarla
task <- Task$new(parent.script = script,
                 func.name = function_name,
                 packages = list.of.packages)

# Informar inicio de ejecución 
script$info("Estimando balance de agua en el suelo")
# Ejecutar tarea distribuida
balance.resultados <- task$run(input.value = datos_corridas,
                                 climate.data = series.hibridas.resultados.tibble, 
                                 parametros.cultivo = parametros.balance, 
                                 fenologia.serie = fenologia.resultados.tibble, 
                                 initial.conditions = condiciones_iniciales, 
                                 soil_file = soil_file, 
                                 script = script, 
                                 config = config, 
                                 number.of.processes = config$max.procesos)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.balance.errors <- task$getErrors()
if (length(task.balance.errors) > 0) {
  for (error.obj in task.balance.errors) {
    id_column <- IdentificarIdColumn(ubicaciones_a_procesar %>% dplyr::top_n(1))
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}): {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# Transformar resultados a un objeto de tipo tibble
balance.resultados.tibble <- balance.resultados %>% purrr::map_dfr(~.x)

# Guardar resultados en un archivo fácil de compartir
results_filename <- glue::glue("./{config$dir$outputs}/{config$params$run_id}/balance.agua.csv")
script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
data.table::fwrite(balance.resultados.tibble, file = results_filename, nThread = config$max.procesos)
# ---------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------#
# --- PASO 9. Finalizar script cerrando conexion a base de datos ----
# -----------------------------------------------------------------------------#

# a) Finalizar script
script$stop()

# b) Crear archivo .info
# info_filename <- glue::glue("{config$dir$data}/{config$files$indices_sequia$info_corrida}")
# if (file.exists(info_filename))
#   file.remove(info_filename)
# file.copy(from = script_logfile, to = info_filename)

# ------------------------------------------------------------------------------
