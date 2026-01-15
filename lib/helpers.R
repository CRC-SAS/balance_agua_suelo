IdentificarIdColumn <- function(ubicacion) {
  if ("station_id" %in% colnames(ubicacion)) {
    id_column <- "station_id"
  } else if ("point_id" %in% colnames(ubicacion)) {
    id_column <- "point_id"
  } else {
    id_column <- ubicacion %>% dplyr::select(dplyr::ends_with('_id')) %>% 
      base::colnames() %>% base::sort() %>% dplyr::first()
  }
  return (id_column)
}

generate_water_layers_vector <- function(max_depth) {
  resultado <- c()
  profundidades <- c(5,10,15,15,15,30,30,30,30,30)
  for (i in 1:10) {
    limite_inferior <- sum(profundidades[1:i-1])
    limite_superior <- sum(profundidades[1:i])
    if (max_depth - limite_superior >= 0) {
      if (i<10) {
        resultado <- c(resultado, profundidades[i])
      } else {
        resultado <- c(resultado, max_depth-limite_inferior)
      }
    } else if (max_depth > limite_inferior) {
      resultado <- c(resultado, max_depth-limite_inferior)
    }
  }
  return (resultado)
}

generate_water_layers_depth_vector <- function(max_depth) {
  resultado <- c()
  profundidades <- c(5,15,30,45,60,90,120,150,180,210)
  for (i in 1:10) {
    limite_inferior <- profundidades[i-1]
    limite_superior <- profundidades[i]
    if (max_depth - limite_superior >= 0) {
      if (i<10) {
        resultado <- c(resultado, profundidades[i])
      } else {
        resultado <- c(resultado, max_depth-limite_inferior)
      }
    } else if (max_depth > limite_inferior) {
      resultado <- c(resultado, max_depth)
    }
  }
  return (resultado)
}


buscar_maxima_profundidad <- function(soil_id, soil_file) {
  
  soil_profiles <- DSSAT::read_sol(soil_file, id_soil = soil_id)
  max_depth <- base::unlist(soil_profiles$SLB) %>% dplyr::last()
  
  return (max_depth)
  
}


calcular_condiciones_iniciales <- function(soil_id, soil_file, weather_content) {
  
  soil_profiles <- DSSAT::read_sol(soil_file, id_soil = soil_id)
  
  condiciones_iniciales <- soil_profiles %>%
    dplyr::select(SLB, SLLL, SDUL) %>%
    tidyr::unnest(cols = c(SLB, SLLL, SDUL)) %>%
    dplyr::mutate(SH2O = weather_content * (SDUL - SLLL) + SLLL,
                  N = 1:n(), M = base::ceiling(n()/2),
                  SNH4 = ifelse(N <= M, 0.5, 0.3), 
                  SNO3 = 12 / 2^(base::seq(n())-1),
                  SNO3 = ifelse(SNO3 < 0.5, 0.2, SNO3),
                  SNO3 = base::trunc(SNO3*10)/10) %>%
    dplyr::select(SLB, SH2O, SNH4, SNO3)
  
  return(condiciones_iniciales)
  
}


calcular_contenido_hidrico_inicial <- function(soil_file, water_initial_contidions = c(0.2, 0.5, 1)) {
  
  soil_profiles <- DSSAT::read_sol(soil_file)
  
  if (!is.vector(water_initial_contidions)) {
    stop('Soil water initial content must be a vector')
  }
  
  combinacion_pedon_contenido_inicial <- 
    purrr::cross2(.x  = unique(soil_profiles$PEDON), .y = water_initial_contidions) %>%
    purrr::transpose()
  
  soil_initial_water <- purrr::pmap_dfr(
    .l = combinacion_pedon_contenido_inicial,
    .f = function(pedon, contenido_inicial) {
      
      suelo_individual <- soil_profiles %>%
        dplyr::filter(., PEDON == pedon) %>%
        dplyr::select(SLB, SLLL, SDUL) %>%
        tidyr::unnest(., cols = c(SLB, SLLL, SDUL)) %>%
        dplyr::mutate(SOILID = pedon,
                      IC = contenido_inicial, 
                      SH2O = contenido_inicial * (SDUL - SLLL) + SLLL) %>%
        dplyr::select(SOILID, IC, SLLL, SDUL, SH2O) %>%
        dplyr::mutate(N = 1:n(), M = base::ceiling(n()/2),
                      SNH4 = ifelse(N <= M, 0.5, 0.3), 
                      SNO3 = 12 / 2^(base::seq(n())-1),
                      SNO3 = ifelse(SNO3 < 0.5, 0.2, SNO3),
                      SNO3 = base::trunc(SNO3*10)/10) %>%
        dplyr::select(SOILID, IC, SLB, SH2O, SNH4, SNO3)
      
    }
  )
  
  return(soil_initial_water)
  
}

agua_disponible_perfil <- function(soil_file, soil_id) {
  
  # Leer datos de suelo
  suelo <- DSSAT::read_sol(soil_file, id_soil = soil_id)
  
  
  suelos_balance <- suelo  %>%
    dplyr::select(PEDON, SLB, SLMH, SLLL, SDUL, SSAT) %>% 
    tidyr::unnest_longer(col = c("PEDON", "SLB", "SLMH", "SLLL", "SDUL", "SSAT")) %>%
    dplyr::arrange(SLB) 
  
  # Maxima profundidad de suelo
  max_soil_depth <- base::unlist(suelo$SLB) %>% dplyr::last()
  
  soil_layers_depth <- suelos_balance %>%
    dplyr::pull(SLB) 
  
  soil_layers_wilting_point <- suelos_balance %>%
    dplyr::pull(SLLL) 
  
  soil_layers_field_capacity <- suelos_balance %>%
    dplyr::pull(SDUL) 
  
  soil_layers_saturation <- suelos_balance %>%
    dplyr::pull(SSAT) 
  
  soil_layers_thickness <- 
    soil_layers_depth - dplyr::lag(soil_layers_depth, default = 0)
  
  soil_layers_count <- length(soil_layers_thickness)
  
  soil_hydraulic_data <- 
    tibble::tibble(soil_id = soil_id,
                   horizon_id = 1:soil_layers_count,
                   horizon_depth = soil_layers_thickness,
                   SLB = soil_layers_depth,
                   SLLL = soil_layers_wilting_point,
                   SDUL = soil_layers_field_capacity,
                   SSAT = soil_layers_saturation) %>%
    dplyr::mutate(AW = SDUL - SLLL)
  
  return(soil_hydraulic_data)
}
