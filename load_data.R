import_data <- function() {
  data <- readr::read_csv("data.csv", col_types = c("i", "d", "d"))
  data %>% 
    dplyr::filter(dplyr::between(age, 23000, 93000)) %>% 
    dplyr::group_by(age) %>% 
    dplyr::summarise(
      d180 = mean(d18O),
      Ca2 = mean(Ca2)
    ) %>% 
    dplyr::mutate(
      age = age / 1000,
      logCa2 = -log(Ca2) - mean(-log(Ca2)),
      Y = 2 * sqrt(Ca2)
    )
}


