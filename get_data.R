rm(list = ls())
# to install simplefmi:
# devtools::install_github("paasim/simplefmi")
library(tidyverse)
library(forcats)
library(stringr)
library(feather)
library(lubridate)
library(simplefmi)


# Read the cycling data
data_url <- "http://www.hel.fi/hel2/tietokeskus/data/helsinki/ksv/Helsingin_pyorailijamaarat.csv"

df1 <- read_csv2(data_url, col_names = TRUE,
                 col_types = cols("Päivämäärä" = "c", .default = "i"),
                 locale = locale(encoding = "ISO-8859-1", decimal_mark = ","))

# only include baana and exclude every observation before 2014 as the format
# before that is not consistent with the format after that & remove duplicates
df_onlybaana <- select(df1, one_of(c("Päivämäärä","Baana"))) %>%
  setNames(c("date", "count")) %>%
  filter(as.numeric(str_extract(date, "(?=\\D*\\d+\\D*)\\d{4}")) >= 2014) %>%
  distinct()


# functions to transform the date string into numeric variables
gen_map <- function(vals, targs) {
  y <- as.character(vals) %>% setNames(targs)
  function(x) str_replace_all(x, y)
}
days_en <- locale()$date_names$day_ab

loc <- locale(date_names = date_names_lang("fi"))
mon_map <- gen_map(1:12, str_sub(loc$date_names$mon, end = -6))
day_map <- gen_map(days_en, loc$date_names$day_ab)
varnames <- c("wkday", "day", "mon", "year", "hr")

# get a tibble with dates in separate variables
df_dates <- separate(df_onlybaana, "date", varnames, sep = " ") %>%
  transmute(date = make_date(year, mon_map(mon), day),
            wkday = parse_factor(day_map(wkday), levels = days_en[c(2:7,1)]),
            count = count) %>%
  group_by(date, wkday) %>%
  summarise(count = sum(count)) %>%
  ungroup()

n <- nrow(df_dates)
df_dates$day_ind <- 1:n - round(n/2) # add zero centered index for days


# get the respective weather data:
fmi_apikey <- readLines("apik") # fmi-apikey
# station id for kaisaniemi from
# http://en.ilmatieteenlaitos.fi/observation-stations
station_id <- "100971"

weather <- fmi_download(fmi_apikey, df_dates$date[1], df_dates$date[n],
                        station_id, hourly = FALSE)

# combine the data frames
df_all <- mutate(weather, rain = pmax(rain, 0)) %>% # rain = -1 => no rain.
  left_join(df_dates, "date") %>%
  select(date, wkday, day_ind, temp, rain, count) #reorder the columns

write_feather(df_all, "data.feather")
