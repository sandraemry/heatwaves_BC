library(terra)
library(lubridate)
library(zoo)

# Outputs: For each year 2021–2025, this writes:
#   
# hw_days_ctx90_YEAR.tif
# daily raster stack; each layer is 1 if that day belongs to a heatwave event, 0 otherwise
# 
# threshold_ctx90_YEAR.tif
# daily threshold raster stack for that year
# 
# hw_day_count_ctx90_YEAR.tif
# number of heatwave days per pixel in that year
# 
# hw_event_count_ctx90_YEAR.tif
# number of distinct heatwave events per pixel in that year
# 
# hw_max_duration_ctx90_YEAR.tif
# longest event duration per pixel in that year
# 
# hw_mean_intensity_ctx90_YEAR.tif
# mean exceedance above threshold during heatwave days
# =========================================================
# USER SETTINGS
# =========================================================

data_dir   <- "./00_data/era5_t2m_bc"
out_dir    <- "./02_outdata/atmospheric_heatwaves_ctx90"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_years <- 2021:2025
baseline_len <- 40
window_half  <- 7   # 15-day window = doy +/- 7
p_thresh     <- 0.90
min_duration <- 3

# =========================================================
# HELPERS
# =========================================================

# Leap year check
is_leap <- function(year) {
  lubridate::leap_year(year)
}

# Remove Feb 29 so all years have 365 days for climatology
drop_feb29 <- function(dates, rast_obj) {
  keep <- !(month(dates) == 2 & mday(dates) == 29)
  list(
    dates = dates[keep],
    rast  = rast_obj[[which(keep)]]
  )
}

# Day-of-year on a 365-day calendar (after Feb 29 removed)
doy365 <- function(dates) {
  yday(dates)
}

# Circular 15-day window around doy
get_window_doys <- function(d, half_window = 7) {
  x <- ((d - half_window):(d + half_window))
  ((x - 1) %% 365) + 1
}

# Read one yearly NetCDF, convert Kelvin -> Celsius if needed
read_year_tmax <- function(year, data_dir) {
  
  files <- list.files(
    data_dir,
    pattern = paste0("^era5_t2m_", year, "_[0-9]{2}_bc\\.nc$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    stop("No files found for year: ", year)
  }
  
  files <- sort(files)
  
  daily_list <- list()
  date_list  <- list()
  
  for (i in seq_along(files)) {
    f <- files[i]
    fname <- basename(f)
    
    ym <- sub("^era5_t2m_([0-9]{4})_([0-9]{2})_bc\\.nc$", "\\1-\\2", fname)
    start_time <- as.POSIXct(paste0(ym, "-01 00:00:00"), tz = "UTC")
    
    r <- rast(f)
    
    # convert K -> C 
    vals <- global(r[[1]], "mean", na.rm = TRUE)[1, 1]
    if (!is.na(vals) && vals > 100) {
      r <- r - 273.15
    }
    
    # assign hourly timestamps
    times <- seq(
      from = start_time,
      by = "1 hour",
      length.out = nlyr(r)
    )
    
    # convert hourly -> daily Tmax
    day_index <- as.Date(times, tz = "UTC")
    r_daily <- tapp(r, index = day_index, fun = max, na.rm = TRUE)
    
    daily_list[[i]] <- r_daily
    date_list[[i]]  <- as.Date(names(table(day_index)))
  }
  
  r_year <- rast(daily_list)
  dates  <- do.call(c, date_list)
  
  if (length(dates) != nlyr(r_year)) {
    stop("Date vector length does not match number of daily layers for year ", year)
  }
  
  out <- drop_feb29(dates, r_year)
  out
}

# Run-length detector for one cell's logical vector
# Returns:
# 1 on all days that are part of a HW event, 0 otherwise
mark_heatwave_days <- function(x, min_duration = 3) {
  # x should be 0/1 or FALSE/TRUE, may include NA
  x[is.na(x)] <- 0
  r <- rle(x == 1)
  out <- inverse.rle(list(
    lengths = r$lengths,
    values  = ifelse(r$values & r$lengths >= min_duration, 1, 0)
  ))
  as.integer(out)
}

# =========================================================
# READ ALL YEARS ONCE
# =========================================================

message("Reading yearly daily Tmax rasters...")
year_list <- lapply(1981:2025, function(y) read_year_tmax(y, data_dir))
names(year_list) <- as.character(1981:2025)

# quick checks
n_layers <- sapply(year_list, function(x) nlyr(x$rast))
if (!all(n_layers == 365)) {
  stop("After dropping Feb 29, every year should have 365 layers.")
}

# combine all years into one big SpatRaster
all_rast  <- rast(lapply(year_list, `[[`, "rast"))
all_dates <- do.call(c, lapply(year_list, `[[`, "dates"))

# check
stopifnot(nlyr(all_rast) == length(all_dates))

# day of year on 365-day calendar
all_doy <- doy365(all_dates)
all_yr  <- year(all_dates)

# =========================================================
# FUNCTION TO BUILD DAILY THRESHOLDS FOR ONE TARGET YEAR
# =========================================================
# For each target year:
# baseline = previous 40 years
# threshold for each doy = 90th percentile over all baseline days
# matching doy +/- 7 days

build_thresholds_for_year <- function(target_year,
                                      all_rast,
                                      all_dates,
                                      all_doy,
                                      all_yr,
                                      baseline_len = 40,
                                      half_window = 7,
                                      p_thresh = 0.90) {
  
  baseline_years <- (target_year - baseline_len):(target_year - 1)
  
  idx_base <- which(all_yr %in% baseline_years)
  if (length(idx_base) != baseline_len * 365) {
    stop("Baseline for ", target_year, " is incomplete.")
  }
  
  base_rast <- all_rast[[idx_base]]
  base_doy  <- all_doy[idx_base]
  
  message("Building thresholds for ", target_year, "...")
  
  thr_list <- vector("list", 365)
  
  for (d in 1:365) {
    wd <- get_window_doys(d, half_window)
    idx <- which(base_doy %in% wd)
    
    # 40 years * 15 days = ~600 layers
    r_sub <- base_rast[[idx]]
    
    thr_d <- app(
      r_sub,
      fun = function(x) quantile(x, probs = p_thresh, na.rm = TRUE, names = FALSE)
    )
    
    names(thr_d) <- paste0("thr_doy_", sprintf("%03d", d))
    thr_list[[d]] <- thr_d
  }
  
  thr <- rast(thr_list)
  thr
}

# =========================================================
# FUNCTION TO DETECT HEATWAVES FOR ONE YEAR
# =========================================================

detect_heatwaves_for_year <- function(target_year,
                                      all_rast,
                                      all_dates,
                                      all_doy,
                                      all_yr,
                                      thr_rast,
                                      min_duration = 3) {
  
  idx_target <- which(all_yr == target_year)
  tmax_year  <- all_rast[[idx_target]]
  dates_year <- all_dates[idx_target]
  doy_year   <- all_doy[idx_target]
  
  # Keep only June–September
  keep <- month(dates_year) %in% 6:9
  
  tmax_year  <- tmax_year[[which(keep)]]
  dates_year <- dates_year[keep]
  doy_year   <- doy_year[keep]
  
  # align thresholds to actual days in target year
  thr_year <- thr_rast[[doy_year]]
  
  # exceedance
  exceed <- tmax_year > thr_year
  names(exceed) <- paste0("exc_", dates_year)
  
  # mark days belonging to >=3-day runs
  hw_days <- app(exceed, fun = function(x) mark_heatwave_days(x, min_duration = min_duration))
  names(hw_days) <- paste0("hw_", dates_year)
  
  # annual summaries
  hw_day_count <- sum(hw_days, na.rm = TRUE)
  names(hw_day_count) <- paste0("hw_day_count_", target_year)
  
  # number of events
  hw_event_count <- app(hw_days, fun = function(x) {
    x[is.na(x)] <- 0
    starts <- c(x[1] == 1, diff(x) == 1)
    sum(starts, na.rm = TRUE)
  })
  names(hw_event_count) <- paste0("hw_event_count_", target_year)
  
  # longest event length
  hw_max_duration <- app(hw_days, fun = function(x) {
    x[is.na(x)] <- 0
    r <- rle(x == 1)
    if (!any(r$values)) return(0)
    max(r$lengths[r$values])
  })
  names(hw_max_duration) <- paste0("hw_max_duration_", target_year)
  
  # mean exceedance intensity across heatwave days only
  intensity <- ifel(hw_days == 1, tmax_year - thr_year, NA)
  hw_mean_intensity <- mean(intensity, na.rm = TRUE)
  names(hw_mean_intensity) <- paste0("hw_mean_intensity_", target_year)
  
  # maximum exceedance intensity across heatwave days only
  hw_max_intensity <- max(intensity, na.rm = TRUE)
  names(hw_max_intensity) <- paste0("hw_max_intensity_", target_year)
  
  list(
    tmax_year          = tmax_year,
    thr_year           = thr_year,
    exceed             = exceed,
    hw_days            = hw_days,
    hw_day_count       = hw_day_count,
    hw_event_count     = hw_event_count,
    hw_max_duration    = hw_max_duration,
    hw_mean_intensity  = hw_mean_intensity,
    hw_max_intensity   = hw_max_intensity,
    dates_year         = dates_year
  )
}

# =========================================================
# RUN FOR 2021-2025
# =========================================================

annual_summary_layers <- list()

for (yy in target_years) {
  
  thr_yy <- build_thresholds_for_year(
    target_year   = yy,
    all_rast      = all_rast,
    all_dates     = all_dates,
    all_doy       = all_doy,
    all_yr        = all_yr,
    baseline_len  = baseline_len,
    half_window   = window_half,
    p_thresh      = p_thresh
  )
  
  res_yy <- detect_heatwaves_for_year(
    target_year   = yy,
    all_rast      = all_rast,
    all_dates     = all_dates,
    all_doy       = all_doy,
    all_yr        = all_yr,
    thr_rast      = thr_yy,
    min_duration  = min_duration
  )
  
  # write daily outputs if wanted
  writeRaster(
    res_yy$hw_days,
    file.path(out_dir, paste0("hw_days_ctx90_", yy, ".tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    res_yy$thr_year,
    file.path(out_dir, paste0("threshold_ctx90_", yy, ".tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    res_yy$hw_day_count,
    file.path(out_dir, paste0("hw_day_count_ctx90_", yy, ".tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    res_yy$hw_event_count,
    file.path(out_dir, paste0("hw_event_count_ctx90_", yy, ".tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    res_yy$hw_max_duration,
    file.path(out_dir, paste0("hw_max_duration_ctx90_", yy, ".tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    res_yy$hw_mean_intensity,
    file.path(out_dir, paste0("hw_mean_intensity_ctx90_", yy, ".tif")),
    overwrite = TRUE
  )
  
  
  writeRaster(
    res_yy$hw_max_intensity,
    file.path(out_dir, paste0("hw_max_intensity_ctx90_", yy, ".tif")),
    overwrite = TRUE
  )
  
  annual_summary_layers[[as.character(yy)]] <- c(
    res_yy$hw_day_count,
    res_yy$hw_event_count,
    res_yy$hw_max_duration,
    res_yy$hw_mean_intensity,
    res_yy$hw_max_intensity
  )
}

# stack annual summaries together
annual_summary <- rast(annual_summary_layers)
writeRaster(
  annual_summary,
  file.path(out_dir, "hw_summary_ctx90_stack.tif"),
  overwrite = TRUE
)

############### summary rasters #################
out_dir_sum <- "./02_outdata/atmospheric_heatwaves_ctx90/summary"
dir.create(out_dir_sum, recursive = TRUE, showWarnings = FALSE)

#stack the event days files 
daycount_files <- list.files("./02_outdata/atmospheric_heatwaves_ctx90/",
                             pattern = "hw_day_count_.*\\.tif$",
                             full.names = TRUE)

daycount_sum <- sum(rast(daycount_files), na.rm = TRUE)

writeRaster(
  daycount_sum,
  file.path(out_dir_sum, "hw_summary_ctx90_daycount_sum.tif"),
  overwrite = TRUE
)

#stack the event count files 

count_files <- list.files("./02_outdata/atmospheric_heatwaves_ctx90/",
                          pattern = "hw_event_count_.*\\.tif$",
                          full.names = TRUE)

eventcount_sum <- sum(rast(count_files), na.rm = TRUE)

writeRaster(
  eventcount_sum,
  file.path(out_dir_sum, "hw_summary_ctx90_eventcount_sum.tif"),
  overwrite = TRUE
)

# stack intensity files
int_files <- list.files(
  "./02_outdata/atmospheric_heatwaves_ctx90/",
  pattern = "hw_mean_intensity_.*\\.tif$",
  full.names = TRUE
)

intensity_mean <- mean(rast(int_files), na.rm = TRUE)

writeRaster(
  intensity_mean,
  file.path(out_dir_sum, "hw_summary_ctx90_intensity_mean.tif"),
  overwrite = TRUE
)


# stack max intensity files
maxint_files <- list.files(
  "./02_outdata/atmospheric_heatwaves_ctx90/",
  pattern = "hw_max_intensity_.*\\.tif$",
  full.names = TRUE
)

intensity_max <- mean(rast(maxint_files), na.rm = TRUE)

writeRaster(
  intensity_max,
  file.path(out_dir_sum, "hw_summary_ctx90_intensity_max.tif"),
  overwrite = TRUE
)



# Checking a few specific coordinates/dates -------------------------------

library(tidyverse)

# ----------------------------
# User settings
# ----------------------------
data_dir <- "./00_data/era5_t2m_bc"
lat <- 50.230864
lon <- -121.582451

start_date <- as.POSIXct("2021-06-01 00:00:00", tz = "UTC")
end_date   <- as.POSIXct("2021-07-31 23:00:00", tz = "UTC")

# ----------------------------
# Find the June and July 2021 files
# ----------------------------
files <- list.files(
  data_dir,
  pattern = "^era5_t2m_2021_(06|07)_bc\\.nc$",
  full.names = TRUE
)

files <- sort(files)

if (length(files) == 0) {
  stop("No June/July 2021 ERA5 files found.")
}

# ----------------------------
# Read files, extract hourly temperature at point
# ----------------------------
ts_list <- vector("list", length(files))

for (i in seq_along(files)) {
  
  f <- files[i]
  fname <- basename(f)
  
  ym <- sub("^era5_t2m_([0-9]{4})_([0-9]{2})_bc\\.nc$", "\\1-\\2", fname)
  start_time_file <- as.POSIXct(paste0(ym, "-01 00:00:00"), tz = "UTC")
  
  r <- rast(f)
  
  # Convert Kelvin to Celsius if needed
  vals <- global(r[[1]], "mean", na.rm = TRUE)[1, 1]
  if (!is.na(vals) && vals > 100) {
    r <- r - 273.15
  }
  
  # Build hourly timestamps
  times <- seq(
    from = start_time_file,
    by = "1 hour",
    length.out = nlyr(r)
  )
  
  # Extract point values
  pt <- vect(data.frame(x = lon, y = lat), geom = c("x", "y"), crs = "EPSG:4326")
  ext <- terra::extract(r, pt)
  
  ts_list[[i]] <- tibble(
    datetime_utc = times,
    temperature_c = as.numeric(ext[1, -1])
  )
}

temp_ts <- bind_rows(ts_list) %>%
  filter(datetime_utc >= start_date,
         datetime_utc <= end_date)

# Optional: convert to local time for plotting in BC summer time
temp_ts <- temp_ts %>%
  mutate(datetime_vancouver = with_tz(datetime_utc, tzone = "America/Vancouver"))

# ----------------------------
# Plot
# ----------------------------
ggplot(temp_ts, aes(x = datetime_vancouver, y = temperature_c)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Air temperature (°C)",
    title = "ERA5 air temperature time series",
    subtitle = "49.272726, -123.192864 | June 1 to July 31, 2021"
  ) +
  theme_minimal()
