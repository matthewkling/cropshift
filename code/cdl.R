library(tidyverse)
library(sf)
library(rgdal)
library(furrr)

plan(multisession, workers = 8)



# load Earth Engine dataset generated at 
# https://code.earthengine.google.com/824794f20ce4320ac8625b49c9a49fec
d <- read_csv("data/county_crop_freq_ts_raw.csv") %>%
   select(counts, county_fips, state_fips, year)

# load crop ID dictionary (from same GEE script as above)
id <- xlsx::read.xlsx("data/CDL_codes_names_colors.xls", 1, startRow = 4) %>%
   select(id = Codes, crop = Class_Names)
# id <- read_csv("data/cdl_bands.csv") %>%
#    rename(id = value,
#           crop = description)

# function to convert json dictionary string to data frame
unpack_dict <- function(x){
   x %>% 
      str_remove("\\{") %>% str_remove("\\}") %>%
      str_split(", ") %>%
      "[["(1) %>%
      str_split("=", simplify = T) %>%
      as.data.frame() %>%
      setNames(c("id", "freq"))
}

# convert dataset to tidy table
f <- 1:nrow(d) %>%
   map_df(function(i){
      if(d$counts[i] == "{}") return(NULL)
      unpack_dict(d$counts[i]) %>% 
         mutate(fips = paste0(d$state_fips[i], d$county_fips[i]),
                year = d$year[i])}) %>%
   as_tibble() %>%
   mutate(id = as.integer(id),
          freq = as.numeric(freq))


## rebalance frequencies using confusion matrices #####

states <- tigris::fips_codes %>%
   select(state, state_code) %>%
   distinct() %>%
   filter(! state %in% c("AK", "HI")) %>%
   slice(1:49)

years <- 2008:2021

# load confusion matrices
matrices <- years %>%
   map(function(year){
      dirs <- list.dirs("big_data")
      dir <- dirs[str_detect(dirs, as.character(year))]
      
      states$state %>%
         future_map(function(state){
            if(state == "DC") state <- "MD"
            file <- list.files(dir, full.names = T,
                               pattern = paste0("_", state))
            cm <- readxl::read_xlsx(file[1], "All - Matrix") %>% 
               select("1":"256") %>% 
               slice(1:256) %>%
               as.matrix() %>%
               apply(1, function(x) x/sum(x)) %>% 
               t()
            cm[!is.finite(cm)] <- 0
            return(cm)
         }) %>%
         setNames(states$state)
   }) %>%
   setNames(years)

# adjust frequencies using confusion matrices
f <- f %>%
   full_join(expand_grid(fips = unique(.$fips), 
                         year = unique(.$year),
                         id = 1:256)) %>%
   arrange(fips, year, id) %>%
   mutate(freq = ifelse(is.na(freq), 0, freq)) %>%
   split(paste(.$fips, .$year)) %>%
   future_map(function(x){
      year <- x$year[1]
      fips <- x$fips[1]
      state <- states$state[states$state_code == str_sub(fips, 1, 2)]
      cm <- matrices[[as.character(year)]][[state]]
      x %>% mutate(cfreq = as.vector(matrix(x$freq, 1) %*% cm))
   }) %>%
   bind_rows()

f %>%
   sample_n(100000) %>%
   ggplot(aes(freq, cfreq)) +
   geom_abline(color = "gray") +
   geom_point(size = .25) +
   scale_x_log10() +
   scale_y_log10()


# add crop names
f <- f %>%
   left_join(id) %>%
   mutate(crop = tolower(crop))









nonag <- c(0, 63:65, 81:195)

r <- f %>%
   filter(! id %in% nonag) %>%
   group_by(crop) %>% summarize(n = sum(freq)) %>% arrange(desc(n))


counties <- st_read("data/Counties", "cb_2018_us_county_500k") %>%
   mutate(area = st_area(.))

cd <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(filter(f, !is.na(crop), crop %in% r$crop[1:25], year == "2019")) %>%
   group_by(fips) %>%
   mutate(p = freq / sum(freq, na.rm = T)) %>% # normalize as prop cropland
   group_by(crop) %>%
   mutate(value = p / max(p, na.rm = T)) %>% 
   filter(!is.na(crop))

p <- ggplot(cd, aes(fill = value)) +
   facet_wrap(~crop, ncol = 5) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_gradientn(colors = c("darkolivegreen2", "forestgreen", "black"), 
                        na.value = "white") +
   theme_void() +
   theme(legend.position = "none")
# ggsave("figures/cdl_crop_maps.png", p, width = 10, height = 6, units = "in")







# marginal climate means for crops and counties ##############

centroids <- st_read("data/Counties", "cb_2018_us_county_500k") %>%
   st_centroid() %>%
   mutate(fips = paste0(STATEFP, COUNTYFP))
centroids <- bind_cols(fips = centroids$fips, 
                       as.data.frame(st_coordinates(centroids))) %>%
   rename(lon = X, lat = Y) %>%
   as_tibble()

climate <- read_csv("big_data/climate.csv") %>%
   filter(season != "growing_season",
          start == 1950) %>%
   select(fips, 
          ppt = prec, gdd, aet, cwd, edd) %>%
   mutate(fips = str_pad(fips, 5, "left", "0"))



g <- f %>%
   filter(! id %in% nonag) %>%
   left_join(climate) %>%
   left_join(centroids) %>%
   mutate(acres = freq / (10/3)^2 * 2.47105) %>%
   filter(is.finite(ppt),
          is.finite(acres)) %>%
   gather(var, county_value, ppt:lat)



p <- g %>% 
   ungroup() %>%
   filter(fips %in% sample(fips, 20)) %>%
   spread(var, county_value) %>%
   ggplot(aes(year, freq, color = crop)) +
   facet_wrap(~fips, scales = "free_y") +
   geom_line() +
   theme_bw() +
   theme(legend.position = "none",
         axis.text.y = element_blank()) +
   scale_x_continuous(breaks = range(g$year)) +
   labs(y = "pixels")
ggsave("figures/cdl/county_crop_ts.png", 
       p, width = 10, height = 6, units = "in")



# summarization

weighted_quantile <- function(x, w, q){
   # x = rnorm(100)
   # w = 0:99
   # q = .333
   bind_cols(x = x, w = w) %>%
      na.omit() %>%
      arrange(x) %>%
      mutate(cw = cumsum(w),
             cq = cw / sum(w),
             diff = q - cq,
             absdiff = abs(diff)) %>%
      arrange(absdiff) %>%
      slice(1:2) %>%
      summarize(y = weighted.mean(x, w / absdiff)) %>%
      pull(y)
}

county <- g %>%
   group_by(crop, var) %>%
   mutate(mean = weighted.mean(county_value, acres, na.rm = T),
          q10 = weighted_quantile(county_value, acres, 0.10),
          q90 = weighted_quantile(county_value, acres, 0.90)) %>%
   gather(crop_stat, crop_value, mean, q10, q90) %>%
   group_by(fips, year, var, crop_stat) %>%
   summarize(crop_value = weighted.mean(crop_value, acres),
             acres = sum(acres)) %>%
   ungroup()

crop <- g %>% 
   group_by(crop, year, var) %>%
   summarize(mean = weighted.mean(county_value, acres, na.rm = T),
             q10 = weighted_quantile(county_value, acres, 0.10),
             q90 = weighted_quantile(county_value, acres, 0.90),
             acres = sum(acres)) %>%
   gather(crop_stat, county_value, mean, q10, q90) %>%
   ungroup()



   

# temporal trends ############

# ggplot(county, aes(year, crop_cwd, color = fips)) +
#    geom_line() +
#    theme_minimal() +
#    theme(legend.position = "none")
# 
# ggplot(crop, aes(year, county_cwd, color = crop)) +
#    geom_line() +
#    theme_minimal() +
#    theme(legend.position = "none")
# 
# lm(crop_gdd ~ year, data = county, weights = acres)


# ggplot(crop, aes(year, county_cwd, group = crop, alpha = log(acres))) +
#    geom_line() +
#    theme_minimal() +
#    theme(legend.position = "none")



# marginal trends ##

d1 <- county %>% 
   group_by(year, var, crop_stat) %>% 
   summarize(value = weighted.mean(crop_value, acres)) %>%
   mutate(margin = "county")

d2 <- crop %>% 
   group_by(year, var, crop_stat) %>% 
   summarize(value = weighted.mean(county_value, acres)) %>%
   mutate(margin = "crop")

p <- bind_rows(d1, d2) %>%
   mutate(crop_stat = factor(crop_stat, levels = c("q90", "mean", "q10"))) %>%
   filter(! var %in% c("lon", "lat"),
          var %in% c("aet", "cwd")) %>%
   ggplot(aes(year, value, color = margin)) +
   facet_wrap(~ crop_stat + margin + var, scales = "free", nrow = 3) +
   geom_line() +
   geom_smooth(method = lm) +
   geom_point() +
   scale_x_continuous(breaks = range(crop$year)) +
   theme_bw() +
   theme(legend.position = "none") +
   labs(y = NULL, x = NULL)
ggsave("figures/cdl/global_marginal_trends.png", 
       p, width = 12, height = 10, units = "in")




# local LMs ##

county_lm <- county %>% 
   rename(value = crop_value) %>%
   group_by(var, crop_stat) %>%
   mutate(value = (value - mean(value)) / sd(value),
          year = year - mean(year)) %>%
   group_by(fips, var, crop_stat) %>%
   do(broom::tidy(lm(value ~ year, ., weights = acres))) %>%
   select(fips, var, crop_stat, term, estimate) %>%
   spread(term, estimate) %>%
   janitor::clean_names()

crop_lm <- crop %>% 
   rename(value = county_value) %>%
   group_by(var, crop_stat) %>%
   mutate(value = (value - mean(value)) / sd(value),
          year = year - mean(year)) %>%
   group_by(crop, var, crop_stat) %>%
   do(broom::tidy(lm(value ~ year, ., weights = acres))) %>%
   select(crop, var, crop_stat, term, estimate) %>%
   spread(term, estimate) %>%
   janitor::clean_names()


county_lm <- county %>% 
   rename(value = crop_value) %>%
   group_by(var, crop_stat) %>%
   mutate(value = (value - mean(value)) / sd(value),
          year = year - mean(year)) %>%
   group_by(fips, var, crop_stat) %>%
   do(broom::tidy(lm(value ~ year, ., weights = acres))) %>%
   select(fips, var, crop_stat, term, estimate) %>%
   spread(term, estimate) %>%
   janitor::clean_names()

corn_lm <- g %>% 
   filter(crop == "corn",
          var == var[1]) %>%
   mutate(year = year - mean(year)) %>%
   group_by(fips) %>%
   do(broom::tidy(lm(acres ~ year, .))) %>%
   select(fips, term, estimate) %>%
   spread(term, estimate) %>%
   janitor::clean_names()


# trend maps

p <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(county_lm) %>%
   filter(!is.na(var)) %>%
   mutate(year = ifelse(abs(year) > .1, sign(year) * .1, year)) %>%
   ggplot(aes(fill = year)) +
   facet_wrap(~ var + crop_stat) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        na.value = "white",
                        limits = c(-.1, .1)) +
   theme_void() +
   theme(legend.position = c(.75, .1),
         legend.direction = "horizontal")
ggsave("figures/cdl/county_trend_maps.png", 
       p, width = 10, height = 6, units = "in")

p <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(county_lm) %>%
   filter(!is.na(var),
          var %in% c("cwd", "aet")) %>%
   mutate(year = ifelse(abs(year) > .1, sign(year) * .1, year)) %>%
   ggplot(aes(fill = year)) +
   facet_wrap(~ var + crop_stat) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        na.value = "white",
                        limits = c(-.1, .1)) +
   theme_void() +
   theme(legend.position = "bottom",
         legend.direction = "horizontal")
ggsave("figures/cdl/county_trend_maps_select.png", 
       p, width = 10, height = 6, units = "in")


p <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(corn_lm) %>%
   # mutate(year = ifelse(abs(year) > .1, sign(year) * .1, year)) %>%
   ggplot(aes(fill = year / intercept)) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        na.value = "white",
                        limits = c(-.1, .1)) +
   theme_void() +
   theme(legend.position = "bottom",
         legend.direction = "horizontal") +
   labs(fill = "corn acreage trend\n(slope / intercept)")
ggsave("figures/cdl/county_corn_trend_maps.png", 
       p, width = 10, height = 6, units = "in")

# county paths
set.seed(5)
ss <- sample(county$fips, 20)
pd <- county %>%
   filter(fips %in% ss,
          crop_stat == "mean") %>%
   spread(var, crop_value)
p <- pd %>%
   ggplot(aes(aet, cwd, size = year,
              color = fips, group = fips)) +
   geom_path() +
   geom_point(data = pd %>% filter(year == max(year))) +
   scale_size_continuous(range = c(.5, 1)) +
   theme_minimal() +
   theme(legend.position = "none")
ggsave("figures/cdl/county_paths.png", 
       p, width = 6, height = 6, units = "in")



# climate trends vs crop trends

trends <- read_csv("e:/ca2cc/ca2cc-county/data/ignore/climate.csv") %>%
   filter(season != "growing_season",
          start %in% c(1950, 1982, 2007, 2012)) %>%
   select(fips, year = start, aet, cwd) %>%
   mutate(fips = str_pad(fips, 5, "left", "0")) %>%
   gather(var, value, aet, cwd) %>%
   mutate(year = paste0("y", year)) %>%
   spread(year, value) %>%
   mutate(delta_1982 = y1982 - y1950,
          delta_2007 = (y2007+y2012)/2 - y1950)

county_size <- county %>% group_by(fips) %>% summarize(acres = mean(acres))

td <- trends %>%
   right_join(county_lm) %>%
   right_join(county_size) %>%
   filter(var %in% c("cwd", "aet")) %>%
   mutate(crop_stat = factor(crop_stat, levels = c("q10", "mean", "q90")))

tsd <- td %>%
   group_by(var, crop_stat) %>%
   summarize(delta_1982 = weighted.mean(delta_1982, acres, na.rm = T),
             delta_2007 = weighted.mean(delta_2007, acres, na.rm = T),
             year = weighted.mean(year, acres, na.rm = T))

p <- td %>%
   ggplot(aes(delta_2007, year, weight = acres)) +
   facet_grid(var ~ crop_stat) +
   geom_point(color = "gray40", size = .5) +
   geom_hline(yintercept = 0, color = "black") +
   geom_vline(xintercept = 0, color = "black") +
   geom_vline(data = tsd, aes(xintercept = delta_2007), color = "red") +
   geom_hline(data = tsd, aes(yintercept = year), color = "red") +
   geom_smooth(method = "lm") +
   coord_cartesian(ylim = c(-.05, .05)) +
   theme_bw() +
   labs(y = "slope", x = "climate trend")
ggsave("figures/cdl/county_trend_vs_trend.png", 
       p, width = 8, height = 6, units = "in")


# spatial homogenization
p <- td %>%
   ggplot(aes(intercept, year, weight = acres)) +
   facet_grid(var ~ crop_stat) +
   geom_point(size = .5) +
   geom_hline(yintercept = 0, color = "gray") +
   geom_vline(xintercept = 0, color = "gray") +
   geom_smooth(method = "lm") +
   coord_cartesian(ylim = c(-.1, .1)) +
   theme_bw() +
   labs(y = "slope")
ggsave("figures/cdl/county_trend_vs_intercept.png", 
       p, width = 8, height = 6, units = "in")


# climate trend maps
pd <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(trends) %>%
   filter(var %in% c("cwd", "aet")) 
p <- pd %>%
   ggplot(aes(fill = delta_2007)) +
   facet_wrap(~var, ncol = 1) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        na.value = "white",
                        limits = max(abs(pd$delta_2007)) * c(-1, 1)) +
   theme_void() +
   theme(legend.position = "right") +
   labs(fill = "delta (mm)\n")
ggsave("figures/cdl/county_climate_trend_maps.png", 
       p, width = 5, height = 6, units = "in")


# climate vs crop trend maps
p <- pd %>%
   left_join(td) %>%
   filter(!is.na(crop_stat)) %>%
   ggplot(aes(fill = paste(sign(delta_2007), sign(year)))) +
   facet_wrap(var ~ crop_stat) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_manual(values = c("darkred", "tomato", "dodgerblue", "darkblue")) +
   theme_void() +
   theme(legend.position = "bottom") +
   labs(fill = "signs of climate and crop trends")
ggsave("figures/cdl/county_trend_comparison_maps.png", 
       p, width = 10, height = 5, units = "in")



# climate maps
p <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(climate) %>%
   gather(var, value, ppt:edd) %>%
   filter(!is.na(value),
          var %in% c("cwd", "aet")) %>%
   group_by(var) %>%
   mutate(value = scales::rescale(value)) %>%
   ggplot(aes(fill = value)) +
   facet_wrap(~var, ncol = 1) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_viridis_c() +
   theme_void() +
   theme(legend.position = "none")
ggsave("figures/cdl/county_climate_maps.png", 
       p, width = 5, height = 6, units = "in")


# county CROP climate maps
p <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(county %>% 
                filter(var %in% c("cwd", "aet"),
                       year == year[1])) %>%
   filter(is.finite(crop_value),
          !is.na(var)) %>%
   group_by(var, crop_stat) %>%
   mutate(value = scales::rescale(crop_value)) %>%
   ggplot(aes(fill = value)) +
   facet_wrap(~var, ncol = 1) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_viridis_c() +
   theme_void() +
   theme(legend.position = "none")
ggsave("figures/cdl/county_crop_climate_maps.png", 
       p, width = 5, height = 6, units = "in")



# difference between crop and actual climates 
xd <- counties %>%
   mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
   left_join(county %>% 
                filter(var %in% c("cwd", "aet"),
                       year == year[1]) %>%
                select(-year)) %>%
   left_join(climate %>%
                gather(var, climate_value, -fips) %>%
                filter(var %in% c("cwd", "aet"))) %>%
   filter(is.finite(crop_value),
          !is.na(var)) %>%
   left_join(county_lm) %>%
   mutate(diff = crop_value - climate_value) %>%
   mutate(crop_stat = factor(crop_stat, levels = c("q10", "mean", "q90")))
   
p <- xd %>%
   ggplot(aes(diff, year, weight = acres)) +
   facet_grid(var ~ crop_stat, scales = "free") +
   geom_point(size = .5) +
   geom_hline(color = "gray", yintercept = 0) +
   geom_vline(color = "gray", xintercept = 0) +
   geom_smooth(method = lm, color = "red") +
   theme_bw() +
   coord_cartesian(ylim = c(-.1, .1)) +
   labs(x = "CROP climate minus ACTUAL climate",
        y = "trend in crop climate")
ggsave("figures/cdl/county_tension_trend_scatters.png", 
       p, width = 8, height = 6, units = "in")

p <- xd %>%
   ggplot(aes(fill = diff)) +
   facet_grid(var ~ crop_stat) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        values = c(0, .4, .5, .6, 1),
                        na.value = "white",
                        limits = max(abs(xd$diff)) * c(-1, 1)) +
   theme_void() +
   theme(legend.position = "bottom") +
   labs(fill = "CROP climate minus ACTUAL climate")
ggsave("figures/cdl/county_tension_trend_maps.png", 
       p, width = 8, height = 6, units = "in")


stop("bananas")


# adaptation rate: change in crop climate, divided by disequilibrium
xd <- xd %>%
   mutate(ratio = year / diff,
          ratio = abs(ratio) ^ (1/3) * sign(ratio),
          ratio = ifelse(abs(ratio) > .1, sign(ratio) * .1, ratio))
p <- xd %>%
   ggplot(aes(diff, year, color = ratio,  weight = acres)) +
   facet_grid(var ~ crop_stat, scales = "free") +
   geom_point(size = .5) +
   geom_hline(color = "gray", yintercept = 0) +
   geom_vline(color = "gray", xintercept = 0) +
   geom_smooth(method = lm, color = "black") +
   scale_color_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        values = c(0, .4, .5, .6, 1),
                        na.value = "white",
                        limits = max(abs(xd$ratio)) * c(-1, 1)) +
   coord_cartesian(ylim = c(-.1, .1)) +
   theme_bw() +
   theme(legend.position = "top") +
   labs(x = "CROP climate minus ACTUAL climate",
        y = "trend in crop climate",
        color = "crop trend rate / disequilibrium     ")
ggsave("figures/cdl/county_tension_trend_scatters_rate.png", 
       p, width = 8, height = 6, units = "in")

p <- xd %>%
   ggplot(aes(fill = ratio)) +
   facet_grid(var ~ crop_stat) +
   geom_sf(data = counties, fill = "gray90", color = NA) +
   geom_sf(color = NA) +
   scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
   scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
   scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        values = c(0, .4, .5, .6, 1),
                        na.value = "white",
                        limits = max(abs(xd$ratio)) * c(-1, 1)) +
   theme_void() +
   theme(legend.position = "bottom") +
   labs(fill = "crop trend rate / disequilibrium     ")
ggsave("figures/cdl/county_tension_trend_maps_rate.png", 
       p, width = 8, height = 6, units = "in")



# can climate trends and tension jointly explain crop trends?
# we expect interactions

txd <- xd %>%
   left_join(trends) %>%
   group_by(var, crop_stat) %>%
   mutate(year = scale(year),
          diff = scale(diff, center = F),
          delta = scale(delta_2007, center = F))

xx <- txd %>%
   split(paste(.$var, .$crop_stat)) %>%
   map_df(function(x){
      fit <- gam(year ~ te(diff, delta), data = x, weights = acres)
      x$year_pred_gam <- predict(fit, x)
      fit <- lm(year ~ diff * delta, x, weights = acres)
      x$year_pred_lm <- predict(fit, x)
      return(x)
   })

p <- ggplot(xx, aes(diff, delta, color = year_pred_gam)) +
   facet_grid(var ~ crop_stat) +
   geom_point() +
   geom_vline(xintercept = 0) +
   geom_hline(yintercept = 0) +
   scale_color_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                        values = c(0, .4, .5, .6, 1),
                        na.value = "white",
                        limits = max(abs(xx$year_pred_gam)) * c(-1, 1)) +
   theme_bw() +
   theme(legend.position = "top") +
   labs(x = "CROP climate minus COUNTY climate",
        y = "COUNTY climate trend",
        color = "fitted trend in CROP climate  "  )
ggsave("figures/cdl/diff_delta_trend_gam_scatters.png", 
       p, width = 8, height = 6, units = "in")


## overall, averaging across stats and vars ##
x <- txd
fit <- gam(year ~ te(diff, delta), data = x, weights = acres)
x$year_pred_gam <- predict(fit, x)
fit <- lm(year ~ diff * delta, x, weights = acres)
x$year_pred_lm <- predict(fit, x)


p <- ggplot(x, aes(diff, delta, color = year_pred_gam)) +
   geom_point() +
   geom_vline(xintercept = 0) +
   geom_hline(yintercept = 0) +
   scale_color_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                         values = c(0, .45, .5, .55, 1),
                         na.value = "white",
                         limits = max(abs(x$year_pred_gam)) * c(-1, 1)) +
   coord_cartesian(xlim = max(abs(x$diff)) * c(-1, 1),
                   ylim = max(abs(x$delta)) * c(-1, 1)) +
   theme_bw() +
   theme(legend.position = "top") +
   labs(x = "CROP climate minus COUNTY climate",
        y = "COUNTY climate trend",
        color = "fitted trend in CROP climate  "  )
ggsave("figures/cdl/diff_delta_trend_gam_scatters_pooled.png", 
       p, width = 4.5, height = 5, units = "in")


xxx <- x %>%
   mutate(diffs = sign(diff),
          deltas = sign(delta)) %>%
   group_by(diffs, deltas) %>%
   mutate(year_pred_gam = (mean(year_pred_gam) + median(year_pred_gam)) / 2)
p <- ggplot(xxx, 
            aes(diff, delta, color = year_pred_gam)) +
   geom_point() +
   geom_vline(xintercept = 0) +
   geom_hline(yintercept = 0) +
   scale_color_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                         # values = c(0, .4, .5, .6, 1),
                         na.value = "white",
                         limits = max(abs(xxx$year_pred_gam)) * c(-1, 1)) +
   coord_cartesian(xlim = max(abs(x$diff)) * c(-1, 1),
                   ylim = max(abs(x$delta)) * c(-1, 1)) +
   theme_bw() +
   theme(legend.position = "top") +
   labs(x = "CROP climate minus COUNTY climate",
        y = "COUNTY climate trend",
        color = "fitted trend in CROP climate  "  )
ggsave("figures/cdl/diff_delta_trend_gam_scatters_pooled_sign.png", 
       p, width = 4.5, height = 5, units = "in")


p <- ggplot(x, aes(diff, delta, color = year_pred_gam)) +
   geom_vline(xintercept = 0) +
   geom_hline(yintercept = 0) +
   scale_color_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred"),
                         values = c(0, .45, .5, .55, 1),
                         na.value = "white",
                         limits = max(abs(x$year_pred_gam)) * c(-1, 1)) +
   coord_cartesian(xlim = max(abs(x$diff)) * c(-1, 1),
                   ylim = max(abs(x$delta)) * c(-1, 1)) +
   theme_bw() +
   theme(legend.position = "top") +
   labs(x = "CROP climate minus COUNTY climate",
        y = "COUNTY climate trend",
        color = "fitted trend in CROP climate  "  )
ggsave("figures/cdl/diff_delta_trend_gam_scatters_pooled_empty.png", 
       p, width = 4.5, height = 5, units = "in")



#########################################################################

# compare CDL to AgCensus #####################

# load and format ag census data
ac <- read_csv("data/ag_census_processed.csv") %>%
   filter(year %in% c(2012, 2017)) %>%
   select(fips, crop = commodity, year, acres = acres_all) %>%
   mutate(u = is.na(acres)) %>%
   full_join(expand(., fips, crop, year)) %>%
   arrange(fips, crop, year) %>%
   mutate(acres = ifelse(is.na(acres), 0, acres),
          crop = tolower(crop))

cd <- f %>% 
   filter(year %in% c(2012, 2017)) %>%
   full_join(expand(., fips, crop, year)) %>%
   mutate(freq = ifelse(is.na(freq), 0, freq),
          crop = tolower(crop),
          crop = ifelse(crop == "alfalfa", "hay", crop),
          crop = ifelse(crop == "chick peas", "chickpeas", crop))

# # export data for manual creation of crosswalk in excel
# ac %>% group_by(crop) %>% summarize(acres = sum(acres)) %>%
#    write_csv("data/agcensus_crop_types.csv")
# cd %>% group_by(crop) %>% summarize(freq = sum(freq)) %>%
#    write_csv("data/cdl_crop_types.csv")


compare <- function(x = "sorghum"){
   xa <- ac %>% filter(str_detect(crop, x)) %>%
      group_by(fips, year) %>%
      summarize(acres = sum(acres),
                u = any(na.omit(u))) %>%
      mutate(crop = x)
   xc <- cd %>% filter(str_detect(crop, x)) %>%
      group_by(fips, year) %>%
      summarize(freq = sum(freq)) %>%
      mutate(crop = x,
             cdl_acres = freq / (10/3)^2 * 2.47105) # 30m to hectare to acres
   xd <- full_join(xa, xc) %>%
      select(year, fips, crop, agcensus_acres = acres, cdl_acres, u)
   return(xd)
   # cor(xd$ac_acres, xd$freq, use = "pairwise.complete.obs")
   # cor(log(xd$ac_acres), log(xd$freq), use = "pairwise.complete.obs")
}


crops <- c("barley", "buckwheat", "canola", "chickpeas", "corn", "cotton",
           "flaxseed", "hay", "hops", "lentils", "mustard", "oats", "peanuts", 
           "rye", "safflower", "sorghum", "soybeans", "wheat")

comp <- map_df(crops, compare)

comp_r2 <- comp %>%
   group_by(crop) %>%
   summarize(r2 = cor(agcensus_acres, cdl_acres, use = "pairwise.complete.obs")^2,
             r2 = round(r2, 2),
             agcensus_acres = min(agcensus_acres, na.rm = T),
             cdl_acres = max(cdl_acres, na.rm = T))

p <- comp %>%
   filter(!u) %>%
   ggplot(aes(agcensus_acres, cdl_acres)) + 
   facet_wrap(~crop, nrow = 3, scales = "free") +
   geom_abline(slope = 1, intercept = 0, color = "red") +
   geom_point(size = .5, alpha = .6, shape = 16) +
   geom_text(data = comp_r2, 
             aes(label = r2),
             hjust = 0, vjust = 1, color = "blue") +
   theme_bw() +
   theme(axis.text = element_blank()) +
   labs(x = "AgCensus acreage",
        y = "CDL acreage")
ggsave("figures/cdl_agcensus_scatters.png", p, 
       width = 10, height = 5, units = "in")
ggsave("figures/cdl_agcensus_scatters_loglog.png", 
       p + scale_x_log10() + scale_y_log10() +
          labs(x = "AgCensus acreage (log scale)",
               y = "CDL acreage (log scale)"), 
       width = 10, height = 5, units = "in")

trends <- comp %>%
   filter(!u) %>%
   rename(agcensus = agcensus_acres,
          cdl = cdl_acres) %>%
   gather(dataset, value, cdl, agcensus) %>%
   mutate(year = paste0("y", year)) %>% 
   spread(year, value) %>%
   mutate(diff = y2017 - y2012,
          ratio = y2017 / y2012) %>%
   select(-y2012, -y2017) %>%
   gather(stat, value, diff, ratio) %>%
   spread(dataset, value)

trend_r2 <- trends %>%
   group_by(stat, crop) %>%
   summarize(r2 = cor(agcensus, cdl, use = "pairwise.complete.obs")^2,
             r2 = round(r2, 2),
             agcensus = min(agcensus, na.rm = T),
             cdl = max(cdl, na.rm = T))

p <- trends %>% 
   filter(stat == "diff") %>%
   ggplot(aes(agcensus, cdl)) + 
   facet_wrap(~crop, nrow = 3, scales = "free") +
   geom_vline(xintercept = 0, color = "gray") +
   geom_hline(yintercept = 0, color = "gray") +
   geom_abline(slope = 1, intercept = 0, color = "red") +
   # geom_smooth(method = "lm", se = F, color = "blue") +
   geom_point(size = .5, alpha = 1, shape = 16) +
   geom_text(data = trend_r2 %>% filter(stat == "diff"), 
             aes(label = r2),
             hjust = 0, vjust = 1, color = "blue") +
   theme_bw() +
   theme(axis.text = element_blank()) +
   labs(x = "AgCensus acreage difference, 2012-2017",
        y = "CDL acreage difference, 2012-2017")
ggsave("figures/cdl_agcensus_diff_scatters.png", p, 
       width = 10, height = 5, units = "in")

p <- trends %>% 
   filter(stat == "ratio") %>%
   ggplot(aes(agcensus, cdl)) + 
   facet_wrap(~crop, nrow = 3, scales = "free") +
   geom_vline(xintercept = 1, color = "gray") +
   geom_hline(yintercept = 1, color = "gray") +
   geom_abline(slope = 1, intercept = 0, color = "red") +
   # geom_smooth(method = "lm", se = F, color = "blue") +
   geom_point(size = .5, alpha = 1, shape = 16) +
   geom_text(data = trend_r2 %>% filter(stat == "ratio"), 
             aes(label = r2),
             hjust = 0, vjust = 1, color = "blue") +
   scale_x_log10() +
   scale_y_log10() +
   theme_bw() +
   theme(axis.text = element_blank()) +
   labs(x = "AgCensus acreage ratio, 2012-2017 (log scale)",
        y = "CDL acreage ratio, 2012-2017 (log scale)")
ggsave("figures/cdl_agcensus_ratio_scatters.png", p, 
       width = 10, height = 5, units = "in")

