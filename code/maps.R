
library(tidyverse)
library(rgdal)
library(sf)

d <- read_csv("data/output/wide_data.csv")

dd <- d %>%
  select(fips,
         AgCensus_cover_crop_acres_2017,
         AgCensus_tile_drainage_acres_2017,
         AgCensus_no_till_acres_2017,
         # AgCensus_cropland_acres_2017,
         AgCensus_acres_operated_2017,
         climate_tmean_full_year_1950_1980) %>%
  rename(cover_crops = AgCensus_cover_crop_acres_2017,
         tile_drainage = AgCensus_tile_drainage_acres_2017,
         no_till = AgCensus_no_till_acres_2017,
         # cl = AgCensus_cropland_acres_2017,
         cl = AgCensus_acres_operated_2017,
         temp = climate_tmean_full_year_1950_1980) %>%
  gather(var, value, cover_crops, tile_drainage, no_till) %>%
  mutate(value = value / cl) %>%
  group_by(var) %>%
  mutate(rank = ecdf(value)(value))

counties <- readOGR("data/Counties", "cb_2018_us_county_500k") %>%
  st_as_sf()

counties <- counties %>%
  mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
  left_join(dd) %>%
  filter(!is.na(var))

p <- ggplot(counties, aes(fill = rank)) +
  facet_wrap(~var, ncol = 1) +
  geom_sf(color = NA) +
  xlim(-130, -60) +
  ylim(25, 51) +
  scale_fill_viridis_c() +
  theme_void() +
  labs(fill = "rank")
ggsave("../for meredith/adapt_maps_2017.png", width = 6, height = 8, units = "in")

p <- ggplot(counties %>% filter(var != "no_till"), 
            aes(fill = value * 100)) +
  facet_wrap(~var, ncol = 1) +
  geom_sf(color = NA) +
  xlim(-130, -60) +
  ylim(25, 51) +
  scale_fill_viridis_c(trans = "log10", 
                       breaks = c(.01, .1, 1, 10, 100),
                       labels = c(.01, .1, 1, 10, 100),
                       na.value = "gray20") +
  theme_void() +
  labs(fill = "% ag land")
ggsave("../for meredith/adapt_maps_2017_v2.png", width = 6, height = 6, units = "in")


ggplot(counties %>% filter(var == var[1]), 
       aes(fill = temp)) +
  geom_sf(color = NA) +
  xlim(-130, -60) +
  ylim(25, 51) +
  scale_fill_viridis_c() +
  theme_void() +
  labs(fill = "rank")
