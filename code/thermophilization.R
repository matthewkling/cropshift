

library(tidyverse)
library(data.table)

source("e:/ca2cc/ca2cc-county/code/production/ag_census.R")



## munge Ag Census data ##################


d <- fread("e:/ca2cc/ca2cc-county/data/ignore/nass/qs.crops_20210930.txt")

d <- d %>%
  filter(SOURCE_DESC == "CENSUS",
         GROUP_DESC == "FIELD CROPS",
         # GROUP_DESC %in% c("FIELD CROPS", "FRUIT & TREE NUTS", "VEGETABLES"),
         UNIT_DESC %in% c("ACRES", "OPERATIONS")) %>%
  select(SECTOR_DESC:AGG_LEVEL_DESC, 
         STATE_FIPS_CODE, COUNTY_CODE,
         YEAR, VALUE) %>%
  janitor::clean_names()


template <- d %>%
  mutate(year = as.integer(year)) %>%
  filter(
    statisticcat_desc == "AREA HARVESTED",
    prodn_practice_desc %in% c("IRRIGATED", "ALL PRODUCTION PRACTICES"),
    class_desc == "ALL CLASSES",
    agg_level_desc %in% c("NATIONAL", "STATE", "COUNTY")[3],
    domain_desc == "TOTAL",
    year %in% c(2002, 2007, 2012, 2017)[2:4]) %>% ###
  mutate(commodity_desc = ifelse(util_practice_desc == "ALL UTILIZATION PRACTICES",
                                 commodity_desc,
                                 paste0(commodity_desc, "_", util_practice_desc)))# %>%
# group_by(commodity_desc) %>%
# filter(length(unique(util_practice_desc)) == 1 |
#          util_practice_desc == "ALL UTILIZATION PRACTICES") %>%

# group_by(commodity_desc) %>%
# filter(length(unique(year[agg_level_desc == "NATIONAL"])) == 3) # was 4



f <- left_join(template, d) %>%
  select(commodity = commodity_desc, unit = unit_desc, prodn = prodn_practice_desc,
         agg = agg_level_desc, state = state_fips_code, county = county_code,
         year, value) %>%
  filter(agg %in% c("NATIONAL", "STATE", "COUNTY")) %>%
  distinct() %>%
  mutate(unit = str_to_lower(unit))

f <- f %>%
  mutate(prodn = substr(prodn, 1, 3) %>% tolower) %>%
  unite(var, unit, prodn) %>%
  spread(var, value)


# convert implicit missing to explicit missing
f <- f %>% ungroup() %>%
  select(agg, state, county) %>%
  distinct() %>%
  expand_grid(commodity = unique(f$commodity), year = unique(f$year)) %>%
  full_join(f)

# convert values to numeric, handling missing codes
to_numeric <- function(x){
  x[is.na(x)] <- "0"
  x[x == "(Z)"] <- "0"
  x[x == "(D)"] <- NA
  x <- str_remove_all(x, ",")
  x <- as.numeric(x)
  x
}
f <- f %>% mutate(acres_all= to_numeric(acres_all),
                  acres_irr = to_numeric(acres_irr),
                  operations_all = to_numeric(operations_all),
                  operations_irr = to_numeric(operations_irr))


# impute missing values
# f <- f %>%
#   filter(! commodity %in% commodity[agg == "NATIONAL" & is.na(acres_all)]) %>%
#   # mutate(acres_raw = acres_all) %>%
#   group_by(commodity, year) %>%
#   mutate(acres_all = fill_missing(acres_all, operations_all,
#                                   stat = "sum", estimate = "mean",
#                                   level = agg, parent = "NATIONAL", child = "STATE"),
#          acres_irr = fill_missing(acres_irr, operations_irr,
#                                   stat = "sum", estimate = "mean",
#                                   level = agg, parent = "NATIONAL", child = "STATE")) %>%
#   group_by(commodity, year, state) %>%
#   mutate(acres_all = fill_missing(acres_all, operations_all,
#                                   stat = "sum", estimate = "mean",
#                                   level = agg, parent = "STATE", child = "COUNTY"),
#          acres_irr = fill_missing(acres_irr, operations_irr,
#                                   stat = "sum", estimate = "mean",
#                                   level = agg, parent = "STATE", child = "COUNTY"))

f <- f %>%
  filter(agg == "COUNTY") %>%
  select(-agg) %>%
  ungroup() %>% 
  mutate(acres = acres_all - acres_irr,
         acres = ifelse(acres < 0, 0, acres))



## add climate data ##################


# climate <- read_csv("data/ignore/climate.csv") %>%
#   filter(season == "growing_season",
#          start %in% c(1950, 2012)) %>%
#   mutate(gdd = gdd0,
#          start = ifelse(start == 1950, 0, 1)) %>%
#   select(fips, start, prec, gdd) %>%
#   gather(var, value, prec, gdd) %>%
#   unite(var, var, start) %>%
#   spread(var, value) %>%
#   mutate(fips = str_pad(fips, 5, "left", "0"))




climate <- read_csv("e:/ca2cc/ca2cc-county/data/ignore/climate.csv") %>%
  filter(season == "growing_season",
         start == 1950) %>%
  select(fips, 
         county_ppt = prec, county_gdd = gdd, county_aet = aet, 
         county_cwd = cwd, county_edd = edd) %>%
  mutate(fips = str_pad(fips, 5, "left", "0"))

trends <- read_csv("e:/ca2cc/ca2cc-county/data/ignore/climate.csv") %>%
  filter(season == "growing_season",
         start %in% c(1950, 1982)) %>%
  select(fips, year = start, county_ppt = prec_z, county_gdd = gdd_z) %>%
  mutate(fips = str_pad(fips, 5, "left", "0"),
         year = paste0("y", year)) %>%
  gather(var, value, county_ppt, county_gdd) %>%
  spread(year, value) %>%
  mutate(delta = y1982 - y1950) %>%
  select(fips, var, delta) %>%
  mutate(var = paste0(var, "_delta")) %>%
  spread(var, delta)


f <- f %>%
  mutate(fips = paste0(str_pad(state, 2, "left", "0"),
                       str_pad(county, 3, "left", "0"))) %>%
  left_join(climate)



## calculate ag climate indices ##################

write_csv(f, "data/ag_census_processed.csv")
stop()


# marginal climate means for crops and counties
f <- f %>%
  filter(is.finite(county_ppt)) %>%
  filter(is.finite(acres)) %>%
  group_by(commodity) %>%
  mutate(crop_gdd = weighted.mean(county_gdd, acres, na.rm = T),
         crop_ppt = weighted.mean(county_ppt, acres, na.rm = T),
         crop_aet = weighted.mean(county_aet, acres, na.rm = T),
         crop_cwd = weighted.mean(county_cwd, acres, na.rm = T),
         crop_edd = weighted.mean(county_edd, acres, na.rm = T))


f <- f %>% filter(!is.na(acres),
                  !is.na(county_ppt),
                  !is.na(crop_ppt))



library(ggridges)
crop_ridges <- function(f, var, mods = NULL){
  f$value <- f[[var]]
  ff <- f %>% 
    
    # sort
    group_by(commodity) %>% 
    mutate(acres = acres / max(acres),
           mean = weighted.mean(value, acres)) %>%
    ungroup() %>% arrange(mean) %>%
    mutate(commodity = factor(commodity, levels = unique(.$commodity))) %>%
    
    group_by(commodity) %>% # sample, so that density geom can be used
    sample_n(10000, replace = T, weight = acres)
  
  p <- ff %>%
    ggplot(aes(value, commodity, group = commodity)) +
    geom_density_ridges(bandwidth = diff(range(ff$value))/ 20) +
    labs(y = NULL, x = var) + 
    theme_minimal()
  
  if(!is.null(mods)) p <- p + mods
  
  ggsave(paste0("figures/ridges_", var, ".png"), 
         p, width = 8, height = 5, units = "in")
}

ff <- f %>% filter(! str_detect(commodity, "TARO"), 
                   !str_detect(commodity, "GRASSES"))

ff %>% crop_ridges("county_aet", xlim(0, 900))
ff %>% crop_ridges("county_cwd", xlim(0, 1200))
ff %>% crop_ridges("county_ppt", xlim(0, NA))
ff %>% crop_ridges("county_gdd", xlim(0, NA))
ff %>% crop_ridges("county_edd", xlim(0, 200))




county <- f %>% 
  group_by(fips, year) %>%
  summarize(crop_gdd = weighted.mean(crop_gdd, acres),
            crop_ppt = weighted.mean(crop_ppt, acres),
            crop_aet = weighted.mean(crop_aet, acres),
            crop_cwd = weighted.mean(crop_cwd, acres),
            crop_edd = weighted.mean(crop_edd, acres),
            acres = sum(acres))

crop <- f %>% 
  group_by(commodity, year) %>%
  summarize(county_gdd = weighted.mean(county_gdd, acres),
            county_ppt = weighted.mean(county_ppt, acres),
            county_aet = weighted.mean(county_aet, acres),
            county_cwd = weighted.mean(county_cwd, acres),
            county_edd = weighted.mean(county_edd, acres),
            acres = sum(acres))


crop_pc <- crop %>%
  ungroup() %>%
  select(county_gdd:county_edd) %>%
  # select(county_aet, county_cwd) %>%
  filter(is.finite(county_aet)) %>%
  scale() %>%
  prcomp()

library(ggfortify)
p <- autoplot(crop_pc, loadings = T, loadings.label = TRUE)
ggsave("figures/crop_pca.png", 
       p, width = 6, height = 6, units = "in")


p <- crop %>%
  group_by(commodity) %>%
  summarize(county_ppt = mean(county_aet, na.rm = T),
            county_edd = mean(county_cwd, na.rm = T)) %>%
  filter(!str_detect(commodity, "OTHER")) %>%
  ggplot(aes(county_ppt, county_edd, label = commodity)) +
  geom_text() +
  scale_y_log10() +
  scale_x_continuous(expand = c(.25, .25)) +
  theme_classic() +
  labs(x = "mean growing season AET", y = "mean growing season CWD")
ggsave("figures/crop_scatter.png", 
       width = 5, height = 5, units = "in")


p <- f %>%
  group_by(fips) %>%
  sample_n(1) %>%
  ggplot(aes(county_aet, county_cwd)) +
  geom_point() +
  scale_y_log10() +
  # scale_x_continuous(expand = c(.25, .25)) +
  theme_classic() +
  labs(x = "mean growing season AET", y = "mean growing season CWD")
ggsave("figures/county_scatter.png", 
       width = 5, height = 5, units = "in")


p <- county %>%
  ggplot(aes(crop_aet, crop_cwd)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  # scale_x_continuous(expand = c(.25, .25)) +
  theme_classic() +
  labs(x = "mean growing season AET", y = "mean growing season CWD")
ggsave("figures/county_crop_scatter.png", 
       width = 5, height = 5, units = "in")




## temporal trends ##############


d1 <- county %>% 
  group_by(year) %>% 
  summarize(crop_gdd = weighted.mean(crop_gdd, acres),
            crop_ppt = weighted.mean(crop_ppt, acres),
            crop_aet = weighted.mean(crop_aet, acres),
            crop_cwd = weighted.mean(crop_cwd, acres),
            crop_edd = weighted.mean(crop_edd, acres)) %>%
  gather(var, value, -year)

d2 <- crop %>% 
  group_by(year) %>% 
  summarize(county_gdd = weighted.mean(county_gdd, acres),
            county_ppt = weighted.mean(county_ppt, acres),
            county_aet = weighted.mean(county_aet, acres),
            county_cwd = weighted.mean(county_cwd, acres),
            county_edd = weighted.mean(county_edd, acres)) %>%
  gather(var, value, -year)

bind_rows(d1 %>% mutate(margin = "county"), 
          d2 %>% mutate(margin = "crop")) %>%
  mutate(var = str_remove(var, "crop_|county_")) %>%
  ggplot(aes(year, value, color = margin)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ var, scales = "free", nrow = 1) +
  scale_x_continuous(breaks = range(crop$year)) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(y = NULL, x = NULL)
ggsave("figures/global_marginal_trends.png", 
       width = 8, height = 3, units = "in")



paths <- county %>%
  ungroup() %>%
  filter(year %in% range(year)) %>%
  # filter(fips %in% sample(unique(fips), 500)) %>%
  group_by(fips) %>%
  arrange(year) %>%
  mutate(dx = crop_aet[2] - crop_aet[1],
         dy = crop_cwd[2] - crop_cwd[1]) %>% ungroup() %>%
  filter(is.finite(dx))
paths$hex <- colormap::colors2d(select(paths, dx, dy) %>% as.data.frame(),
                                xtrans = "rank", ytrans = "rank")


p <- ggplot(paths, aes(crop_aet, crop_cwd, group = fips,
                       color = hex)) +
  geom_path(arrow = grid::arrow(length = unit(.1, "in"),
                                type = "closed", angle = 15)) +
  scale_color_identity() +
  theme(panel.background = element_rect(fill = "black"),
        panel.grid = element_blank())
ggsave("figures/county_vectors.png", 
       p, width = 8, height = 8, units = "in")



county_lm <- county %>% ungroup() %>%
  gather(var, value, crop_gdd:crop_edd) %>%
  filter(is.finite(value)) %>%
  group_by(var) %>%
  mutate(value = value / sd(value)) %>%
  group_by(fips, var) %>%
  do(broom::tidy(lm(value ~ year, ., weights = acres))) %>%
  select(fips, var, term, estimate) %>%
  spread(term, estimate) %>%
  janitor::clean_names()

crop_lm <- crop %>% ungroup() %>%
  gather(var, value, county_gdd:county_edd) %>%
  filter(is.finite(value)) %>%
  group_by(var) %>%
  mutate(value = value / sd(value)) %>%
  group_by(commodity, var) %>%
  do(broom::tidy(lm(value ~ year, ., weights = acres))) %>%
  select(commodity, var, term, estimate) %>%
  spread(term, estimate) %>%
  janitor::clean_names()



trans_cbrt <- scales::trans_new(
  name = "trans_cbrt",
  transform = function(x) abs(x) ^ (1/3) * sign(x),
  inverse = function(x) x ^ 3
)



library(sf)
cd <- readRDS("e:/ca2cc/county_explorer/codex/data/counties.rds") %>%
  st_as_sf() %>%
  left_join(climate) %>%
  left_join(county_lm) %>%
  filter(!is.na(var))

p <- cd %>%
  filter(var %in% c("crop_aet", "crop_cwd")) %>%
  ggplot(aes(fill = year * 10)) +
  facet_wrap(~ var, ncol = 1) +
  geom_sf(color = NA) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white")) +
  scale_fill_gradientn(colors = c("orange", "red", "black", "blue", "dodgerblue"),
                       trans = trans_cbrt,
                       breaks = c(-10, -1, 0, 1, 10)) +
  labs(fill = "sd / decade")
ggsave("figures/county_lm_maps.png", 
       p, width = 6, height = 8, units = "in")



p <- readRDS("e:/ca2cc/county_explorer/codex/data/counties.rds") %>%
  st_as_sf() %>%
  left_join(county_lm %>%
              select(-intercept) %>%
              spread(var, year)) %>%
  filter(is.finite(crop_aet)) %>%
  mutate(hex = colormap::colors2d(bind_cols(.$crop_aet, .$crop_cwd) %>% as.data.frame(),
                            xtrans = "rank", ytrans = "rank")) %>%
  ggplot(aes(fill = hex)) +
  geom_sf(color = NA) +
  scale_fill_identity() +
  theme_void()
ggsave("figures/county_lm_maps_2d.png", 
       p, width = 6, height = 4, units = "in")



p <- cd %>%
  filter(var == var[1]) %>% select(-var) %>%
  gather(var, value, county_ppt:county_edd) %>%
  group_by(var) %>%
  filter(is.finite(value)) %>%
  mutate(value = (value - mean(value)) / sd(value)) %>%
  
  ggplot(aes(fill = value)) +
  facet_wrap(~ var, ncol = 2) +
  geom_sf(color = NA) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white")) +
  scale_fill_viridis_c() +
  labs(fill = "growing\nseason\nclimate")
ggsave("figures/county_climate_maps.png", 
       p, width = 8, height = 6, units = "in")




cd <- readRDS("e:/ca2cc/county_explorer/codex/data/counties.rds") %>%
  st_as_sf() %>%
  left_join(county %>% group_by(fips) %>%
              summarise(aet = weighted.mean(crop_aet, acres, na.rm = T),
                        cwd = weighted.mean(crop_cwd, acres, na.rm = T)))


p <- cd %>%
  gather(var, value, aet, cwd) %>%
  group_by(var) %>%
  filter(is.finite(value)) %>%
  mutate(value = (value - mean(value)) / sd(value)) %>%
  
  ggplot(aes(fill = value)) +
  facet_wrap(~ var, ncol = 1) +
  geom_sf(color = NA) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white")) +
  scale_fill_viridis_c() +
  labs(fill = "growing\nseason\nclimate")
ggsave("figures/county_crop_climate_maps.png", 
       p, width = 8, height = 6, units = "in")



crop_ord <- f %>% group_by(commodity) %>% 
  summarize(aet = weighted.mean(crop_aet, acres, na.rm = T)) %>%
  ungroup() %>%
  arrange(aet)
county_ord <- f %>% group_by(fips) %>% 
  summarize(aet = weighted.mean(county_aet, acres, na.rm = T)) %>%
  ungroup() %>%
  arrange(aet)
x <- f %>%
  group_by(commodity, fips) %>% 
  summarize(county_aet = mean(county_aet, na.rm = T),
            crop_aet = mean(crop_aet, na.rm = T),
            acres = sum(acres, na.rm = T)) %>%
  ungroup() %>%
  mutate(commodity = factor(commodity, levels = crop_ord$commodity),
         fips = factor(fips, levels = county_ord$fips)) %>%
  filter(acres > 0)

p <- x %>%
  ggplot(aes(fips, commodity, fill = acres)) +
  geom_raster() +
  scale_fill_gradientn(colors = c("gray80", "orange", "darkred"),
                       trans = "log10") +
  theme_void() +
  theme(legend.position = "none")
ggsave("figures/county_crop_heatmap.png", 
       p, width = 8, height = 5, units = "in")



x <- f %>%
  mutate(commodity = factor(commodity, levels = crop_ord$commodity),
         fips = factor(fips, levels = county_ord$fips)) %>%
  filter(acres > 0)

p <- x %>%
  group_by(commodity) %>%
  mutate(acres = acres / max(acres, na.rm = T)) %>%
  ggplot(aes(fips, paste(commodity, year), 
             fill = factor(year), alpha = acres)) +
  geom_raster() +
  scale_alpha_continuous(range = c(.5, 1), trans = "sqrt") +
  scale_fill_manual(values = c("black", "red", "dodgerblue")) +
  theme_void() +
  theme(legend.position = "none")
ggsave("figures/county_crop_heatmap_yr.png", 
       p, width = 8, height = 5, units = "in")



crop_lm %>%
  mutate(var = str_remove(var, "crop_|county_")) %>%
  ggplot(aes(year, commodity, color = sign(year) == 1)) +
  facet_grid(. ~ var, scales = "free") +
  geom_vline(xintercept = 0) +
  geom_point()

