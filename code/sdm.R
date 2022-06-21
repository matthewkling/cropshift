



clim <- read_csv("data/climate_annual.csv") %>%
  mutate(fips = str_pad(fips, 5, "left", "0"))

clim_bl <- clim %>%
  filter(between(year, 1990, 2010)) %>%
  select(-year) %>%
  group_by(fips) %>%
  summarize_all(mean) %>%
  rename_at(vars(-fips), function(x) paste0(x, "_bl"))

clim_ts <- clim %>%
  rename_at(vars(-fips, -year), function(x) paste0(x, "_ts"))


# starts with creation of g object in cdl.r

gg <- g %>%
  select(fips, year, crop, acres) %>%
  distinct()

exp <- expand_grid(fips = unique(g$fips),
                   crop = unique(g$crop),
                   year = unique(g$year))


gg <- exp %>%
  left_join(gg) %>%
  left_join(clim_bl) %>%
  left_join(clim_ts) %>%
  mutate(acres = ifelse(is.na(acres), 0, acres)) %>%
  group_by(fips, year) %>%
  mutate(ag_acres = sum(acres)) %>%
  group_by(crop) %>%
  mutate(total_crop_acres = sum(acres)) %>%
  ungroup() %>%
  arrange(total_crop_acres) %>%
  select(-total_crop_acres) %>%
  ungroup() %>%
  split(.$crop)

# gg <- gg[1:6]

library(furrr)
plan(multisession, workers = 6)

d <- future_map_dfr(gg, function(x){
  require(mgcv)
  
  message(x$crop[1])
  # x <- gg$corn
  
  x <- rename(x, tmean = tmean_ts, prec = prec_ts)
  
  # time series means (used to fit models)
  m <- x %>%
    group_by(fips) %>%
    summarize(acres = mean(acres),
              ag_acres = mean(ag_acres),
              prop = acres / ag_acres,
              tmean = mean(tmean_bl),
              prec = mean(prec_bl))
  
  # fit model
  fit <- gam(prop ~ te(tmean, prec),
             family = binomial,
             weights = ag_acres,
             data = m) %>%
    suppressWarnings() %>%
    try()
  if(class(fit) == "try-error") return(x)
  
  # predicted proportions
  x$pred <- predict(fit, x, type = "response")
  
  m$pred <- predict(fit, m, type = "response")
  m <- m %>% 
    gather(stat, value, pred, prop) %>%
    mutate(stat = factor(stat, levels = c("prop", "pred"),
                         labels = c("observed", "predicted"))) %>%
    arrange(value)
  p <- ggplot(m, aes(tmean, prec, color = value)) +
    facet_wrap(~ stat, nrow = 1) +
    geom_point() +
    scale_color_viridis_c() +
    theme_bw() +
    labs(color = paste0(x$crop[1], "\nprevalence"))
  tg <- str_replace(x$crop[1], "\\/", "_")
  ggsave(paste0("figures/sdm/cropwise/prevalence_model_", tg, ".png"), 
         p, width = 8, height = 4, units = "in")
  
  return(x)
})



# average suitability of local crop portfolio, weighted by acreage
s <- d %>%
  mutate(pred = as.vector(pred)) %>%
  group_by(crop) %>%
  mutate(pred = pred / max(pred, na.rm = T)) %>%
  group_by(fips, crop) %>%
  mutate(acres_bl = acres[year == min(year)][1]) %>% 
  mutate(year = year - mean(year)) %>%
  group_by(fips, year) %>%
  filter(sum(acres_bl) > 0) %>%
  summarize(pred_bl = weighted.mean(pred, acres_bl),
            pred = weighted.mean(pred, acres),
            acres = sum(acres))

s_lm <- s %>%
  group_by(fips) %>%
  do(broom::tidy(glm(pred ~ year, family = binomial, data = .))) %>%
  select(fips, term, estimate) %>%
  spread(term, estimate) %>%
  janitor::clean_names()

counties <- st_read("data/Counties", "cb_2018_us_county_500k")

s_lm <- counties %>%
  mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
  left_join(s_lm)
  
p <- s_lm %>%
  mutate(year = ifelse(abs(year) > .2, sign(year) * .2, year)) %>%
  ggplot(aes(fill = year)) +
  geom_sf(data = counties, fill = "gray90", color = NA) +
  geom_sf(color = NA) +
  scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
  scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred") %>%
                         rev(),
                       values = c(0, .4, .5, .6, 1),
                       na.value = "white",
                       limits = c(-.2, .2)) +
  theme_void() +
  theme(legend.position = c(.2, .1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 15)) +
  labs(fill = "suitability  \ntrend")
ggsave("figures/sdm/county_trend_map.png", 
       p, width = 10, height = 6, units = "in")


p <- s_lm %>%
  ggplot(aes(fill = intercept)) +
  geom_sf(data = counties, fill = "gray90", color = NA) +
  geom_sf(color = NA) +
  scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
  scale_fill_gradientn(colors = viridis::viridis_pal()(5),
                       values = c(0, .5, .6, .7, 1)) +
  theme_void() +
  theme(legend.position = c(.2, .1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 15)) +
  labs(fill = "suitability  \nintercept")
ggsave("figures/sdm/county_intercept_map.png", 
       p, width = 10, height = 6, units = "in")


# suitability relative to counterfactual
cf <- s %>%
  group_by(fips) %>%
  filter(is.finite(pred),
         is.finite(pred_bl)) %>%
  summarize(diff = sum(pred - pred_bl),
            ratio = mean(log(pred / pred_bl)),
            acres = mean(acres))
cf <- counties %>%
  mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
  left_join(cf)

p <- cf %>%
  mutate(diff = ifelse(abs(diff) > 1, sign(diff) * 1, diff)) %>%
  ggplot(aes(fill = diff)) +
  geom_sf(data = counties, fill = "gray90", color = NA) +
  geom_sf(color = NA) +
  scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
  scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred") %>%
                         rev(),
                       values = c(0, .4, .5, .6, 1),
                       na.value = "white",
                       limits = c(-1, 1)) +
  theme_void() +
  theme(legend.position = c(.2, .1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 15)) +
  labs(fill = "integrated\ndifference    ")
ggsave("figures/sdm/county_counterfactual_map.png", 
       p, width = 10, height = 6, units = "in")

p <- cf %>%
  mutate(ratio = ifelse(abs(ratio) > 1, sign(ratio) * 1, ratio)) %>%
  ggplot(aes(fill = ratio)) +
  geom_sf(data = counties, fill = "gray90", color = NA) +
  geom_sf(color = NA) +
  scale_x_continuous(limits = c(-125, -66), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 50), expand = c(0, 0)) +
  scale_fill_gradientn(colors = c("darkblue", "cornflowerblue", "gray80", "orange", "darkred") %>%
                         rev(),
                       values = c(0, .4, .5, .6, 1),
                       na.value = "white",
                       limits = c(-1, 1)) +
  theme_void() +
  theme(legend.position = c(.2, .1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 15)) +
  labs(fill = "mean\nlog ratio ")
ggsave("figures/sdm/county_counterfactual_map_ratio.png", 
       p, width = 10, height = 6, units = "in")


weighted.mean(cf$diff, cf$acres)

p <- s %>%
  ungroup() %>%
  filter(fips %in% sample(unique(fips), 20)) %>%
  gather(stat, value, pred_bl, pred) %>%
  mutate(stat = factor(stat, levels = c("pred", "pred_bl"),
                       labels = c("actual", "counterfactual"))) %>%
  ggplot(aes(year, value, color = stat)) +
  facet_wrap(~fips, scales = "free", nrow = 5) +
  geom_line() +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top") +
  labs(y = "suitability",
       color = NULL)
ggsave("figures/sdm/county_counterfactual_lines.png", 
       p, width = 10, height = 6, units = "in")





## fit a single multicrop multinomial -- not yet working, SLOW ##

x <- bind_rows(gg) %>%
  filter(fips %in% sample(unique(fips), 500)) %>%
  mutate(crp = as.integer(factor(crop)) - 1)

# time series means (used to fit models)
m <- x %>%
  group_by(fips, crop, crp) %>%
  summarize(acres = mean(acres),
            ag_acres = mean(ag_acres),
            prop = acres / ag_acres,
            aet = mean(aet),
            cwd = mean(cwd),
            gdd = mean(gdd),
            edd = mean(edd),
            ppt = mean(ppt)) %>%
  ungroup()

# m <- filter(m, crp < 10)

# fit model (see ?mgcv::multinom)
K <- max(m$crp)
formulae <- rep("~ s(ppt) + s(gdd)", K)
formulae[1] <- paste("crp", formulae[1])
formulae <- lapply(formulae, as.formula)
fit <- gam(formulae,
           family = multinom(K = K),
           weights = ag_acres,
           data = m)

# predicted proportions
m_pred <- predict(fit, m, type = "response")
x_pred <- predict(fit, x, type = "response")




