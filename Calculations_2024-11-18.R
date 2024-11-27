# Data load ---------------------------------------------------------------
library(tidyverse)
theme_set(
    theme_bw() +
        theme(
            text = element_text(family = "serif", size = 15),
            legend.position = "bottom"
        )
)
L <- c("фон_2009", "фон_2014", "буф_2009", "буф_2014", "имп_2009", "имп_2014")
long <- readxl::read_excel("Data/Carabidae_25.01.2023.xlsx", 
                           sheet = "main_data") %>% 
    select(-duration, -traps) %>% 
    arrange(desc(site), taxa) %>% # -taxa
    mutate(tur = factor(tur), 
           site = fct_inorder(site),
           zone = substr(zone, 1, 3),
           zone4 = factor(zone, levels = c("fon", "buf", "imp", "sup")),
           zone3 = fct_collapse(zone4, imp = c("imp", "sup")),
           km = str_extract(site, "[:digit:]{1,}"),
           km = as.numeric(km), 
           year = as.factor(year), 
           .after = "site")

# 1 = turs 1 and 2 are united
wide1 <- long %>% # div2 - turs are united
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    select(-no_insects)
# 1 = turs 1 and 2 are united
div1 <- tibble(wide1[,1:7], 
               abu = apply(wide1[,8:ncol(wide1)], 1, sum),
               nsp = apply(wide1[,8:ncol(wide1)], 1, function(a){length(a[a>0])}),
               shan= vegan::diversity(wide1[,8:ncol(wide1)], 
                                      MARGIN = 1, index = "shannon", base = exp(1))
) %>% 
    mutate(km2 = km^2, kmLog = log(km), abuLog = log(abu+1), .after = km)


# Models template ---------------------------------------------------------
models_fit <- function(formulas){
fits <- tibble(ff = formulas) %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){ 
        if(str_detect(a$ff, "Segmented")){
            #######
            segmented::segmented(lm(
                str_replace_all(a$ff, "Segmented", ""), 
                data = div1), seg.Z = ~km)
            #######
        } else {
            lm(formula = a$ff, data = div1)
        }
    }) %>% 
    tibble(ff = formulas, fit = .)
pred <- expand_grid(year = factor(c(2009, 2014)), km = seq(1, 32, by = 0.5)) %>%
    mutate(
        km2 = km^2, 
        kmLog = log(km), 
        zone = case_when(km < 3 ~ "sup", km < 7 ~ "imp", 
                         km < 22 ~ "buf", TRUE ~ "fon"), 
        zone4 = factor(zone, levels = c("fon", "buf", "imp", "sup")),
        zone3 = fct_collapse(zone4, imp = c("imp", "sup")),
    ) %>% 
    mutate(
        zonal3    = predict(fits$fit[[1]], .),
        zonal4    = predict(fits$fit[[2]], .),
        linear    = predict(fits$fit[[3]], .),
        segmented = predict(fits$fit[[4]], .),
        nonlinear = predict(fits$fit[[5]], .), 
        logarithm = predict(fits$fit[[6]], .)
    ) %>% 
    select(-km2:-zone3)
if(str_detect(toupper(formulas[1]), "LOG")) {
    pred <- pred %>% 
        mutate_at(3:ncol(.), exp)
}
pred <- pred %>% 
    `colnames<-`(c("year", "km", formulas)) %>% 
    pivot_longer(names_to = "model", values_to = "abu", -c("year", "km"))
pred$model <-  map_chr(str_split(pred$model, "year \\+ "), ~.x[2])
pred <- pred %>% 
    mutate(model = factor(model, 
        levels = c("zone3", "zone4", "km", "kmSegmented", "km + km2", "kmLog"))
    ) %>% 
    split(.$model)

fits %>% 
    pull(fit) %>% 
    lapply(function(a){
        tibble(
            aic = AIC(a), 
            r2 = summary(a)$adj.r.squared, 
            sh.test = round(shapiro.test(a$residuals)$p.value, 4)
        ) 
    }) %>% 
    map_df(rbind) %>% 
    tibble(fits, .) %>% 
    mutate(d = pred)
}

model_viz <- function(df, yy){
    div1 %>% 
        rename(yy = which(colnames(div1)==yy)) %>% 
        ggplot(aes(x = km, y = yy)) + 
        geom_line(aes(km, abu, color = model), data = map_dfr(df$d, rbind), 
                  linewidth = 1) + 
        geom_point(shape = 21, size = 3) + 
        facet_wrap(~year, scales = "fixed", ncol = 2) +
        theme(panel.grid = element_blank()) +
        guides(fill="none") 
}

res <- list()

# Abundance ---------------------------------------------------------------
res$abundance <- c("zone3", "zone4", "km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("abu ~ year + ", .) %>% 
    models_fit()
model_viz(res$abundance, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели исходных показателей обилия")

ggsave("1a. Abundance.png", height = 8, width = 11, dpi = 600)
res$abundance

# Abundance LLOG -----------------------------------------------------------
res$abundance_log <- c("zone3", "zone4", "km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("abuLog ~ year + ", .) %>% 
    models_fit()
model_viz(res$abundance_log, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели логарифмированных показателей обилия")

ggsave("1b. Abundance_log.png", height = 8, width = 11, dpi = 600)
res$abundance_log

# N_species ---------------------------------------------------------------
res$nsp <- c("zone3", "zone4", "km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("nsp ~ year + ", .) %>% 
    models_fit()
model_viz(res$nsp, "nsp") + 
    labs(x = NULL, y = "Количество видов",
         subtitle = "Видовое богатство")

ggsave("2. n_species.png", height = 8, width = 11, dpi = 600)
res$nsp

# Shannon -----------------------------------------------------------------
res$shan <- c("zone3", "zone4", "km", "kmSegmented", "km + km2", "kmLog") %>% 
    paste0("shan ~ year + ", .) %>% 
    models_fit()
model_viz(res$shan, "shan") + 
    labs(x = NULL, y = "Индекс Шеннона",
         subtitle = "Видовое разнообразие")

ggsave("3. Shannon.png", height = 8, width = 11, dpi = 600)
res$shan

# final export ------------------------------------------------------------
res %>% 
    map(~select(.x, -d, -fit)) %>% 
    writexl::write_xlsx(paste0("models_", Sys.Date(), ".xlsx"))




