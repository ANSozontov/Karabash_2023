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
pred <- seq(1, 32, by = 0.5) %>% 
    c(fits$fit[[2]]$psi[2]) %>% 
    sort() %>% 
    unique() %>%
    expand_grid(
        year = factor(c(2009, 2014)), 
        km = .) %>%
    mutate(
        km2 = km^2, 
        kmLog = log(km)
        # zone = case_when(km < 3 ~ "sup", km < 7 ~ "imp", 
        #                  km < 22 ~ "buf", TRUE ~ "fon"), 
        # zone4 = factor(zone, levels = c("fon", "buf", "imp", "sup")),
        # zone3 = fct_collapse(zone4, imp = c("imp", "sup")),
    ) %>% 
    mutate(
        linear    = predict(fits$fit[[1]], .),
        segmented = predict(fits$fit[[2]], .),
        nonlinear = predict(fits$fit[[3]], .), 
        logarithm = predict(fits$fit[[4]], .)
    ) %>% 
    select(-km2, -kmLog)
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
        levels = c("km", "kmSegmented", "km + km2", "kmLog")) # "zone3", "zone4", 
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
res$abundance <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4", 
    paste0("abu ~ year + ", .) %>% 
    models_fit()
res$abundance
model_viz(res$abundance, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели исходных показателей обилия")

ggsave("1a. Abundance.png", height = 8, width = 11, dpi = 600)

# Abundance LLOG -----------------------------------------------------------
res$abundance_log <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4", 
    paste0("abuLog ~ year + ", .) %>% 
    models_fit()
res$abundance_log
model_viz(res$abundance_log, "abu") + 
    labs(x = NULL, y = "Обилие (особей на 100 лов.-сут.)", 
         subtitle = "Модели логарифмированных показателей обилия")
ggsave("1b. Abundance_log.png", height = 8, width = 11, dpi = 600)

# N_species ---------------------------------------------------------------
res$nsp <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4", 
    paste0("nsp ~ year + ", .) %>% 
    models_fit()
res$nsp
model_viz(res$nsp, "nsp") + 
    labs(x = NULL, y = "Количество видов",
         subtitle = "Видовое богатство")

ggsave("2. n_species.png", height = 8, width = 11, dpi = 600)

# Shannon -----------------------------------------------------------------
res$shan <- c("km", "kmSegmented", "km + km2", "kmLog") %>% # "zone3", "zone4",
    paste0("shan ~ year + ", .) %>% 
    models_fit()
res$shan
model_viz(res$shan, "shan") + 
    labs(x = NULL, y = "Индекс Шеннона",
         subtitle = "Видовое разнообразие")
ggsave("3. Shannon.png", height = 8, width = 11, dpi = 600)

# final export ------------------------------------------------------------
res %>% 
    map(~select(.x, -d, -fit)) %>% 
    writexl::write_xlsx(paste0("models_", Sys.Date(), ".xlsx"))




