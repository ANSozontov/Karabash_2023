# Data load ---------------------------------------------------------------
library(tidyverse)
theme_set(
    theme_bw() +
        theme(
            text = element_text(family = "serif", size = 15),
            legend.position = "bottom"
        )
)
manual.aic <- function(a){
    if(!(a$family$family %in% c("quasipoisson", "poisson"))){return(NA)}
    loglik <- sum(dpois(a$model[,1], fitted.values(a), log = TRUE))
    phi <- summary(a)$dispersion
    -2 * loglik + 2 * summary(a)$df[3] * phi
}

long <- readxl::read_excel("Data/Carabidae_25.01.2023.xlsx", 
                           sheet = "main_data") %>% 
    select(-duration, -traps) %>% 
    arrange(desc(site), taxa) %>% # -taxa
    mutate(tur = factor(tur), 
           site = fct_inorder(site),
           zone = case_when(zone == "fon" ~ "фоновая", 
                            zone == "bufer" ~ "буферная",
                            zone == "superimpact" ~ "суперимпактная", 
                            TRUE ~ "импактная"),
           zone = factor(zone, levels = c("фоновая", "буферная", "импактная", "суперимпактная")),
           km = str_extract(site, "[:digit:]{1,}"),
           km = as.numeric(km), 
           year = as.factor(year), 
           .after = "site")

# 2 = tur 1 and tur 2 are separated
wide2 <- long %>% 
    select(-num) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    select(-no_insects)

# 2 = tur 1 and tur 2 are separated
div2 <- tibble(wide2[,1:6], 
              abu = apply(wide2[,7:ncol(wide2)], 1, sum),
              nsp = apply(wide2[,7:ncol(wide2)], 1, function(a){length(a[a>0])}),
              shan= vegan::diversity(wide2[,7:ncol(wide2)], 
                                     MARGIN = 1, index = "shannon", base = exp(1))
)

# 1 = turs 1 and 2 are united
wide1 <- long %>% # div2 - turs are united
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    select(-no_insects)
# 1 = turs 1 and 2 are united
div1 <- tibble(wide1[,1:5], 
               abu = apply(wide1[,6:ncol(wide1)], 1, sum),
               nsp = apply(wide1[,6:ncol(wide1)], 1, function(a){length(a[a>0])}),
               shan= vegan::diversity(wide1[,6:ncol(wide1)], 
                                      MARGIN = 1, index = "shannon", base = exp(1))
)

# Models selection ---------------------------------------------------
fits <- expand_grid(
            yy = c("abu ", "nsp ", "shan "),
            xx = paste0("year"), 
            zz = c(" + km", " + zone"), # " + site", 
            zero_filter = c("0_incl", "0_excl"),
            distr = c("poisson", "quasipoisson",  "gaussian"),
            tt = c("", "-t", " * tur", " + tur"),
    ) %>% 
    filter(yy == "shan " | zero_filter == "0_incl", 
           yy != "shan "  | distr == "gaussian") %>% 
    mutate(
          lg = case_when(distr != "gaussian" ~ FALSE, 
                         yy == "shan " ~ FALSE,
                         TRUE ~ TRUE),
          xx = paste0(xx, tt, zz), #.keep = "unused",
          turs.united = str_detect(tt, "-t")
          ) %>% 
    unite(ff, yy, xx, sep = "~ ") %>%  
    transmute(
        ff, zero_filter, distr, lg, turs.united,
        nm = paste0(ff, "_", zero_filter, "_", distr),
    )

fits <- fits %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){
            if(a$turs.united){
                df <- div1
                a$ff <- str_replace_all(a$ff, "-t", "")
            } else {
                df <- div2
            }
            if(a$zero_filter == "0_excl"){ 
                  df <- filter(df, shan > 0)
            }
            if(a$lg){df <- mutate(df, abu = log(abu+1))}
            glm(formula = a$ff, data = df, family = a$distr)
      }) %>% 
      `names<-`(fits$nm) %>% 
      mutate(fits, fit = ., info = lapply(., summary)) 

fits <- fits %>% 
      pull(fit) %>% 
      lapply(function(a){ 
            resid = residuals(a)
            ftted = fitted.values(a)
            if(a$family$family == "poisson"){
                  suppressMessages(overdispersion <- a |> 
                                         performance::check_overdispersion() |>
                                         capture.output() |> 
                                         str_subset("value") |> 
                                         str_replace_all("[:alpha:]|[:blank:]|[:symbol:]|-", "") |> 
                                         as.numeric())
            } else {overdispersion <- NA}
            
            tibble(aic = AIC(a), 
                   aic2 = manual.aic(a),
                   shapiro = shapiro.test(resid)$p.value,
                   overdispersion = overdispersion,
                   var.res = var(resid), 
                   var.ftd = var(ftted), 
                   rows = nrow(a$data),
                   coef = a %>% 
                       coefficients %>% 
                       data.frame(val = .) %>% 
                       rownames_to_column("var") %>% 
                       filter(str_detect(var, "year|tur")) %>% 
                       mutate(val = round(val, 2)) %>% 
                       unite("coef", sep = " = ") %>% 
                       pull(coef) %>% 
                       paste0(collapse = "; ")
            ) %>% 
                  mutate(r2 = var.ftd/(var.res + var.ftd)) %>% 
                  return()
      }) %>% 
      map_df(rbind, .id = "nm") %>% 
      left_join(fits, by = "nm") %>% 
      mutate(aic = case_when(is.na(aic) ~ aic2, TRUE ~ aic), .keep = "unused")

# Models export ---------------------------------------------------
fits %>% 
    filter(
        ff %in% c("abu ~ year + tur + site", "abu ~ year + site", 
                  "nsp ~ year-t + site", "shan ~ year-t + site"), 
        str_detect(nm, "shan", negate = TRUE) | zero_filter == "0_excl",
        distr == "gaussian"
    ) %>% 
    pull(fit) %>% 
    map(~.x %>% 
            summary %>% 
            `$`(coefficients) %>% 
            as.data.frame %>% 
            rownames_to_column("predictor") %>% 
            mutate(Коэффициент = round(Estimate, 2),
                   Ошибка = round(`Std. Error`, 2),
                   t = round(`t value`, 1), 
                   p.value = round(`Pr(>|t|)`, 3),
                   p.value = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ as.character(p.value)),
                   .keep = "unused"
            )
    ) %>% 
    `names<-`(str_replace_all(names(.), "-t", "")) %>% 
    writexl::write_xlsx(paste0("Models_selected_", Sys.Date(), ".xlsx"))

fits %>% 
    mutate(pred = substr(ff, 1, 4)) %>% 
    split(.$pred) %>% 
    map(~.x %>% transmute(
        `Туры объединены` = turs.united, 
        `Нулевые пробы` = zero_filter,
        `Cтрок (наблюдений) в исходных данных` = rows,
        `Распределение` = distr, 
        `Логарифмирование` = lg, 
        Formula = ff, 
        aic, 
        r2, 
        `Нормальность остатков` = round(shapiro, 3), 
        overdispersion,
        Coefficients = coef,
        )
    ) %>% 
    `names<-`(c("Модели обилия_ВСЕ", "Модели вид.бог._ВСЕ", "Модели разн.Шен._ВСЕ")) %>%
    # append(fits.abu.sel) %>% # "Модели по обилию_Изб.1", "Модели по обилию_Изб.2"
    writexl::write_xlsx(paste0("Models_all_", Sys.Date(), ".xlsx"))

# Models Viz --------------------------------------------------------------
custom_viz <- function(otkl, pred, yy, tt, df){
    if(df == 2){
        df <- div2
    } else if(df == 1){
        df <- mutate(div1, tur = 0)
    }
    p <- df %>% 
        mutate(
            km2 = case_when(substr(site, 4, 4) == "N" ~ substr(site, 2, 3), 
                            TRUE ~ paste0("-", substr(site, 2, 3))),
            km2 = as.numeric(km2),
            km1 = as.numeric(substr(site, 2, 3)), 
            type = case_when(
                tur == 1 ~ paste0("1. Начало лета, ", year),
                tur == 2 ~ paste0( "2. Конец лета, ", year),
                TRUE ~ year)
        ) %>% 
        select(type, otkl = all_of(otkl), pred = all_of(pred)) %>% 
        group_by(type, pred) %>% 
        summarise(MED = median(otkl), 
                  MAX = max(otkl),
                  MIN = min(otkl),
                  Q25 = quantile(otkl, 0.25), 
                  Q75 = quantile(otkl, 0.75), 
                  .groups = "drop") %>% 
        ggplot(aes(x = pred, y = MED, 
                   color = pred, )) + 
        geom_vline(xintercept = 0, color = "coral", 
                   linewidth = 2, alpha = 0.5, # linetype = "dashed"
        ) + 
        geom_linerange(aes(ymin = MIN, ymax = MAX), linewidth = 0.7) + 
        geom_crossbar(aes(ymin = Q25, ymax = Q75), fill = "white", linewidth = 0.7) +
        facet_wrap(~type, scales = "free", ncol = 2) +
        theme(panel.grid = element_blank()) +
        labs(x = NULL, y = yy, #"Обилие (особей на 100 лов.-сут.)", 
             subtitle = tt)+ #"Обилие в разные годы и туры"
        guides(color="none") + 
        theme(strip.background=element_rect(fill= "white"))
    
    if(pred == "km1"){
        p + scale_colour_gradient2(low = "red", mid = "darkgoldenrod4", midpoint = 15, high = "green")
    } else {
        p + scale_colour_gradient2(low = "green", mid = "red", midpoint = 0, high = "green")
    }
}

p <- expand_grid(
        otkl = c("abu", "nsp", "shan"), 
        pred = c("km1", "km2")) %>% 
    mutate(yy = case_when(otkl == "abu" ~ "Обилие (особей на 100 лов.-сут.)",
                          otkl == "nsp" ~ "Богатство (количество видов)", 
                          TRUE ~ "Разнообразие (индекс Шеннона)"), 
           tt = case_when(otkl == "abu" ~ "Обилие в разные годы и туры",
                          otkl == "nsp" ~ "Видовое богатство в разные годы", 
                          TRUE ~ "Видовое разнообразие в разные годы"), 
           df = case_when(otkl == "abu" ~ 2, TRUE ~ 1),
           nm = paste0(otkl, "_", pred, ".png")) %>% 
    split(.$nm) %>% 
    map(~custom_viz(.x$otkl, .x$pred, .x$yy, .x$tt, .x$df))

map2(p[1:2], names(p)[1:2], 
     ~ggsave(filename = .y, plot = .x, width = 10, height = 7, dpi = 900), 
     .progress = TRUE
)
map2(p[3:6], names(p)[3:6], 
     ~ggsave(filename = .y, plot = .x, width = 10, height = 3.5, dpi = 900), 
     .progress = TRUE
)


# N.species models selection ---------------------------------------------------
fits.nsp <- expand_grid(
    yy = "nsp ",
    xx = paste0("year"), 
    zz = c(" + site", " + km", " + zone"),
    super_imp = c("S.I.", ""),
    distr = c("poisson", "quasipoisson",  "gaussian"),
    tt = c("", "-t", " * tur", " + tur"),
) %>% 
    filter(super_imp == "S.I.") %>% 
    mutate(
        lg = case_when(distr != "gaussian" ~ FALSE, TRUE ~ TRUE),
        xx = paste0(xx, tt, zz), #.keep = "unused",
        turs.united = str_detect(tt, "-t")
    ) %>% 
    unite(ff, yy, xx, sep = "~ ") %>%  
    transmute(
        ff, super_imp, distr, lg, turs.united,
        nm = paste0(ff, "_", super_imp, "_", distr),
    )

fits.nsp <- fits.nsp %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){
        if(a$turs.united){
            df <- div1
            a$ff <- str_replace_all(a$ff, "-t", "")
        } else {
            df <- div2
        }
        if(a$super_imp != "S.I."){ 
            df <- filter(df, zone != "суперимпактная")
        }
        if(a$lg){df <- mutate(df, nsp = log(nsp+1))}
        glm(formula = a$ff, data = df, family = a$distr)
    }) %>% 
    `names<-`(fits.nsp$nm) %>% 
    mutate(fits.nsp, fit = ., info = lapply(., summary)) 

fits.nsp <- fits.nsp %>% 
    pull(fit) %>% 
    lapply(function(a){ 
        resid = residuals(a)
        ftted = fitted.values(a)
        if(a$family$family == "poisson"){
            suppressMessages(overdispersion <- a |> 
                                 performance::check_overdispersion() |>
                                 capture.output() |> 
                                 str_subset("value") |> 
                                 str_replace_all("[:alpha:]|[:blank:]|[:symbol:]|-", "") |> 
                                 as.numeric())
        } else {overdispersion <- NA}
        
        tibble(aic = AIC(a), 
               aic2 = manual.aic(a),
               shapiro = shapiro.test(resid)$p.value,
               overdispersion = overdispersion,
               var.res = var(resid), 
               var.ftd = var(ftted), 
               rows = nrow(a$data),
               coef = a %>% 
                   coefficients %>% 
                   data.frame(val = .) %>% 
                   rownames_to_column("var") %>% 
                   filter(str_detect(var, "year|tur")) %>% 
                   mutate(val = round(val, 2)) %>% 
                   unite("coef", sep = " = ") %>% 
                   pull(coef) %>% 
                   paste0(collapse = "; ")
        ) %>% 
            mutate(r2 = var.ftd/(var.res + var.ftd)) %>% 
            return()
    }) %>% 
    map_df(rbind, .id = "nm") %>% 
    left_join(fits.nsp, by = "nm") %>% 
    mutate(aic = case_when(is.na(aic) ~ aic2, TRUE ~ aic), .keep = "unused")

# writexl::write_xlsx(fits.nsp, paste0("Модели_все_", Sys.Date(), ".xlsx"))

# N.species models export ---------------------------------------------------
fits.nsp.sel <- fits.nsp %>% 
    filter(
        ff %in% c("nsp ~ year + tur + site", "nsp ~ year + site"), 
        super_imp == "S.I.", 
        distr == "gaussian"
    ) %>% 
    pull(fit) %>% 
    map(~.x %>% 
            summary %>% 
            `$`(coefficients) %>% 
            as.data.frame %>% 
            rownames_to_column("predictor") %>% 
            mutate(Коэффициент = round(Estimate, 2),
                   Ошибка = round(`Std. Error`, 2),
                   t = round(`t value`, 1), 
                   p.value = round(`Pr(>|t|)`, 3),
                   p.value = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ as.character(p.value)),
                   .keep = "unused"
            )
    )

div2 %>% 
    select(year, tur, zone, site, nsp) %>% 
    mutate(
        km = case_when(substr(site, 4, 4) == "N" ~ substr(site, 2, 3), 
                       TRUE ~ paste0("-", substr(site, 2, 3))),
        km = as.numeric(km),
        tur = case_when(tur == 1 ~ "1. Начало лета", TRUE ~ "2. Конец лета"),
        # km = as.numeric(substr(site, 2, 3)), 
        type = paste0(tur, ", ", year), 
        .keep = "unused") %>% 
    group_by(type, km, zone) %>% 
    summarise(MED = median(nsp), 
              MAX = max(nsp),  # MEAN = mean(abu),
              MIN = min(nsp), # SD = sd(abu), 
              Q25 = quantile(nsp, 0.25), 
              Q75 = quantile(nsp, 0.75), 
              .groups = "drop") %>% 
    ggplot(aes(x = km, y = MED, 
               color = km, )) + 
    geom_vline(xintercept = 0, color = "coral", 
               linewidth = 2, alpha = 0.5, # linetype = "dashed"
    ) + 
    geom_linerange(aes(ymin = MIN, ymax = MAX), linewidth = 0.7) + 
    geom_crossbar(aes(ymin = Q25, ymax = Q75), fill = "white", linewidth = 0.7) +
    facet_wrap(~type, scales = "free", ncol = 2) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Количество видов (N)", 
         subtitle = "Количество видов (N) в разные годы и туры") + 
    # scale_colour_gradient2(low = "red", mid = "darkgoldenrod4", midpoint = 15, high = "green")+
    scale_colour_gradient2(low = "green", mid = "red", midpoint = 0, high = "green")+
    guides(color="none") + 
    theme(strip.background=element_rect(fill= "white"))

ggsave("diversity_nsp.png", width = 10, height = 7, dpi = 900)

# Shannon models selection ---------------------------------------------------
fits.nsp <- expand_grid(
    yy = "nsp ",
    xx = paste0("year"), 
    zz = c(" + site", " + km", " + zone"),
    super_imp = c("S.I.", ""),
    distr = c("poisson", "quasipoisson",  "gaussian"),
    tt = c("", "-t", " * tur", " + tur"),
) %>% 
    filter(super_imp == "S.I.") %>% 
    mutate(
        lg = case_when(distr != "gaussian" ~ FALSE, TRUE ~ TRUE),
        xx = paste0(xx, tt, zz), #.keep = "unused",
        turs.united = str_detect(tt, "-t")
    ) %>% 
    unite(ff, yy, xx, sep = "~ ") %>%  
    transmute(
        ff, super_imp, distr, lg, turs.united,
        nm = paste0(ff, "_", super_imp, "_", distr),
    )

fits.nsp <- fits.nsp %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){
        if(a$turs.united){
            df <- div1
            a$ff <- str_replace_all(a$ff, "-t", "")
        } else {
            df <- div2
        }
        if(a$super_imp != "S.I."){ 
            df <- filter(df, zone != "суперимпактная")
        }
        if(a$lg){df <- mutate(df, nsp = log(nsp+1))}
        glm(formula = a$ff, data = df, family = a$distr)
    }) %>% 
    `names<-`(fits.nsp$nm) %>% 
    mutate(fits.nsp, fit = ., info = lapply(., summary)) 

fits.nsp <- fits.nsp %>% 
    pull(fit) %>% 
    lapply(function(a){ 
        resid = residuals(a)
        ftted = fitted.values(a)
        if(a$family$family == "poisson"){
            suppressMessages(overdispersion <- a |> 
                                 performance::check_overdispersion() |>
                                 capture.output() |> 
                                 str_subset("value") |> 
                                 str_replace_all("[:alpha:]|[:blank:]|[:symbol:]|-", "") |> 
                                 as.numeric())
        } else {overdispersion <- NA}
        
        tibble(aic = AIC(a), 
               aic2 = manual.aic(a),
               shapiro = shapiro.test(resid)$p.value,
               overdispersion = overdispersion,
               var.res = var(resid), 
               var.ftd = var(ftted), 
               rows = nrow(a$data),
               coef = a %>% 
                   coefficients %>% 
                   data.frame(val = .) %>% 
                   rownames_to_column("var") %>% 
                   filter(str_detect(var, "year|tur")) %>% 
                   mutate(val = round(val, 2)) %>% 
                   unite("coef", sep = " = ") %>% 
                   pull(coef) %>% 
                   paste0(collapse = "; ")
        ) %>% 
            mutate(r2 = var.ftd/(var.res + var.ftd)) %>% 
            return()
    }) %>% 
    map_df(rbind, .id = "nm") %>% 
    left_join(fits.nsp, by = "nm") %>% 
    mutate(aic = case_when(is.na(aic) ~ aic2, TRUE ~ aic), .keep = "unused")

# writexl::write_xlsx(fits.nsp, paste0("Модели_все_", Sys.Date(), ".xlsx"))

# кусочно-линейная регрессия

lm_shan <- lm(shan ~ km, data = div1)
summary(lm_shan)
x2_shan <- lm(shan ~ km + km2, data = mutate(div1, km2 = km^2))
summary(x2_shan)
seg_shan <- segmented::segmented(lm_shan, seg.Z = ~km)
summary(seg_shan)


ggplot(div1, aes(x = km, y = shan)) + 
    geom_vline(xintercept = seg_shan$psi[2], color = "green", 
               linetype = "dashed", linewidth = 1)+
    # geom_smooth(formula = 'y ~ x', method = "lm", linewidth = 2) + 
    # geom_abline(slope = coefficients(lm_shan)[2], intercept = coefficients(lm_shan)[1], 
                # color = "grey", linewidth = 1.5) + 
    geom_line(data = data.frame(km = 0:33) %>% 
                  mutate(shan = predict(lm_shan, .)), 
              color = "grey", linewidth = 1.5) +
    geom_line(data = data.frame(km = 0:33, km2 = (0:33)^2) %>% 
                  mutate(shan = predict(x2_shan, .)), 
              color = "lightblue", linewidth = 1.5) + 
    geom_line(data = data.frame(km = c(0, seg_shan$psi[2], 33)) %>% 
                  mutate(shan = predict(seg_shan, .)), 
              color = "green", linewidth = 1.3) + 
    geom_point() + 
    geom_label(x = 20, y = 0.5, label = paste0("Models parameters:\n", lb)) + 
    labs(x = "Удаление от завода, км", y = "индекс Шеннона")

lb <- data.frame(
        mod = c("linear", "polynom", "segment"),
        r2  = round(map_dbl(lst(lm_shan, x2_shan, seg_shan), ~summary(.x)$adj.r.squared), 2),
        aic = round(sapply(lst(lm_shan, x2_shan, seg_shan), AIC), 1)
) %>% 
    transmute(i = paste0(mod, ": ", "r2=", r2, "; aic=", aic)) %>% 
    pull(i) %>% 
    paste0(collapse = "\n")


# Shannon models export ---------------------------------------------------
fits.nsp.sel <- fits.nsp %>% 
    filter(
        ff %in% c("nsp ~ year + tur + site", "nsp ~ year + site"), 
        super_imp == "S.I.", 
        distr == "gaussian"
    ) %>% 
    pull(fit) %>% 
    map(~.x %>% 
            summary %>% 
            `$`(coefficients) %>% 
            as.data.frame %>% 
            rownames_to_column("predictor") %>% 
            mutate(Коэффициент = round(Estimate, 2),
                   Ошибка = round(`Std. Error`, 2),
                   t = round(`t value`, 1), 
                   p.value = round(`Pr(>|t|)`, 3),
                   p.value = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ as.character(p.value)),
                   .keep = "unused"
            )
    )

div2 %>% 
    select(year, tur, zone, site, nsp) %>% 
    mutate(
        km = case_when(substr(site, 4, 4) == "N" ~ substr(site, 2, 3), 
                       TRUE ~ paste0("-", substr(site, 2, 3))),
        km = as.numeric(km),
        tur = case_when(tur == 1 ~ "1. Начало лета", TRUE ~ "2. Конец лета"),
        # km = as.numeric(substr(site, 2, 3)), 
        type = paste0(tur, ", ", year), 
        .keep = "unused") %>% 
    group_by(type, km, zone) %>% 
    summarise(MED = median(nsp), 
              MAX = max(nsp),  # MEAN = mean(abu),
              MIN = min(nsp), # SD = sd(abu), 
              Q25 = quantile(nsp, 0.25), 
              Q75 = quantile(nsp, 0.75), 
              .groups = "drop") %>% 
    ggplot(aes(x = km, y = MED, 
               color = km, )) + 
    geom_vline(xintercept = 0, color = "coral", 
               linewidth = 2, alpha = 0.5, # linetype = "dashed"
    ) + 
    geom_linerange(aes(ymin = MIN, ymax = MAX), linewidth = 0.7) + 
    geom_crossbar(aes(ymin = Q25, ymax = Q75), fill = "white", linewidth = 0.7) +
    facet_wrap(~type, scales = "free", ncol = 2) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Количество видов (N)", 
         subtitle = "Количество видов (N) в разные годы и туры") + 
    # scale_colour_gradient2(low = "red", mid = "darkgoldenrod4", midpoint = 15, high = "green")+
    scale_colour_gradient2(low = "green", mid = "red", midpoint = 0, high = "green")+
    guides(color="none") + 
    theme(strip.background=element_rect(fill= "white"))

ggsave("diversity_nsp.png", width = 10, height = 7, dpi = 900)



# Ordination --------------------------------------------------------------
dis <- wide1 %>% 
    mutate(zone = fct_collapse(
        zone, `импактная` = c("импактная", "суперимпактная"))) %>%
    select(-site) %>% 
    # group_by(year, zone, km, plot) %>% 
    # summarise_all(mean) %>% 
    # ungroup %>% 
    # # select_if(~ !is.numeric(.) || sum(.) != 0)
    unite("ID", year, zone, km, plot, sep = "_") %>% 
    filter(rowSums(.[,2:ncol(.)]) > 0 
           # str_detect(ID, "супер", negate = TRUE)
    ) %>%
    # filter(str_detect(ID, "супер", negate = TRUE))
    column_to_rownames("ID") %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
pc <- ape::pcoa(dis)
eig <- pc$values$Eigenvalues
eig <- round(eig/sum(eig)*100, 1)
pc <- pc$vectors %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    as_tibble()

pc %>% 
    # select(1:3) %>% 
    separate(ID, into = c("year", "zone", "site", "plot"), 
             sep = "_") %>% 
    mutate(#zone = str_replace_all(zone, "ая", "ый"),
           zone = factor(zone, levels = c("фоновая", "буферная", "импактная")),
           zone = fct_relabel(zone, ~paste0(.x, " территория"))) %>% 
ggplot(aes(x = Axis.1, y = Axis.2, linetype = year,
               fill = zone, color = zone, shape = year)) +
    geom_point(color = "black", size = 2.5) +  
    stat_ellipse() + 
        # mapping = aes(color = year), linewidth = 0.5, linetype = "dashed") +
    # facet_grid(cols = vars(zone)) +
    scale_shape_manual(values = c(21, 22))+ 
    scale_linetype_manual(values = c("dashed", "solid")) + 
    scale_fill_manual(values = c("red", "orange", "green")) +
    scale_color_manual(values = c("red", "orange", "green")) +
    labs(subtitle = "Ординация по динамической плотности\nТуры объединены", 
         x = paste0("Ось 1 (", eig[1], " %)"), 
         y = paste0("Ось 2 (", eig[2], " %)"), 
         # fill = "Год", color = "Год", shape = "Год"
         ) +
    ggplot2::guides(size = "none", color = "none") + 
    theme(panel.grid = element_blank())
ggsave("Plot5.svg", width = 22, height = 9, units = "cm")

m <- as.matrix(dis)
m[upper.tri(m)] <- NA
diag(m) <- NA
m2 <- m %>% 
    as.data.frame() %>% 
    rownames_to_column("id1") %>% 
    as_tibble() %>% 
    pivot_longer(names_to = "id2", values_to = "dis", -id1) %>% 
    filter(!is.na(dis)) %>% 
    separate(id1, c("year1", "zone1"), sep = "_", extra = "drop") %>% 
    separate(id2, c("year2", "zone2"), sep = "_", extra = "drop")

# zone + year
res <- list()
res$within_raw <- m2 %>% 
    filter(year1 == year2, zone1 == zone2) %>% 
    transmute(i = paste0(substr(zone1, 1,3), "_", year1), dis) %>% 
    group_by(i) %>% 
    summarise(dis = mean(dis), .groups = "drop") %>% 
    mutate(i = factor(i, levels = c("фон_2009", "фон_2014", "буф_2009", "буф_2014", "имп_2009", "имп_2014" ))) %>% 
    arrange(i)

m2 %>% 
    filter(year1 != year2 | zone1 != zone2) %>% 
    group_by(year1, zone1, year2, zone2) %>% 
    summarise(dis = mean(dis), .groups = "drop") %>% 
    mutate(id1 = paste0(year1, "_", zone1), 
           id2 = paste0(year2, "_", zone2), 
           .keep = "unused") %>% 
    mutate_at(2:3, function(a){substr(a, 1, 8)}) %>% 
    pivot_wider(names_from = id2, values_from = dis)

# zone
m2 %>% 
    filter(zone1 == zone2) %>% 
    group_by(zone1) %>% 
    summarise(dis = mean(dis), .groups = "drop")

m2 %>% 
    filter(zone1 != zone2) %>% 
    group_by(zone1, zone2) %>% 
    summarise(dis = mean(dis), .groups = "drop")


# dis other way -----------------------------------------------------------
res$within_2.axes <- pc %>% 
    group_by(year, zone) %>% 
    mutate(Axis.1m = mean(Axis.1), Axis.2m = mean(Axis.2)) %>% 
    mutate(zone = substr(zone, 1, 3), 
           d = sqrt((Axis.1-Axis.1m)^2 + (Axis.2-Axis.2m)^2)) %>% 
    summarise(d = mean(d), .groups = "drop") %>% 
    unite("i", zone, year, sep = "_") %>% 
    mutate(i = factor(i, levels = c("фон_2009", "фон_2014", "буф_2009", "буф_2014", "имп_2009", "имп_2014" ))) %>% 
    arrange(i)

m3 <- pc %>% 
    column_to_rownames("ID") %>% 
    dist() %>% 
    as.matrix()

m3[upper.tri(m3)] <- NA
diag(m3) <- NA
m3 <- m3 %>% 
    as.data.frame() %>% 
    rownames_to_column("id1") %>% 
    as_tibble() %>% 
    pivot_longer(names_to = "id2", values_to = "dis", -id1) %>% 
    filter(!is.na(dis)) %>% 
    separate(id1, c("year1", "zone1"), sep = "_", extra = "drop") %>% 
    separate(id2, c("year2", "zone2"), sep = "_", extra = "drop")




    mutate(ID = paste0(substr(ID, 6, 9)), substr(ID, 6, 9)), 

pc

    ungroup() %>% 
    mutate(i = paste0(substr(zone, 1,3), "_", year), .keep = "unused") %>% 
    column_to_rownames("i") %>% 
    dist()



# export ------------------------------------------------------------------
list(abundance = fits.abu, 
     number.of.species = fits.nsp, 
     shannon.index = fits.shn) %>% 
    lapply(function(a){select(a, -fit, -info)}) %>% 
    writexl::write_xlsx(paste0("diversity.models_", Sys.Date(), ".xlsx"))


# Correlation -------------------------------------------------------------
corr <- function(df, x, y, grp){
    # x = "abu"; # y = "shan"; # grp = c("year", "zone")
    # unify arguments 
    if(is.numeric(x)){x <- colnames(df)[x]}
    if(is.numeric(y)){x <- colnames(df)[y]}
    if(is.numeric(grp)){grp <- colnames(df)[grp]}
    # run
    df %>% 
        unite(all_of(grp), col = type, sep = "_") %>% 
        select(type, all_of(c(x, y))) %>% 
        rename(x = 2, y = 3) %>% 
        split(.$type) %>% 
        map(~cor.test(.x$x, .x$y, method = "spearman")) %>% 
        suppressWarnings() %>% 
        suppressMessages() %>% 
        map_dfr(~tibble(r = .x$estimate, p = .x$p.value), .id = "type") %>% 
        separate(type, sep = "_", into = grp)
}

corr_result <- rbind(mutate(corr(div2, "km", "abu", c("year", "tur")), 
             year = paste0(year, "_t", tur), .keep = "unused"), 
      corr(div1, "km", "nsp", c("year")), 
      corr(div1, "km", "shan", c("year"))) %>% 
    mutate(x = c(rep("abu", 4), rep("nsp", 2), rep("shan", 2)), 
           y = "km", 
           date = year, 
           p.val = as.character(round(p, 3)),
           p.val = case_when(p.val == "0" ~ "<0.001", 
                             TRUE ~ paste0(" =", p.val)),
           lab = paste0("ρ = ", round(r, 2), "\np.value ", p.val)
    )
corr_result %>% 
    transmute(x, y, year, r = round(r, 2), p = round(p, 6)) %>% 
    writexl::write_xlsx(paste0("correlation_", Sys.Date(), ".xlsx"))

div2 %>% 
    group_by(year, tur, km) %>% 
    summarise(abu = mean(abu), 
              .groups = "drop_last") %>% 
    mutate(km2 = rank(km)) %>% 
    ungroup() %>% 
    mutate(date = paste0(year,"_t", tur), .keep = "unused") %>% 
    ggplot(aes(x = km2, y = abu, color = date)) + 
    geom_smooth(method = "lm", formula = "y~x", se = FALSE) + 
    geom_point(size = 2.7, color = "white") +
    geom_point(size = 2) +
    geom_label(aes(x = 3, y = 550, label = lab), 
               data = filter(corr_result, x == "abu"),
               alpha = 0.6, size = 3) +
    facet_wrap(~date) +
    scale_x_continuous(breaks = 1:10,
                       labels = c(1, 4, 5, 9, 11, 12, 18, 26, 27, 32)) +
    labs(x = "Distance from smelter (not to scale)", 
         subtitle = "all plots are united",
         y = "Abundance") +
    guides(color=guide_legend(title="years & turs")) +
    theme(
        strip.background = element_rect(fill = "ghostwhite"),
        panel.grid = element_blank(),
        # panel.grid.minor = element_blank(),
        legend.position = "none")
# ggsave("corr_abu.png", height = 6, width = 7, dpi = 900)

div1 %>% 
    group_by(year, km) %>%
    summarise(nsp = mean(nsp),
              shan= mean(shan),
              .groups = "drop_last") %>%
    mutate(km2 = rank(km)) %>% 
    ungroup() %>% 
    select(year, km2, nsp, shan) %>% 
    pivot_longer(names_to = "x", values_to = "val", -1:-2) %>% 
    mutate(type = paste0(year,x), x = case_when(
        x == "nsp" ~ "Number of species", 
        TRUE ~ "Shannon diversity")) %>% 
    ggplot(aes(x = km2, y = val, color = paste0(year, x))) + # , color = date
    geom_smooth(method = "lm", formula = "y~x", se = FALSE) +
    geom_point(size = 2.7, color = "white") +
    geom_point(size = 2) +
    geom_label(aes(x = 8, y = c(3, 3, 0.3, 0.3), label = lab),
               data = corr_result |> 
                    filter( x != "abu") |> 
                    mutate(x = case_when(
                       x == "nsp" ~ "Number of species", 
                       TRUE ~ "Shannon diversity")),
               alpha = 0.6, size = 3) +
    facet_grid(cols = vars(year), rows = vars(x),
               scales = "free_y") +
    scale_x_continuous(breaks = 1:10,
                       labels = c(1, 4, 5, 9, 11, 12, 18, 26, 27, 32)) +
    labs(x = "Distance from smelter (not to scale)", 
         subtitle = "all plots are united",
         y = NULL) +
    guides(color=guide_legend(title="years & turs")) +
    theme(
        strip.background = element_rect(fill = "ghostwhite"),
        panel.grid = element_blank(), 
        # panel.grid.minor = element_blank(), 
        legend.position = "none")
# ggsave("corr_nsp.shan.png", height = 6, width = 7, dpi = 900)








