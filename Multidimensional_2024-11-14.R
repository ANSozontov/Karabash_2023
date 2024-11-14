# Data load ---------------------------------------------------------------
library(tidyverse)
theme_set(
    theme_bw() +
        theme(
            text = element_text(family = "serif", size = 15),
            legend.position = "bottom"
        )
)

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

# 1 = turs 1 and 2 are united
wide1 <- long %>% # div2 - turs are united
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    select(-no_insects)

dis <- wide1 %>% 
    mutate(zone = fct_collapse(
        zone, `импактная` = c("импактная", "суперимпактная"))) %>%
    select(-site) %>% 
    # group_by(year, zone, km, plot) %>% 
    # summarise_all(mean) %>% 
    # ungroup %>% 
    # # select_if(~ !is.numeric(.) || sum(.) != 0)
    unite("ID", year, zone, km, plot, sep = "_") %>% 
    # filter(rowSums(.[,2:ncol(.)]) > 0 
    #        # str_detect(ID, "супер", negate = TRUE)
    # ) %>%
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

# Ordination viz ----------------------------------------------------------
pc_data <- pc %>% 
    separate(ID, into = c("year", "zone", "site", "plot"), 
             sep = "_") %>% 
    mutate(#zone = str_replace_all(zone, "ая", "ый"),
           zone = factor(zone, levels = c("фоновая", "буферная", "импактная")),
           zone = fct_relabel(zone, ~paste0(.x, " территория"))) 

pc_data %>% 
ggplot(aes(x = Axis.1, y = Axis.2, linetype = year,
               fill = zone, color = zone, shape = year)) +
    geom_point(color = "black", size = 2.5) +  
    stat_ellipse() + 
    scale_shape_manual(values = c(21, 22))+ 
    scale_linetype_manual(values = c("dashed", "solid")) + 
    scale_fill_manual(values = c("green", "orange", "red")) +
    scale_color_manual(values = c("green", "orange", "red")) +
    labs(subtitle = "Ординация по динамической плотности\nТуры объединены", 
         x = paste0("Ось 1 (", eig[1], " %)"), 
         y = paste0("Ось 2 (", eig[2], " %)"), 
         # fill = "Год", color = "Год", shape = "Год"
         ) +
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
    # group_by(year1, zone1, year2, zone2) %>% 
    # summarise(dis = mean(dis), .groups = "drop") %>% 
    mutate(id1 = paste0(substr(zone1, 1, 3), "_", year1), 
           id2 = paste0(substr(zone2, 1, 3), "_", year2), 
           .keep = "unused") %>% 
    mutate(ID1 = min(id1, id2), ID2 = max(id1, id2), .keep = "unused") %>% 
    mutate_at(2:3, function(a)factor(a, levels = c(
        "фон_2009", "фон_2014", "буф_2009", "буф_2014", "имп_2009", "имп_2014"
    ))) 
    # arrange(desc(id1), id2)
    # mutate_at(2:3, function(a){substr(a, 1, 8)}) %>% 
    pivot_wider(names_from = ID2, values_from = dis, 
                values_fill = NA, values_fn = mean) 
    arrange(ID1)

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
    separate(ID, into = c("year", "zone", "site", "plot"), 
             sep = "_") %>% 
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

res$within_all.axes <- m3 %>% 
    filter(year1 == year2, zone1 == zone2) %>% 
    transmute(i = paste0(substr(zone1, 1,3), "_", year1), dis) %>% 
    group_by(i) %>% 
    summarise(dis = mean(dis), .groups = "drop") %>% 
    mutate(i = factor(i, levels = c("фон_2009", "фон_2014", "буф_2009", "буф_2014", "имп_2009", "имп_2014" ))) %>% 
    arrange(i)



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








