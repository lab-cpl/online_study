pacman::p_load(
    tidyverse,
    ggplot2,
    ggcorrplot
)

# https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
# get path of source file
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file)==0)
    {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}
# sets path on source file location
script_path <- getCurrentFileLocation()
setwd(script_path)

# load csv
files <-
    list.files(
        pattern = "[0-9]+_[a-z]+.csv",
        full.names = TRUE
    )

# load metadata
metadata <-
    read_csv("final_dataset_filtered.csv") %>% 
    group_by(ID) %>% 
    slice_head(n=1) %>% 
    select(
        ID,
        Sex
    ) %>% 
    mutate(
        ID = as.factor(ID)
    )

data <-
    files %>% 
    map(
        ., .progress = TRUE, function(f){
            raw <- 
                read_csv(f, show_col_type = FALSE)
            ID <-
                str_extract(f, "[0-9]+")
            context <-
                str_extract(f, "[a-z]+")
            out <-
                raw %>% 
                mutate(
                    ID = ID,
                    context = context
                )
            return(out)
        }
    )

d <- bind_rows(data)

write_csv(bind_rows(data), "ddm_fits.csv")

data_stitch <-
    bind_rows(data) %>% 
    mutate(
        ID = as.factor(ID)
    ) %>% 
    left_join(
        ., metadata, by = c("ID")
    )

p1_hist <-
    data_stitch %>% 
    pivot_longer(
        cols = c(d_t:bias),
        names_to = "parameter",
        values_to = "val"
    ) %>% 
    ggplot(aes(
        val, color = Sex
    )) +
    geom_density() +
    facet_wrap(~parameter, scales = "free")
p1_hist

p2_boxplot <-
    data_stitch %>% 
    pivot_longer(
        cols = c(d_t:bias),
        names_to = "parameter",
        values_to = "val"
    ) %>% 
    ggplot(aes(
        Sex, val
    )) +
    stat_summary(fun.data = "mean_se", geom = "point", size = 2, shape = 15) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", size = 1.5, width = 0.2) +
    facet_wrap(~parameter * context, scales = "free")
p2_boxplot

data_stitch %>% 
    group_by(Sex, context) %>% 
    group_split() %>% 
    map(., function(df){
        comb <- paste(unique(df$context), unique(df$Sex), sep = "_")
        data <- df %>% select(d_t:bias)
        out <- cor(data)
        return(ggcorrplot(out, hc.order = FALSE, outline.col = "white", type = "lower") + ggtitle(comb))
    })

# single subject
s <-
    d %>% 
    filter(ID == "6835953", context == "typical")
s

tibble(
    x = 1:200,
    w_taste = s$d_t * x,
    w_health = s$d_h * x
) %>% 
    pivot_longer(
        cols = w_taste:w_health,
        names_to = "names",
        values_to = "values"
    ) %>% 
    ggplot(aes(
        x, values, color = names
    )) +
    geom_line()

