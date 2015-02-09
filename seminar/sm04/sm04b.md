# sm04b
Eva Y  
January 28, 2015  

**Take-home Exercise: Can you modify the code above to answer this questions: which country gains the most growth in GDP in a 5-year interval? Between which 2 years?**


```r
suppressPackageStartupMessages(library(dplyr))

# load gapminder data
gd_url <- "http://tiny.cc/gapminder"
gtbl <- gd_url %>% read.delim %>% tbl_df
gtbl %>% glimpse
```

```
## Variables:
## $ country   (fctr) Afghanistan, Afghanistan, Afghanistan, Afghanistan,...
## $ year      (int) 1952, 1957, 1962, 1967, 1972, 1977, 1982, 1987, 1992...
## $ pop       (dbl) 8425333, 9240934, 10267083, 11537966, 13079460, 1488...
## $ continent (fctr) Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asi...
## $ lifeExp   (dbl) 28.801, 30.332, 31.997, 34.020, 36.088, 38.438, 39.8...
## $ gdpPercap (dbl) 779.4453, 820.8530, 853.1007, 836.1971, 739.9811, 78...
```

```r
# country that gained the most growth in GDP in a 5-year interval
gtbl %>%
  group_by(country) %>%
  select(country, year, continent, gdpPercap) %>%
  mutate(le_delta = gdpPercap - lag(gdpPercap)) %>%
  summarize(highest_le_delta = max(le_delta, na.rm = TRUE)) %>%
  filter(min_rank(desc(highest_le_delta)) < 2)
```

```
## Source: local data frame [1 x 2]
## 
##   country highest_le_delta
## 1  Kuwait         28452.98
```

```r
# between which 2 years?
gtbl %>%
  filter(country == "Kuwait") %>%
  select(country, year, continent, gdpPercap) %>%
  mutate(le_delta = gdpPercap - lag(gdpPercap)) %>%
  filter(min_rank(desc(le_delta)) < 2)
```

```
## Source: local data frame [1 x 5]
## 
##   country year continent gdpPercap le_delta
## 1  Kuwait 1972      Asia  109347.9 28452.98
```

**Kuwait is the country that gained the most growth in GDP per capita ($28452.98) in a 5-year interval (between 1967-1972).**
