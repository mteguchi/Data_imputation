
# functions and constants used in making the presentation.


color.JM <- "orangered"
fill.color.JM <- "lightpink"
color.W <- "olivedrab4"
fill.color.W <- "mediumseagreen"

impute.color <- "midnightblue"

period.JM <- 12
period.W <- 6
maxN <- 10000

data.extract <- function(location, 
                         year.begin, 
                         year.end, 
                         season.begin = year.begin, 
                         season.end = year.end){
  # In March 2019, we received new data for 2018. So, the raw data file
  # has been updated.  
  # On 16 April 2019, the last few data points for 2019 were received
  # so the data files have been updated. 
  if (is.null(season.begin)) season.begin <- year.begin
  if (is.null(season.end)) season.end <- year.end
  
  if (location == "JM"){
    data.0 <- read.csv("data/JM_nests_April2019.csv")
    
    data.0 %>% 
      select(Year_begin, Month_begin, JM_Nests) %>%
      mutate(Nests = JM_Nests) -> data.0
    
  } else if (location == "W"){
    data.0 <- read.csv("data/W_nests_April2019.csv")
    data.0 %>% 
      select(Year_begin, Month_begin, W_Nests) %>%
      mutate(Nests = W_Nests) -> data.0
  }
  
  # create regularly spaced time series:
  data.2 <- data.frame(Year = rep(min(data.0$Year_begin,
                                      na.rm = T):max(data.0$Year_begin,
                                                     na.rm = T),
                                  each = 12),
                       Month_begin = rep(1:12,
                                         max(data.0$Year_begin,
                                             na.rm = T) -
                                           min(data.0$Year_begin,
                                               na.rm = T) + 1)) %>%
    mutate(begin_date = as.Date(paste(Year,
                                      Month_begin,
                                      '01', sep = "-"),
                                format = "%Y-%m-%d"),
           Frac.Year = Year + (Month_begin-0.5)/12) %>%
    select(Year, Month_begin, begin_date, Frac.Year)
  
  # also make "nesting season" that starts April and ends March
  
  data.0 %>% mutate(begin_date = as.Date(paste(Year_begin,
                                               Month_begin,
                                               '01', sep = "-"),
                                         format = "%Y-%m-%d")) %>%
    mutate(Year = Year_begin,
           Month = Month_begin,
           f_month = as.factor(Month),
           f_year = as.factor(Year),
           Frac.Year = Year + (Month_begin-0.5)/12) %>%
    select(Year, Month, Frac.Year, begin_date, Nests) %>%
    na.omit() %>%
    right_join(.,data.2, by = "begin_date") %>%
    transmute(Year = Year.y,
              Month = Month_begin,
              Frac.Year = Frac.Year.y,
              Nests = Nests,
              Season = ifelse(Month < 4, Year-1, Year),
              Seq.Month = ifelse(Month < 4, Month + 9, Month - 3)) %>%
    reshape::sort_df(.,vars = "Frac.Year") %>%
    filter(Season >= season.begin & Season <= season.end) -> data.1
  
  data.1 %>% filter(Month > 3 & Month < 10) -> data.summer
  data.1 %>% filter(Month > 9 | Month < 4) %>%
    mutate(Seq.Month = Seq.Month - 6) -> data.winter
  
  jags.data <- list(y = log(data.1$Nests),
                    m = data.1$Seq.Month,
                    T = nrow(data.1))
  
  y <- matrix(log(data.1$Nests), ncol = 12, byrow = TRUE)
  
  jags.data2 <- list(y = y,
                     m = matrix(data.1$Seq.Month, 
                                ncol = 12, byrow = TRUE),
                     n.years = nrow(y))
  
  y.summer <- matrix(log(data.summer$Nests),
                     ncol = 6, byrow = TRUE)
  
  y.winter <- matrix(log(data.winter$Nests),
                     ncol = 6, byrow = TRUE)
  
  jags.data2.summer <- list(y = y.summer,
                            m = matrix(data.summer$Seq.Month, 
                                       ncol = 6, byrow = TRUE),
                            n.years = nrow(y.summer))
  
  jags.data2.winter <- list(y = y.winter,
                            m = matrix(data.winter$Seq.Month, 
                                       ncol = 6, byrow = TRUE),
                            n.years = nrow(y.winter))
  
  out <- list(jags.data = jags.data,
              jags.data2 = jags.data2,
              jags.data.summer = jags.data2.summer,
              jags.data.winter = jags.data2.winter,
              data.1 = data.1,
              data.summer = data.summer,
              data.winter = data.winter)
  return(out)
}



Girondot_fcn <- function(d, S, K, P, min, max){
  K <- abs(K)
  S <- abs(S)
  S1 <- -S
  M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
  M2 <- (1 + (2 * exp(K) - 1) * exp((1/S) * (P - d))) ^ (-1/exp(K))
  N <- min + (max - min) * (M1 * M2)
  return(N)
}

# Extracting posterior samples of deviance or any other variable from jags output:
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}
