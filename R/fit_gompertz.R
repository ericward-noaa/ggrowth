#' fit_gompertz is the primary function for fitting growth data
#' @param dataframe The dataframe with columns: Year,
#'
#' @return list
#' @export
#'
fit_gompertz = function(dataframe) {
  dataframe = dplyr::arrange(dataframe, Year, Month, Day)
  dataframe = dplyr::group_by(dataframe, ID) %>%
    mutate(refMonth = Month[1], refYear = Year[1], refDay = Day[1]) %>%
    mutate(Jday = date::mdy.date(Month,Day,Year)-date::mdy.date(refMonth,refDay,refYear)) %>%
    filter(-refMonth, -refYear, -refDay)

  # taken from stackoverflow: http://stackoverflow.com/questions/35996877/fitting-multiple-nls-functions-with-dplyr
  starts = data.frame(ID = unique(dataframe$ID),
    Sinf_start = rnorm(length(unique(dataframe$ID)), 1000, sd = 10),
    k_start = rnorm(length(unique(dataframe$ID)), 0.1, sd = 0.03),
    t0_start = rnorm(length(unique(dataframe$ID)), -1, sd = 0.1))

  mods = dataframe %>%
    left_join(starts,by="ID") %>%
    group_by(ID) %>%
    do(tidy(nls(Weight ~ Sinf * exp(-exp(-K * (Jday - t0))),
      data = .,
      start = list(K = first(.$k_start),
        Sinf = first(.$Sinf_start),
        t0 = first(.$t0_start)))))

  reshape_pars <- reshape2::melt(mods[,c("ID","term","estimate")], id.vars = c("ID", "term"))
  pars = reshape2::dcast(reshape_pars, ID ~ term + variable)

  df_w_pars = left_join(dataframe, pars,by="ID")
  df_w_pars$predict = df_w_pars$Sinf_estimate * exp(-exp(-df_w_pars$K_estimate * (df_w_pars$Jday - df_w_pars$t0_estimate)))

  return(select(df_w_pars, -Sinf_estimate, -K_estimate, -t0_estimate))
}
