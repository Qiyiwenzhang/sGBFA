args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (length(args)) as.integer(args[1]) else 1L
stopifnot(!is.na(rep_id) && rep_id >= 1L)
idx = rep_id


source("/ihome/qzhang/qiz201/aim_1c/aim_1c_final.R")

###### What's followed is to generate Gaussian data #####
outcome.no  = c("outcome_1","outcome_5","outcome_10")
noise.level = c("high","med","low")
outcome.sel = outcome.no[3]
noise.sel = noise.level[3]

base_in_data  <- paste0("/ihome/qzhang/qiz201/aim_1c/simulation_data/outcome_sim/", outcome.sel)
in_file_data  <- file.path(
  base_in_data,
  sprintf("save_data_gaussian_%s_ai_outcome.rds",noise.sel)
)

save_data_mask =readRDS(in_file_data)

base_out_result =  "/ix1/qzhang/qiz201/aim_1c/result_outcome"
out_file_result  <- file.path(
  base_out_result,
  outcome.sel,
  sprintf("result_with_outcome_%s",noise.sel)
)

path = out_file_result

T = 2500
start = 500
end = T
g=2
s=1
outcome = 1 
system.time({
  sim_list = sGBFA_sim(idx=idx,outcome=outcome, T=T, start=start, end=end, g=g, path=path)
})


