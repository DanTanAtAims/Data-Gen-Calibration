all_fg_data = get_survival_entries(CSV.read(CORAL_OBS_PATH, DataFrame))
all_fg_data[!, :Cluster]  .= String.(all_fg_data.Cluster)
all_fg_data[!, :Reef]     .= String.(all_fg_data.Reef)
all_fg_data[!, :Site_UID] .= String.(all_fg_data.Site_UID)

tabular_acropora_data       = all_fg_data[is_tabular_Acropora.(all_fg_data.Taxa), :]
corymbose_acropora_data     = all_fg_data[is_corymbose_Acropora.(all_fg_data.Taxa), :]
corymbose_non_acropora_data = all_fg_data[is_corymbose_non_Acropora.(all_fg_data.Taxa), :]
small_massive_data          = all_fg_data[is_small_massive.(all_fg_data.Taxa), :]
large_massive_data          = all_fg_data[is_large_massive.(all_fg_data.Taxa), :]
