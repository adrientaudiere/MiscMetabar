# Session Context

## User Prompts

### Prompt 1

in function tax_bar_pq add a new parameter to set the order of the modality. The idea is to be able to change the order of each bar using e.g.
physeq@sam_data[["modality"]] <- forcats::fct_relevel(physeq@sam_data[["modality"]], c(new_param))
We must help the user if the value in the new_param don't overlap the value in physeq@sam_data[["modality"]], if there is less value in new_param it's ok the user want to see only some modalities (just cli_message), if there is value in new_param not pres...

