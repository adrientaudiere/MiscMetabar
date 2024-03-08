## Benjamin Altmann comments

1. Examples for unexported function
  plot_guild_pq() in:
     merge_samples2.Rd
Please export those functions or omit the example.

2. You still have some \dontrun{} examples. Please change that to \donttest{} or let us know why \dontrun{} is necessary in your case.

3. Note, please wrap examples that need packages in 'Suggests' in if(requireNamespace("pkgname")){} instead.

## Answers

1. Sorry, I don't understand this error. merge_samples2() and plot_guild_pq() are exported. Moreover, there is no reference to plot_guild_pq() in the file merge_samples2.Rd.

2. There are 5 still \dontrun in examples for functions relying on external software (namely blast, cutadapt, mumu and krona x2).

3. Done. There's no more library() calls in examples. Thanks for the comment.