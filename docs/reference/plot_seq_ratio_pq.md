# A diagnostic plot of the number of sequences per samples

A diagnostic plot of the number of sequences per samples

## Usage

``` r
plot_seq_ratio_pq(physeq, min_nb_seq = 1000, annotations = TRUE)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- min_nb_seq:

  (int) The minimum number of sequences per samples to compare the
  ratio.

- annotations:

  (logical, default TRUE). If FALSE, no annotations are plotted

## Value

A ggplot2 object

## Details

The x axis depict the number of sequences per samples and the y axis
depicted the ratio of the number of sequences for a given sample divide
by the number of sequences of the previous sample when ordered by the
number of sequences. A high ratio indicate an important and quick
increase of the number of sequence which may indicate that below this
ratio, samples are suspicious.

The general idea is to first removed all samples with definitively not
enough sequences and then, among the kept samples, find the higher
augmentation (ratio) to possibly detect suspicious samples.

## Author

Adrien Taudière

## Examples

``` r
plot_seq_ratio_pq(data_fungi, min_nb_seq = 200)

data(GlobalPatterns)
plot_seq_ratio_pq(GlobalPatterns, min_nb_seq = 100000)

plot_seq_ratio_pq(data_fungi_mini, min_nb_seq = 10, annotations = FALSE)
```
