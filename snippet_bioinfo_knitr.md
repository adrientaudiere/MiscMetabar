---
title: "snippet_bioinfo_knitr"
output: html_document
---

```{r knitr set up, include=FALSE}
library(knitr)
hooks <- knitr::knit_hooks$get()
hook_foldable <- function(type, background_color) {
  force(type)
  function(x, options) {
    res <- hooks[[type]](x, options)

    if (isFALSE(options[[paste0("fold.", type)]])) {
      return(res)
    }

    paste0(
      "<details><summary style ='background-color: ",
      background_color,
      "'>",
      type,
      "</summary>\n\n",
      res,
      "\n\n</details>"
    )
  }
}
knitr::knit_hooks$set(
  output = hook_foldable("output", background_color = "#a6d3e3"),
  message = hook_foldable("message", background_color = "#fae2df "),
  warning = hook_foldable("warning", background_color = "#dc863b")
)

opts_chunk$set(
  fig.align = "center",
  fig.retina = 2,
  fig.width = 10
)
```

## Plots

```{r Visualise target plan, out.width='80%'}
# targets::tar_visnetwork(targets_only = T)
```

## Tables

```{r results='asis'}
gt::gt(cars)
```

```{r results='asis'}
#res <- adonis_phyloseq(
#  physeq = d_asv_ecm,
#  formula = "Type_echantillon + Type_sol * Sp"
#)
#gt(round(res$aov.tab, 3), rownames_to_stub = TRUE,
#   caption="Résultat de l\'analyse de variance sur les distances (Permanova)")
```


## Diagrammes

```{r}
DiagrammeR::mermaid("
graph LR
  subgraph Résilience
    Che[Chênaie]---Sol
    Che---Apex
    Apex-.-> Qi[Quercus ilex] 
    Apex-.-> Cm[Cistus monspeliensis]
    Apex-.-> Au[Arbutus unedo]  
    Apex-.-> Ph[Pinus halepensis]
    Urb[Urbain]---Sol
    Urb---Apex
    Gar[Garrigues]---Sol
    Gar---Apex

  end
")
```

```{r results="asis"}
DiagrammeR::grViz("
digraph a_nice_graph {

graph [layout = dot,
rankdir = LR]

# node definitions with substituted label text
node [fontname = Helvetica, shape=box, style = filled]

a [label = '@@1']
b [label = '@@2']
c [label = 'The end']

# edge definitions with the node IDs
a -> {c} [label = '+ yOp']
c -> {a} [label = '+ yEp']
b -> {c} [label = '- 20']
}

[1]: paste0('SPEED \\n (',mean(cars$speed), ' +/-, ', sd(cars$speed), ' s.)')
[2]: paste0('DIST \\n (',mean(cars$dist), ' +/-, ', sd(cars$speed), ' s.)')
")
```

## CSS class

::: in_a_nutshell
CSS class in a nutshell for synthesis
:::

::: question
CSS class question for ... question
:::

### Collapsing a part of the results

<button class="btn btn-primary" data-toggle="collapse" data-target="#details">
See the details
</button>  

::: {#details .collapse}
Parts to collapse

```r
summary(cars)
```
:::



## sub-naviguation

## Level ONE   {.tabset .tabset-fade .tabset-pills}
### sublevel A  {.tabset .tabset-fade}

#### SPEED - DIST
```{r}
plot(cars$speed~ cars$dist)
```

#### Linear model
```{r}
plot(cars$speed~ cars$dist)
abline(lm(cars$speed~ cars$dist))
```

### sublevel B  {.tabset .tabset-fade}

#### SPEED - DIST
```{r}
plot(cars$speed*60~ sqrt(cars$dist))
```

#### Linear model
```{r}
plot(cars$speed*60~ sqrt(cars$dist))
abline(lm(cars$speed*60~ sqrt(cars$dist)))
```
