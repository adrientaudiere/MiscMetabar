# rotl wrapper for phyloseq data

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Make a taxonomic tree using the ASV names of a physeq object and the
Open Tree of Life tree.

## Usage

``` r
rotl_pq(
  physeq,
  taxonomic_rank = c("Genus", "Species"),
  context_name = "All life",
  discard_genus_alone = TRUE,
  pattern_to_remove_tip = c("ott\\d+|_ott\\d+"),
  pattern_to_remove_node = c("_ott.*|mrca*")
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxonomic_rank:

  (Character) The column(s) present in the @tax_table slot of the
  phyloseq object. Can be a vector of two columns (e.g. the default
  c("Genus", "Species")). If only one column is set it need to be format
  in this way ("Genus species" for ex. "Quercus robur") with a space.

- context_name:

  : can bue used to select only a part of the Open Tree of Life. See
  `?rotl::tnrs_contexts()` for available values

- discard_genus_alone:

  (logical) If TRUE (default), genus without information at the species
  level are discarded.

- pattern_to_remove_tip:

  (character regex string) A regex to remove unwanted part of tip names.
  If set to null, tip names are left intact.

- pattern_to_remove_node:

  (character regex string) A regex to remove unwanted part of node
  names. If set to null, node names are left intact.

## Value

A plot

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to `rotl` package if you use this function.

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("rotl")) {
  tr <- rotl_pq(data_fungi_mini, pattern_to_remove_tip = NULL)
  plot(tr)

  tr_Asco <- rotl_pq(data_fungi,
    taxonomic_rank = c("Genus", "Species"),
    context_name = "Ascomycetes"
  )
  plot(tr_Asco)
}
#> Loading required namespace: rotl
#> Warning: Dropping singleton nodes with labels: mrcaott109ott9895, mrcaott109ott50573, mrcaott109ott233596, mrcaott109ott206, mrcaott109ott240843, mrcaott109ott11233, mrcaott11233ott46244, mrcaott11233ott47319, mrcaott47319ott135577, mrcaott47319ott317441, mrcaott47319ott55252, mrcaott2427ott25198, mrcaott25198ott140745, Lyophyllaceae ott1048880, mrcaott25198ott343967, mrcaott343967ott434508, Ossicaulis ott343990, mrcaott30821ott46242, Pterulaceae ott637619, Radulomyces ott82394, mrcaott471ott44719, mrcaott471ott5265, mrcaott32867ott86441, mrcaott86441ott477103, mrcaott86441ott159293, mrcaott159293ott212207, Hericium ott1044744, mrcaott106245ott181624, mrcaott106245ott295001, mrcaott106245ott424266, mrcaott106245ott237731, mrcaott1939ott4100, mrcaott1939ott86215, mrcaott1939ott52426, mrcaott1939ott7501, mrcaott1939ott459078, mrcaott1939ott31634, mrcaott1939ott228595, mrcaott1939ott137433, Trametes ott205112, mrcaott6711ott6854, mrcaott6711ott39549, Meruliaceae ott42234, mrcaott42575ott455004, Schizoporaceae ott580418, Xylodon ott640041, Peniophorella ott351728, Marchandiomyces ott985876, mrcaott2361ott40692, Cantharellales ott558119, mrcaott2361ott10072, mrcaott2361ott70669, mrcaott2361ott30094, mrcaott2361ott133843, Hydnaceae ott160852, Sistotrema ott604159, mrcaott3860ott11234, mrcaott11234ott134670, Auricularia ott183795, Exidia ott133079, Aporpiaceae ott938070, Elmerina ott654643, Fomes ott5342097

#> Warning: Stereum ostrea, Xylodon raduloides, Ossicaulis lachnopus, Stereum hirsutum, Sistotrema oblongisporum, Fomes fomentarius, Mycena renati, Radulomyces molaris, Elmerina caryae, Hyphoderma roseocremeum, Hyphoderma setigerum, Trametes versicolor, Exidia glandulosa, Peniophorella pubera, Auricularia mesenterica, Marchandiomyces buckii, Hericium coralloides, Sistotremastrum niveocremeum, Clitopilus hobsonii, Pluteus plautus, Tremella lloydiae-candidae, Calocera cornea, Burgoa verzuoliana, Pleurotus pulmonarius, Fibulobasidium inconspicuum, Sistotremastrum guttuliferum, Chondrostereum purpureum, Auricularia auricula-judae, Sistotrema brinkmannii, Omphalotus olearius, Merismodes fasciculata, Stypella subgelatinosa, Amaurodon viridis, Sistotrema raduloides, Sistotrema coronilla, Tremella mayrhoferi, Mycena abramsii, Mortierella humilis, Scopuloides hydnoides, Coprinellus micaceus, Peniophorella praetermissa, Mycenella trachyspora, Umbelopsis isabellina, Gymnopilus penetrans, Tremella encephala, Stypella grilletii, Phallus impudicus, Absidia glauca are not matched
#> Warning: Dropping singleton nodes with labels: mrcaott235ott67231, mrcaott235ott22423, mrcaott235ott123355, mrcaott235ott6657, mrcaott235ott58888, Eurotiales ott800595, mrcaott235ott26171, Aspergillus ott550772, mrcaott1724ott67014, mrcaott1724ott34914, mrcaott1724ott13442, mrcaott27392ott600504, mrcaott138369ott219756, Cyphellophoraceae ott5345109, mrcaott15397ott21687, mrcaott15397ott138376, mrcaott15397ott1086861, mrcaott15397ott681640, mrcaott15397ott193689, mrcaott15397ott818009, Phialophora ott818011, mrcaott56765ott361027, mrcaott56765ott214858, Veronaea ott393478, mrcaott274726ott1051945, mrcaott25274ott30372, Lecidella ott518627, mrcaott115563ott408216, mrcaott1320ott18865, mrcaott1320ott70963, mrcaott1320ott22976, mrcaott1320ott8438, mrcaott1320ott372181, mrcaott1320ott33787, mrcaott1320ott689966, mrcaott1320ott37385, mrcaott1320ott4382, mrcaott1320ott31674, mrcaott1320ott4380, mrcaott4387ott168922, mrcaott4387ott18228, mrcaott18228ott707011, mrcaott18228ott721356, Coniothyriaceae ott5345468, Coniothyrium ott679340, mrcaott43686ott928515, mrcaott928515ott1002127, Lophiostomataceae ott659815, Platystomum ott76428, Paraphoma ott948553, Sporormiaceae ott335006, Massarina ott36173, Suttonomyces ott7513409, Biatriosporaceae ott5345472, Biatriospora ott295651, Lophiotremataceae ott5345476, Lophiotrema ott700491, Occultibambusaceae ott7513523, Brunneofusispora ott7513524, Acericola ott7513569, mrcaott8428ott33433, mrcaott404584ott3720804, Stenella (genus in Nucletmycea) ott4075947, Dothideales ott406068, Dothioraceae ott1025079, Aureobasidium ott252445, mrcaott1741ott44265, mrcaott1741ott9815, mrcaott1741ott262725, mrcaott1741ott155480, mrcaott1741ott3675, mrcaott4995ott20900, Lasiosphaeriaceae ott234778, Lasiosphaeris ott254782, mrcaott95328ott144112, mrcaott95328ott652855, Chaetosphaeriales ott778723, Chaetosphaeriaceae ott737940, Porosphaerella (genus in Nucletmycea) ott946160, mrcaott5957ott31687, mrcaott5957ott8405, Valsaceae ott103004, Cytospora ott4052334, Togniniaceae ott4053923, mrcaott15510ott324849, mrcaott15510ott324855, mrcaott5476ott89139, mrcaott5476ott113472, mrcaott5476ott219359, Daldinia ott183725, Kretzschmaria ott404550, Nemania ott1075136, Diatrypaceae ott183729, Diatrype ott974055, Leotiomycetes ott346134, Leotiomycetidae ott5344705, mrcaott3396ott5089, mrcaott3396ott15644, mrcaott15644ott132902, Mollisia ott30120, Everhartia ott981918, mrcaott3913ott34649, mrcaott34649ott279567, Lachnaceae ott5345451, Lachnum ott472177, mrcaott328846ott442023, Hyphodiscus ott770426, Chlorociboria ott839008, Strossmayeria ott841104, Xylogramma ott4043168, Cadophora ott34262, mrcaott3264ott77029, mrcaott3264ott41259, mrcaott3264ott8087, mrcaott3264ott9155, mrcaott9155ott9466, mrcaott359518ott633514, mrcaott9669ott129242, mrcaott9669ott15646, mrcaott15646ott688524, Niessliaceae ott59361, Monocillium ott744148, Stictidaceae ott968418, Cryptodiscus ott843698, Phlyctidaceae ott199502, Phlyctis (genus in Opisthokonta) ott969101, Pertusariaceae ott216191, Ochrolechiaceae ott966001, mrcaott4573ott12794, mrcaott4573ott6404, mrcaott4573ott652857, mrcaott4573ott89967, mrcaott4573ott129215, mrcaott4573ott861697, mrcaott4573ott36640, Phialemonium ott669443, Cryptendoxyla ott908974, Natantiella ott666441, mrcaott14840ott112654, Trichosphaeriales ott123358, Trichosphaeriaceae ott123364, Brachysporium ott3705328, Rhodoveronaea ott214865, Rhamphoria ott645069, mrcaott77437ott103601, Coniochaetales ott516157, Coniochaetaceae ott812336, mrcaott91069ott143827, Flavoparmelia ott733707, mrcaott158950ott351331, mrcaott158950ott922716, Punctelia ott794602, mrcaott311283ott3768257, mrcaott311283ott311547, mrcaott311547ott311549, mrcaott311547ott649527, mrcaott311547ott872116, mrcaott33894ott53417, mrcaott53417ott134178, Sphaeropsis (genus in Nucletmycea) ott45369, Ramalinaceae ott691002, Bacidia ott635255, Physciaceae ott216195, Hyperphyscia ott730650, Caliciaceae ott339438, Amandinea ott155861, Orbiliales ott972723, Hyalorbilia ott885176, Myriangiales ott1063917, Myriangiaceae ott295792, Myriangium ott295791, Candelariales ott418461, Candelaria (genus in Nucletmycea) ott1050470, Helminthosphaeriaceae ott183866, Spadicoides ott403131, Rhizocarpales ott5296864, Catillariaceae ott716652, Catillaria ott877127, Hysteriales ott109863, Hysteriaceae ott366416, Hysterobrevium ott67936, Clonostachys ott802728, Scoliciosporaceae ott5345128, Xylonomycetes ott644766, Symbiotaphrinales ott5673779, Symbiotaphrinaceae ott7520906, Symbiotaphrina ott436261, Gibellulopsis ott724689, Gondwanamycetaceae ott5345437, Custingophora ott300281, Torrentispora ott4054988, Ciliciopodium ott4057458, Pseudodiplococcium ott7519919, Saccharomycotina ott971714, Saccharomycetes (class in h2007-2) ott989999, Pichiaceae ott821913, Pichia ott858840, Diddensiella ott5287334, Saccharomycetaceae ott989994, Kuraishia ott140858, Saccharomycodaceae ott207065, Hanseniaspora ott309730, Scheffersomyces ott14464, mrcaott275595ott434852, mrcaott434852ott434857, mrcaott363173ott1070320, mrcaott363173ott882705, Schwanniomyces ott391109, Lipomycetaceae ott834206, Myxozyma ott638678, Taphrinomycetes ott921288, Taphrinomycetidae ott5292180, Taphrinales ott698728, Taphrinaceae ott698727, Taphrina ott958772

# }
```
