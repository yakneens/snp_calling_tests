import pybedtools
from pybedtools.scripts import venn_mpl
from pybedtools.contrib import venn_maker

rheos = pybedtools.BedTool("/Users/siakhnin/data/giab/RMNISTHS_30xdownsample.chr20.rheos_kde_dels.bed")
delly = pybedtools.BedTool("/Users/siakhnin/data/giab/RMNISTHS_30xdownsample.chr20.no_spikein.delly.pass_only.dels.bed")
giab = pybedtools.BedTool("/Users/siakhnin/data/giab/Personalis_1000_Genomes_deduplicated_deletions_chr20_gt500bp.bed")
venn_mpl.venn_mpl(rheos, delly, giab, labels=["rheos", "delly", "giab"])
venn_maker.venn_maker(["/Users/siakhnin/data/giab/RMNISTHS_30xdownsample.chr20.rheos_kde_dels.bed",
                       "/Users/siakhnin/data/giab/RMNISTHS_30xdownsample.chr20.no_spikein.delly.pass_only.dels.bed",
                       "/Users/siakhnin/data/giab/Personalis_1000_Genomes_deduplicated_deletions_chr20_gt500bp.bed"],
                      names=["rheos","delly","giab"],
                      run=True,
                      figure_filename='out.png',
                      additional_args=['euler.d=TRUE',
                                       'scaled=TRUE',
                                       'cat.col=c("red","blue","green")',
                                       'col=c("red","blue","green")',
                                       'imagetype="png"']
                      )
