library(ectotemp)

#####################Thermoregulation effectiveness#####################

salsal <- read.delim("salsal.txt")
ichalp <- read.delim("ichalp.txt")
lishel <- read.delim("lishel.txt")
bufbuf <- read.delim("bufbuf.txt")
rantem <- read.delim("rantem.txt")

tespring <- na.omit(ichalp[,"te"])
teautumn <- na.omit(salsal[,"te"])
tbsal <- na.omit(salsal[,"tb"])
tbich <- na.omit(ichalp[,"tb"])
tblis <- na.omit(lishel[,"tb"])
tbbuf <- na.omit(bufbuf[,"tb"])
tbran <- na.omit(rantem[,"tb"])

## Bootstrapping to estimate confidence interval of E

bootsalsal <- bootstrap_E(teautumn, tbsal,
                          17.40, 20.98,
                          'blouin',
                          10000)

bootichalp <- bootstrap_E(tespring, tbich,
                          14.44, 18.33,
                          'blouin',
                          10000)

bootlishel <- bootstrap_E(tespring, tblis,
                          16.24, 21.26,
                          'blouin',
                          10000)

bootbufbuf <- bootstrap_E(tespring, tbbuf,
                          19.35, 26.44,
                          'blouin',
                          10000)

bootrantem <- bootstrap_E(tespring, tbran,
                          19.34, 26.70,
                          'blouin',
                          10000)

## E comparisons

salsalvsichalp <- compare_E(salsal, ichalp,
                            17.40, 20.98, 
                            14.44, 18.33,
                            'blouin',
                            10000)

salsalvslishel <- compare_E(salsal, lishel,
                            17.40, 20.98, 
                            16.24, 21.26,
                            'blouin',
                            10000)

salsalvsbufbuf <- compare_E(salsal, bufbuf,
                           17.40, 20.98, 
                           19.35, 26.44,
                           'blouin',
                           10000)

salsalvsrantem <- compare_E(salsal, rantem,
                            17.40, 20.98, 
                            19.34, 26.70,
                            'blouin',
                            10000)

ichalpvslishel <- compare_E(ichalp, lishel,
                            14.44, 18.33, 
                            16.24, 21.26,
                            'blouin',
                            10000)

ichalpvsbufbuf <- compare_E(ichalp, bufbuf,
                            14.44, 18.33, 
                            19.35, 26.44,
                            'blouin',
                            10000)

ichalpvsrantem <- compare_E(ichalp, rantem,
                            14.44, 18.33, 
                            19.34, 26.70,
                            'blouin',
                            10000)

lishelvsbufbuf <- compare_E(lishel, bufbuf,
                            16.24, 21.26, 
                            19.35, 26.44,
                            'blouin',
                            10000)

lishelvsrantem <- compare_E(lishel, rantem,
                            16.24, 21.26, 
                            19.34, 26.70,
                            'blouin',
                            10000)

bufbufvsrantem <- compare_E(bufbuf, rantem,
                            19.35, 26.44, 
                            19.34, 26.70,
                            'blouin',
                            10000)
