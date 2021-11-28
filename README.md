# Methylation distribution in NRC regions of _Maconellicoccus hirsutus_
<br/>

## Methylation distribution in over and under represented NRC regions

<img src='Plots/avg. meth in over rep NRC.png' width='700px' height='500px' />
<img src='Plots/Avg. meth in Under rep. NRC.png' width='700px' height='500px' />
<br/>

### Using ```Shapiro Wilk test``` to check for normality distribution
```> shapiro.test(subset(df, NRC_regions == "Over rep.")$Frequency)```
```
Shapiro-Wilk normality test

data:  subset(df, NRC_regions == "Over rep.")$Frequency
W = 0.86094, p-value < 2.2e-16
```
<br/>

```> shapiro.test(subset(df, NRC_regions == "Under rep.")$Frequency)```
```
Shapiro-Wilk normality test

data:  subset(df, NRC_regions == "Under rep.")$Frequency
W = 0.88901, p-value < 2.2e-16

```
From both histograms the methylation data does not appear to be normally distributed. This is further confirmed statistically by Shapiro-Wilk's test which shows p-value < 0.05.
Thus we reject the null hypothesis of normality for both distributions at the 5% significance level.
<br/>

## Methylation comparision between over and under represented regions of NRC

<img src='Plots/Over-Under freq. comparison.png' width='700px' height='500px' />

To statistically validate the differential methylation status between over and under represented NRC regions, we use the non-parametric ```Wilcoxon test```.

```wilcox.test(df$Frequency ~ df$NRC_regions)```

```
Wilcoxon rank sum test with continuity correction

data:  df$Frequency by df$NRC_regions
W = 1632757, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0
```
The p-value is less than 0.05, thus the null hypothesis is rejected and it is concluded that there is significant difference the methylation status of two NRC regions.
<br/>

To understand which genomic regions within NRC caused this differential expression, we analysed methylation in intergenic, exonic and intronic reigons.
## Distribution of methylation frequency in Intergenic regions of NRC

<img src='Plots/intergenic - over rep..png' width='700px' height='500px' />
<img src='Plots/Intergenic - under rep.png' width='700px' height='500px' />

### Using ```Shapiro Wilk test``` to check for normality distribution
```> shapiro.test(subset(df, NRC_regions == "Over rep.")$Frequency)```
```
Shapiro-Wilk normality test

data:  subset(df, NRC_regions == "Over rep.")$Frequency
W = 0.82099, p-value < 2.2e-16
```
<br/>

```> shapiro.test(subset(df, NRC_regions == "Under rep.")$Frequency)```
```
Shapiro-Wilk normality test

data:  subset(df, NRC_regions == "Under rep.")$Frequency
W = 0.87614, p-value < 2.2e-16
```
## Methylation comparision in intergenic regions of NRC

<img src='Plots/Intergenic distribution comparision-NRC.png' width='700px' height='500px' />

```
Wilcoxon rank sum test with continuity correction

data:  df$Frequency by df$NRC_regions
W = 801389, p-value = 2.46e-13
alternative hypothesis: true location shift is not equal to 0
```
The p-value is less than 0.05. Therefore, there is significant difference in the methylation status of two intergenic regions.
<br/>
