# Presence of Methylation markers in Plasmodium vivax sal I
<br/>
We used two different basecalling softwares - Guppy and Megalodon to study the presence and distribution of epigenetic markers across the genome.

## Methylation distribution in exonic regions

<p float="left">
<img src='plots/meth dist. in exon - guppy.png' width='480px' height='400px' />
<img src='plots/meth dist. in exon - mega.png' width='480px' height='400px' />
</p>
<br/>

## Methylation distribution in intronic regions
<p float="left">
<img src='plots/meth dist. in intron - guppy.png' width='480px' height='400px' />
<img src='plots/meth dist. in intron - mega.png' width='480px' height='400px' />
</p>
<br/>


## Methylation distribution in intergenic regions
<p float="left">
<img src='plots/meth dist. in intergenic - guppy.png' width='480px' height='400px' />
<img src='plots/meth dist. in intergenic - mega.png' width='480px' height='400px' />
</p>
<br/>


## Comparison of Methylation distribution across all regions
<p float="left">
<img src='plots/meth freq in FY - guppy.png' width='480px' height='400px' />
<img src='plots/meth freq in FY - megalodon.png' width='480px' height='400px' />
</p>
<br/>

## Chromosome wise Methylation distribution in uncomplicated malaria
<p float="left">
<img src='plots/Mean meth in P.vivax (FY-guppy).png' width='480px' height='400px' />
<img src='plots/Mean meth in P.vivax (EL-guppy).png' width='480px' height='400px' />
</p>
<br/>

## Chromosome wise Methylation distribution in complicated malaria
<p float="left">
<img src='plots/Mean meth in P.vivax (EB-guppy).png' width='480px' height='400px' />
</p>
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

Thus, we have presented the variation and distribution of 5mC methylated CpG sites in intergenic, intronic and exonic regions in Plasmodium vivax sal I.