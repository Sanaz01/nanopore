#! /usr/bin/env python3

import sys
import csv
import argparse
import gzip

class SiteStats:
    def __init__(self, g_size, g_seq):
        self.num_reads = 0
        self.called_sites = 0
        self.called_sites_methylated = 0
        self.group_size = g_size
        self.sequence = g_seq

def update_call_stats(key, num_called_cpg_sites, is_methylated, sequence):
    if key not in sites:
        sites[key] = SiteStats(num_called_cpg_sites, sequence)

    sites[key].num_reads += 1
    sites[key].called_sites += num_called_cpg_sites
    if is_methylated > 0:
        sites[key].called_sites_methylated += num_called_cpg_sites

parser = argparse.ArgumentParser( description='Calculate methylation frequency at genomic CpG sites')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.0)
parser.add_argument('-s', '--split-groups', action='store_true')
args, input_files = parser.parse_known_args()
assert(args.call_threshold is not None)

sites = dict()
# iterate over input files and collect per-site stats
for f in input_files:
    if f[-3:] == ".gz":
        in_fh = gzip.open(f, 'rt')
    else:
        in_fh = open(f)
    csv_reader = csv.DictReader(in_fh, delimiter='\t')
    for record in csv_reader:

        num_sites = int(record['num_motifs'])
        llr = float(record['log_lik_ratio'])

        # Skip ambiguous call
        if abs(llr) < args.call_threshold * num_sites:
            continue
        sequence = record['sequence']

        is_methylated = llr > 0

        # if this is a multi-cpg group and split_groups is set, break up these sites
        if args.split_groups and num_sites > 1:
            c = str(record['chromosome'])
            s = int(record['start'])
            e = int(record['end'])

            # find the position of the first CG dinucleotide
            sequence = record['sequence']
            cg_pos = sequence.find("CG")
            first_cg_pos = cg_pos
            while cg_pos != -1:
                key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                update_call_stats(key, 1, is_methylated, "split-group")
                cg_pos = sequence.find("CG", cg_pos + 1)
        else:
            key = (str(record['chromosome']), int(record['start']), int(record['end']))
            update_call_stats(key, num_sites, is_methylated, sequence)

# header
print("\t".join(["chromosome", "start", "end", "num_motifs_in_group", "called_sites", "called_sites_methylated", "methylated_frequency", "group_sequence"]))

sorted_keys = sorted(list(sites.keys()), key = lambda x: x)

for key in sorted_keys:
    if sites[key].called_sites > 0:
        (c, s, e) = key
        f = float(sites[key].called_sites_methylated) / sites[key].called_sites
        print("%s\t%s\t%s\t%d\t%d\t%d\t%.3f\t%s" % (c, s, e, sites[key].group_size, sites[key].called_sites, sites[key].called_sites_methylated, f, sites[key].sequence))
        
        
        
l='&))(%%-+,#;,&-""'%+/',5(032+.*;:<,--#",331,95-8+"#%"&$#"%+,%%()2-&#:/'&>7&2/.1(*'.7?)J5,P314.+>?)A.4,.AI+-(*+/);'2'/76(1GO,8O5/,*7MA;'36&"""'*./)#J01:--;/)5<5I^]5,1L+*0.>G$;9,/G5)RM@22-+')73,5(4)/.+$&'*6>0;9,)%#.<B'JJQ-*..-$-"&$"$230A9(609K:(1?F,*$+77&%+)57*:,*9/,:>L.(<P30--"#"'#+)3:-N++.8%*6&83+/*?-=1&()00+&58)+8;+../-+11'5(,(5'.D,5+<N/5&:/%)'1+./;**)-/9,'CH,3(%3#('(*)'(6#+"-&$"##,#)-&,-#(#&&$-)2'5+4//++,-)*9,9'>($#/NF**&''%'5,3(%+.))+&2227*))(0562-/(,,,%*5.+/)+'23#%42/9(4*=$J66;2?)%,*?LN8-:*$%1/(..&/#R.2J<")*(3)4*6$+)1))-,)2BN><&C//+0*.?-#"$%/+207-A=)(.>;:'E101-)-3(0%'$&(*#C$04Q+/2(:-.:),00=GC2($(.(#8#"#'*(0@;+L0--CC2/.*2+35.-)/1(+-'-)-I-4LM1:-5:7578)4$#1(;29EGED.-9,0+'%'%#--0/-$8D-5.*%8,-4%)%#($$($-(6A>#//>C9:#,'+0-&2G((6%('""")099;-/4033D.#&7%%0220*&43"((%"15(7-6*8*2/0.5-H1+.LH=D/86,&??$-/0;57&3#0))0)*1+2);$&8WR80C/39)EF,C"/+1B*)"$20IL-1%6*G+NQI-%#,'()%%&749-7,&"%$'.*45++:;1.5:)1>*4+%&$$&.,8CO3+/$'<5)*',AD*1//&)&-+3JST1-.?+8N5.(*(&+'-,>*#$2+.3(#(&""#$%&#)93&1/=6*-3+*7)+?.(*%&+&,)2#"(*&(/7+$)4?(K3)'$,,&5,)85D5(%4';&)*0-(4(*?'/**$)#"#&$,+'0*)5%T+5E*M7*'.9.*0++&++BI.(9-(J3;-7.*JQRE#5+%*&-:4$6</1+""%/3*19319=)='/"("100+&%(-.:R,.(3.+(1+3S8)L*(7*30*6)(4*$$&)%)%(*-1%:'*//9(45*3,=LHF*%3,2%&:2+?(&-%+++D(+@;,&$(,7=+,/P#=-9)7G:E*2*+*6A00(LS47&7#"#*%&+.+-6->1Q5-&+)-C$=,(.84(**)$21.&7*$JN6*#&+')>-M'.>&:0($#'++*,)()"&-+*-43&"$'%8>5=,&%-%002(B0(,27%8?0+.92/##*46.'""  rl:i:0'
                                  
