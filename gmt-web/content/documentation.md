# EXAMPLES

***

## Run pcr format SV calls:

<p class='terminal' markdown='1'>
tigra_sv washu_pilot2_trio_large_deletions.pcr pcr_filelist.txt
</p>

## Run breakdancer format SV calls:

<p class='terminal' markdown='1'>
tigra_sv -b breakdancer.sv /gscuser/tumor.bam /gscuser/normal.bam
</p>

# TIGRA_SV

***

## NAME
TIGRA_SV - a tool that conducts targeted local assembly of structural variants (SV).


## DESCRIPTION

TIGRA_SV is a program that conducts targeted local assembly of structural variants (SV) using
the iterative graph routing assembly (TIGRA) algorithm (L. Chen, unpublished). It takes as input
a list of putative SV calls and a set of bam files that contain reads mapped to a reference genome
such as NCBI build36. For each SV call, it assembles the set of reads that were mapped or partially
mapped to the region of interest (ROI) in the corresponding bam files. Instead of outputing a single
consensus sequence, TIGRA_SV outputs all the alternative alleles in the ROI as long as they received 
sufficient sequence coverage (usually >= 2x). It is shown that TIGRA_SV is quite effective at improving 
the SV prediction accuracy in short reads analysis and can produce accurate breakpoint sequences that 
are valuable to understand the origin, mechanism and pathology underlying the SVs.


## SYNOPSIS

tigra_sv <SV file> <a.bam> <b.bam> ...

tigra_sv <SV file> <bam_list_file>

## OPTIONS

<dl>
<dt>-A INT</dt>
<dd>Esimated maximal insert size [500]</dd>

<dt>-l INT</dt>
<dd>Flanking size [500]</dd>

<dt>-w INT</dt>
<dd>Pad local reference by additional [200] bp on both ends</dd>

<dt>-q INT</dt>
<dd>Only assemble reads with mapping quality > [1]</dd>

<dt>-N INT</dt>
<dd>Number of mismatches required to be tagged as poorly mapped [5]</dd>

<dt>-p INT</dt>
<dd>Ignore cases that have average read depth greater than [1000]</dd>

<dt>-Q INT </dt>
<dd>Minimal BreakDancer score required for analysis [0]</dd>

<dt>-I STRING</dt>
<dd>Read intermediate files from DIR instead of creating them</dd>

<dt>-L STRING</dt>
<dd>Ignore calls supported by libraries that contains (comma separated) STRING</dd>

<dt>-c STRING</dt>
<dd>Specify chromosome for position 2 to parallelzing job</dd>

<dt>-b</dt>
<dd>Check when the format is in breakdancer</dd>

<dt>-r</dt>
<dd>Whether write reference to a file with .ref.fa as the suffix</dd>

<dt>-d</dt>
<dd>Whether to write dumped reads to a file with .fa as the suffix</dd>

<dt>-R</dt>
<dd>Reference file location with the full path</dd>
</dl>

# AUTHOR

***

Xian Fan, Ken Chen, Lei Chen 


# DISCLAIMER

***

This software is provided "as is" without warranty of any kind.
