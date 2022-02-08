# mtbtyper
A tool for genotyping _Mycobacterium tuberculosis_ (Mtb) isolates from whole-genome sequencing (WGS) data

# requirements
python3 with numpy, pandas and [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/)

# input & output

Inputs: vcf files (`<sample>.vcf.gz`) in one directory, one file per sample.

Output: `lineage.csv` has the form

| sample_id | genotype | genotype_specific_snp | .. |
| --------- | -------- | --------------------- | -- |
| sample_1  | L2.2.M4.1 | L2.2 (24/24), L2.2.M4.1 (21/21), L2.2.Modern (6/6), L2 (6/6), L2.2.M4 (3/3) | .. |
| sample_2  | L2.2.M4 | L2.2 (24/24), L2.2.1 (10/10), L2.2.Modern (6/6), L2 (6/6), L2.2.M4 (3/3) | .. |
| sample_3  | L1.2.2.2 | L1.2 (56/56), L1.2.2 (53/53), L1.2.2.2 (19/19), L1 (15/15) | .. |
| ..  | .. | .. | .. |

The column `genotype_specific_snp` contains a list of counts of genotype-specific SNPs found in the sample according to the default scheme (see below). For instance, `L2.2 (24/24)` means the sample has 24 SNPs out of 24 L2.2-specific SNPs in this scheme.

The most specific sublineage level with a high proportion of genotype-specific SNPs present is returned as the predicted genotype in the column `genotype`.

Other columns are specific schemes from published studies (optional).


# usage

```
mtbtyper.py vcf_dir [options]
```

| option            | description |
| ----------------- | ----------- |
| `-o`, `--out`     | output directory (default: current working directory) |
| `-f`, `--fout`    | output file name (default: lineage.csv) |
| `-e`, `--vcf_end` | ending pattern of vcf file (default: vcf.gz) |
| `--all_schemes`   | add prediction from all available SNP schemes (default: false) |
| `--snpdb`         | path to genotyping SNP schemes (default: snpdb) |
| `--quiet`         | suppress screen output (default: false) |


# example

```
mtbtyper.py vcf -o lineage --all_schemes
```

In this example, the input vcf files are placed in a directory called `vcf`. The program outputs to `lineage/lineage.csv`, which include all available schemes (see below).


# SNP-typing schemes

The default scheme is a combination of the best scheme for each group of Mtb. It contains over 130 genotypes at different levels of classification hierarchy; see e.g. [Coll et al. (2014)](https://doi.org/10.1038/ncomms5812).

| genotypes            | source |
| -------------------- | ------ |
| L1 sublineages       | [Netikul et al. (2022)](https://doi.org/10.1038/s41598-022-05524-0) |
| L2 sublineages       | [Thawornwattana et al. (2021)](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000697) |
| L4.5.1               | [Mokrousov et al. (2017)](https://doi.org/10.1016/j.ympev.2017.09.002) |
| L4.5.2, L4.5.3       | [Ajawatanawong et al. (2019)](https://doi.org/10.1038/s41598-019-50078-3) |
| Animal-adapted lineages | [Lipworth et al. (2019)](https://wwwnc.cdc.gov/eid/article/25/3/18-0894_article) |
| Animal-adapted lineages + L6 | Unpublished |
| L8                   | [Napier et al. (2020)](https://doi.org/10.1186/s13073-020-00817-3) |
| Other L1-L7 lineages | [Coll et al. (2014)](https://doi.org/10.1038/ncomms5812) (diagnostic SNPs) |


Other schemes are based on individual published studies. Use `--all_schemes` flag to also output SNP counts from these schemes.

| scheme              | description | reference |
| ------------------- | ----------- | --------- |
| `coll2014`          | L1-L7 from a global collection | [Coll et al. (2014)](https://doi.org/10.1038/ncomms5812), available [here](https://datacompass.lshtm.ac.uk/id/eprint/414/) |
| `coll2014_diag`     | Diagnostic subset of `coll2014` | [Coll et al. (2014)](https://doi.org/10.1038/ncomms5812) |
| `coll2014_barcode`  | Barcoding subset of `coll2014_diag` | [Coll et al. (2014)](https://doi.org/10.1038/ncomms5812), Table S3 |
| `cr1`               | L1-L4 from Chiang Rai, Thailand | [Ajawatanawong et al. (2019)](https://doi.org/10.1038/s41598-019-50078-3) |
| `freschi2021`       | L1-L4; implemented in [fast-lineage-caller](https://github.com/farhat-lab/fast-lineage-caller) | [Freschi et al. (2021)](https://doi.org/10.1038/s41467-021-26248-1) |
| `freschi2021_hierarchical` | Same as `freschi2021` but with different genotype names | [Freschi et al. (2021)](https://doi.org/10.1038/s41467-021-26248-1) |
| `l1`                | L1 from a globally representative collection | [Netikul et al. (2022)](https://doi.org/10.1038/s41598-022-05524-0), Table S6 |
| `l1_barcode`        | Barcoding subset of `l1` | [Netikul et al. (2022)](https://doi.org/10.1038/s41598-022-05524-0), Table S6 |
| `l2`                | L2 from a globally representative collection | [Thawornwattana et al. (2021)](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000697), Table S7 |
| `l2_barcode`        | Barcoding subset of `l2` | [Thawornwattana et al. (2021)](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000697), Table S8 |
| `lipworth2019`      | L1-L6 and animal-adapted strains from a UK collection, implemented in [snpit](https://github.com/philipwfowler/snpit) | [Lipworth et al. (2019)](https://wwwnc.cdc.gov/eid/article/25/3/18-0894_article) |
| `merker2015`        | L2 from a global collection | [Merker et al. (2015)](https://doi.org/10.1038/ng.3195), Table S8 |
| `napier2020`        | Revised scheme of `coll2014` | [Napier et al. (2020)](https://doi.org/10.1186/s13073-020-00817-3), Table S2 |
| `napier2020_barcode`| Barcoding subset of `napier2020` | [Napier et al. (2020)](https://doi.org/10.1186/s13073-020-00817-3), Table S3 |
| `shitikov2017`      | L2 from a global collection | [Shitikov et al. (2017)](https://doi.org/10.1038/s41598-017-10018-5), Table S6 |

