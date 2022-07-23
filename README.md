# The Modern Promethease

The Modern Promethease is a free and open source tool that replicates the functionality of the tool [Promethease](https://promethease.com/). This tool creates a sumamry of a user's data from the popular genotyping service 23andMe. Example output for the user [Lilly Mendel](https://www.snpedia.com/index.php/User:Lilly_Mendel) is provided in the `mendelgenome/` folder.

## Installation Instructions

```
cd <path where you want to install The Modern Promethease>
git clone https://github.com/ThomasKAtkins/TheModernPromethease.git
cd TheModernPromethease
```

## Running Instructions (using pre-built datasets)

Assumes Python3 and Pandas are installed. 

Download your raw data from 23andMe by following [these instructions](https://customercare.23andme.com/hc/en-us/articles/212196868-Accessing-Your-Raw-Genetic-Data). Then, edit the provided file by deleting every line that starts with `#` except for the line 

```
# rsid	chromosome	position	genotype
```

which should be changed to:

```
rsid	chromosome	position	genotype
```

Now, to generate a report, run the command

```
python3 generate_report.py data/snp_df.csv data/geno_df.csv <path to 23andMe file>
```

to generate two files, `report.html`, and `snpedia_data.csv`. Opening `report.html` will display a graphical report of the 23andMe data, while `snpedia_data.csv` contains the data in tabular format.

## Re-Generating the Datasets (advanced users)

One advantage of this tool over Promethease is the ability to update the datasets used to create the tool. We provide somewhat up-to-date versions of these datasets in the `data/` folder, but also provide the code to generate these files with the script `data/getSnpediaPages.R` (requires the `devtools`, `dplyr`, and `readr` libraries). The datasets can be generated using

```
cd data
rm snp_df.csv && rm geno_df.csv
Rscript getSnpediaPages.R <path to 23andMe file>
```

The 23andMe file is necessary so the R script knows which SNPedia pages to scrape. This will produce two new data files, `snp_df.csv` (contains data on the SNP pages) and `geno_df.csv` (contains data on the individual genotype pages).