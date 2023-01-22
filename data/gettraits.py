# Scrapes GWAS data from EBI API and frequency data from 1000 Genomes Project
# Saves to the file conditions.csv, if this file already exists, no need to run this script
# Takes many hours to run

import pandas as pd
import requests
import time
# surpress DtypeWarning
import warnings
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

snps = []

for i in range(576):
    url = "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms?page=" + str(i) + "&size=500"
    response = requests.get(url)
    data = response.json()
    pages = data['_embedded']['singleNucleotidePolymorphisms']
    for page in pages:
        snps.append(page['rsId'] + '\n')


geno_df = pd.read_csv("../mendelgenome/mendelgenome.txt", sep="\t")

union = set(geno_df["rsid"]).intersection(set(snps))
union = sorted(list(union))

alleles = []
ORs = []
betas = []
betadirections = []
betaunits = []
ranges = []
pvals = []
traits = []
RAFs = []
genes = []

start = time.time()

for j, snp in enumerate(union):
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}/associations?projection=associationBySnp"
    response = requests.get(url)
    if response.status_code != 200:
        continue
    data = response.json()
    if data["_embedded"]["associations"] == []:
        continue
    for i in range(len(data["_embedded"]["associations"])):
        if len(data["_embedded"]["associations"][i]["loci"][0]["strongestRiskAlleles"]) > 1:
            continue
        alleles.append(data["_embedded"]["associations"][i]["loci"][0]["strongestRiskAlleles"][0]["riskAlleleName"])
        ORs.append(data["_embedded"]["associations"][i]["orPerCopyNum"])
        betas.append(data["_embedded"]["associations"][i]["betaNum"])
        betadirections.append(data["_embedded"]["associations"][i]["betaDirection"])
        betaunits.append(data["_embedded"]["associations"][i]["betaUnit"])
        ranges.append(data["_embedded"]["associations"][i]["range"])
        RAFs.append(data["_embedded"]["associations"][i]["loci"][0]["strongestRiskAlleles"][0]["riskFrequency"])
        try:
            traits.append(data["_embedded"]["associations"][i]["efoTraits"][0]["trait"])
        except:
            traits.append(None)
        pvals.append(data["_embedded"]["associations"][i]["pvalue"])
        try:
            genes.append(data["_embedded"]["associations"][i]["loci"][0]["authorReportedGenes"][0]["geneName"])
        except:
            genes.append(None)
    # print progress and time remaining
    if j % 10 == 0:
        time_remaining = (time.time() - start) * (len(union) - j) / (j + 1) / 60
        print(f"{j/len(union):.3f}\t{time_remaining:.2f} minutes remaining")
    

rsids = [allele.split("-")[0] for allele in alleles]
# catch cases where there is no risk allele
risk_alleles = []
for allele in alleles:
    if len(allele.split("-")) == 1:
        risk_alleles.append(None)
    else:
        risk_alleles.append(allele.split("-")[1])

snp_df = pd.DataFrame({"rsid": rsids, "risk_allele": risk_alleles, "OR": ORs, "beta": betas, 
                     "beta_direction": betadirections, "beta_unit": betaunits, "range": ranges,
                     "RAF": RAFs, "trait": traits, "pval": pvals, "gene": genes})

snp_df = snp_df[pd.notna(snp_df["OR"])]
snp_df = snp_df[(snp_df.risk_allele == "A") |
                (snp_df.risk_allele == "C") |
                (snp_df.risk_allele == "G") |
                (snp_df.risk_allele == "T")]

start_time = time.time()

populations = ['SAMN10492705', 'SAMN10492695', 'SAMN10492703', 'SAMN10492696', 
               'SAMN10492698', 'SAMN10492704', 'SAMN10492697', 'SAMN10492701',
               'SAMN10492699', 'SAMN10492700', 'SAMN10492702', 'SAMN11605645']

rafs = {population: [] for population in populations}
i = 0
for rsid, risk_allele in zip(snp_df.rsid, snp_df.risk_allele):
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid[2:]}/frequency"
    response = requests.get(url)
    data = response.json()
    try:
        id = list(data["results"].keys())[0]
    except:
        for population in populations:
            rafs[population] += [None]
        continue
    for population in populations:
        try:
            total = sum(data["results"][id]["counts"]["PRJNA507278"]["allele_counts"][population].values())
            rafs[population] += [data["results"][id]["counts"]["PRJNA507278"]["allele_counts"][population][risk_allele] / total]
        except:
            rafs[population] += [None]
    i += 1
    if i % 10 == 0:
        time_remaining = (time.time() - start_time) / i * (len(snp_df) - i) / 60
        print(f"{i} / {len(snp_df)} ({time_remaining:.2f} minutes remaining)")

# add columns to snp_df
for population in populations:
    snp_df[population] = rafs[population]
snp_df.to_csv("conditions.csv", index=False)