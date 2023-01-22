import requests
import pandas as pd
import io

condition_df = pd.read_csv("data/trait_df.csv")

# merge non-unique rsids
condition_df = condition_df.groupby(["rsid", "chromosome", "position"]).size().reset_index(name="count")
# sort by chromosome and position
condition_df = condition_df.sort_values(by=["chromosome", "position"])

k = 0

key = input("Enter your API key: ")

equilibrium_df = pd.DataFrame(columns=["rsid1", "rsid2", "r2"])

def get_r2(current_set):
    url = f"https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?snps={'%0A'.join(list(current_set))}&pop=ALL&r2_d=d&genome_build=grch38_high_coverage&token={key}"
    r = requests.get(url)
    if r.status_code != 200:
        return None
    # get text of response (not json)
    data = r.text
    print(data)
    # convert to dataframe
    data = pd.read_csv(io.StringIO(data), sep="\t")
    if "RS_number" not in data.columns:
        return None
    # convert from wide to long format
    data = pd.melt(data, id_vars=["RS_number"], var_name="r2", value_name="value")
    data.columns = ["rsid1", "rsid2", "r2"]
    data["r2"] = data["r2"].astype(float)
    # remove self comparisons
    data = data[data["rsid1"] != data["rsid2"]]
    # remove comparisons with r2 < 0.2
    data = data[data["r2"] >= 0.2]
    # remove symmetric comparisons
    data = data[~data[["rsid1", "rsid2"]].apply(frozenset, axis=1).duplicated()]
    
    return data

for chrom in list(pd.unique(condition_df["chromosome"]))[::-1]:
    current_set = set()
    chrom_df = condition_df[condition_df["chromosome"] == chrom]
    i = 0
    j = 0
    while j < len(chrom_df):
        print(i, j)
        if i == j:
            i += 1
            continue
        try:
            row1 = chrom_df.iloc[i]
            row2 = chrom_df.iloc[j]
        except IndexError:
            break
        if row1["position"] - row2["position"] > 200000:
            j += 1
            continue
        current_set.add(row1["rsid"])
        current_set.add(row2["rsid"])
        j += 1
        if len(current_set) > 20:
            data = get_r2(current_set)
            
            
            # add to equilibrium_df
            if data is not None:
                equilibrium_df = pd.concat([equilibrium_df, data])

            current_set.clear()
            print(chrom, j/len(chrom_df))
    if len(current_set) > 1:
        data = get_r2(current_set)
        
        # add to equilibrium_df
        equilibrium_df = pd.concat([equilibrium_df, data])

        current_set.clear()
        print(chrom, j/len(chrom_df))
equilibrium_df.to_csv("data/equilibrium_df.csv", index=False)