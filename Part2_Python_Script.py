import os, sys, csv
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

def get_bounds(ar, tol=None):
    if not tol:
        tol = 10
    new = []
    for row in ar:
        mz = row[0]
        mz_lb, mz_ub = mz - (mz * tol / 1000000), mz + (mz * tol / 1000000)
        nrow = [row[0], row[1], mz_lb, mz_ub]
        new.append(nrow)

    return new

def validity_check(ar1, ar2):
    new = []
    ar1 = get_bounds(ar1, tol=10)
    for r1 in ar1:
        for r2 in ar2:
            if r2[0] > r1[2] and r2[0] < r1[3]:
                new.append(r1[0:2])
                break
    return new

def run_comparison(ft, outfile, input, ft_count):

    directory = ''
    count = 0
    for filename in input:
        scores = {}
        full_filename = os.path.join(directory, filename)

        if filename.endswith(".csv"):

            df = pd.read_csv(full_filename)
            new_ar1 = []
            new_ar2 = []
            ar1 = ft.to_numpy()
            print(ar1, "<--before")
            ar2 = df.to_numpy()
            new_ar1 = validity_check(ar1, ar2)
            ar1 = np.array(new_ar1)
            print(ar1, "<--after")
            if len(ar1) != 0:
                score = cosine_similarity(ar1, ar2)
                avg_score = np.mean(score)
                scores["Feature"] = ft_count
                count+=1
                print(count, filename)
                scores["Avg Coef"] = avg_score
                outfile.writerow(scores)
            else:
                scores["Avg Coef"] = 0.0
                outfile.writerow(scores)
            continue
        else:
            continue
    return scores, outfile

### Input handling
args = sys.argv
if len(args) > 1: #Add for loop for specification of multiple arguments
    input = [args[1]]
else:
    input = ["Experiment/MSJ00167.csv"]
    print("Warning: no input specified. Defaulting to", str(input[0]))
outname = str(input[0].split("/")[-1]).split(".")[0] + "_output.csv"
outfile = csv.DictWriter(open(outname, "w", newline=''), fieldnames=["Feature", "Avg Coef"])
outfile.writeheader()
print ("Results appending to:", outname)
results = {}
mastername = "features.csv"
masterfile = csv.DictReader(open(mastername, "r"))
master = [x for x in masterfile]
i = 0
masteroutname = mastername.split(".")[0] + "_results.csv"
headers = ["ID" ,"mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intb", "maxo", "sn", "sample", "coefscore"]
masterout = csv.DictWriter(open(masteroutname, "w"), fieldnames=headers, extrasaction='ignore')
masterout.writeheader()
for row in master:

    i += 1
    fname = "features" + str(i) + ".csv"
    outname = str(fname.split(".")[0]) + "_output.csv"

    df = pd.read_csv(open(fname))

    sub_df = df[['mz','into']]
    max = sub_df['into'].max()
    sub_df['Intensity'] = 100 * sub_df['into'] / max
    sub_df['m/z'] = sub_df['mz']
    sub_df = sub_df[['m/z', 'Intensity']]

    scores, outfile = run_comparison(sub_df, outfile, input, i)
    print(scores, "<------------scores")
    results[fname] = scores
    row["ID"] = row['']
    row["coefscore"] = scores['Avg Coef']
    masterout.writerow(row)
    print(fname)
print("Results appended to:", masteroutname, "\n\n")
