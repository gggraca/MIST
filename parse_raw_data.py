import csv
import sys
csv.field_size_limit(sys.maxsize)

fname = "data.csv"
fdata = csv.DictReader(open(fname, "r"))
data = [x for x in fdata]

headers = ["m/z", "Intensity"]

for x in data:
    print("Writing results to", str(x["MassBankID"]) + ".csv")
    outname = "Experiment/" + str(x["MassBankID"]) + ".csv"
    output = csv.DictWriter(open(outname, "w", newline=""), fieldnames=headers)
    output.writeheader()
    spectrum = x['spectrum']
    spectrum_list = []
    for i in spectrum.split(" "):
        d = {}
        d["m/z"], d["Intensity"] = i.split(":")
        spectrum_list.append(d)
    output.writerows(spectrum_list)

sys.exit()
for x in data:
    print(x)
    break
    spectrum = x['spectrum']
    spectrum_list = []
    for i in spectrum.split(" "):
        d = {}
        d["m/z"], d["Intensity"] = i.split(":")
        spectrum_list.append(d)
    output.writerows(spectrum_list)

    break
