using CSV
using DataFrames
using DelimitedFiles

function convert_data()

    namefile = "path/to/data/cosmic_ic.csv"

    df = CSV.read(namefile, DataFrame, delim=',', header=false)

    tabm = df[:, 3]
    tabx = df[:, 6]
    taby = df[:, 7]
    tabz = df[:, 8]
    tabvx = df[:, 9]
    tabvy = df[:, 10]
    tabvz = df[:, 11]

    namefile_dat = "path/to/data/dat.10"
    writedlm(namefile_dat, [tabm tabx taby tabz tabvx tabvy tabz])

end

convert_data()