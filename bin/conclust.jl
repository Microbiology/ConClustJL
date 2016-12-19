#! /Applications/Julia-0.5.app/Contents/Resources/julia/bin/julia
# Geoffrey Hannigan
# University of Michigan

println("PROGRESS: Starting contig clustering.")

# Set dependencies
using Bio.Seq
using Clustering
using Distances
using Gadfly
using GaussianMixtures
using Cairo
using Fontconfig
using ArgParse

# Setup ArgParse
argx = ArgParseSettings()
@add_arg_table argx begin
    "--fasta", "-f"
        help = "Input fasta file for identifying circular contigs."
    "--abundance", "-a"
        help = "Input abundance table."
    "--output", "-o"
        help = "Output plot file."
    "--kmer", "-k"
        help = "Kmer lengths."
        default = 4
end

parsed_args = parse_args(ARGS, argx)

# Get files
fastafilepath = parsed_args["fasta"]
abundancefilepath = parsed_args["abundance"]
fileout = parsed_args["output"]
outputfile = open(fileout, "w")

println("PROGRESS: Reading in abundance table.")

# Read the abundance file into a dictionary
abundancetab = Dict()
coveratetablefile = open(abundancefilepath)
samplelength = 0

for line in eachline(coveratetablefile)
    covsplit = split(line, "\t")
    covname = shift!(covsplit)
    # Skip the header
    if covname == "contig"
        continue
    end
    covsplit = map(x->parse(Float64,x), covsplit)
    abundancetab[covname] = covsplit
    samplelength = length(covsplit)
end

println("PROGRESS: Creating k-mer table.")

kwindow = parsed_args["kmer"]

outmatrix = Matrix{Int64}(0, 4^kwindow + samplelength)

for s in open(FASTAReader{ReferenceSequence}, fastafilepath)
    result = zeros(Int64, 4^kwindow)
    sequence = DNASequence(s.seq)
    revseq = reverse_complement(sequence)
    for (pos,kmer) in each(DNAKmer{kwindow}, sequence)
        @inbounds result[reinterpret(Int64, kmer)+1] += 1
    end
    # Add reverse complement because we don't know the directionality
    # of the contigs so this is important.
    for (pos,kmer) in each(DNAKmer{kwindow}, revseq)
        @inbounds result[reinterpret(Int64, kmer)+1] += 1
    end
    result = Array{Float64}(result)
    append!(result, abundancetab[s.name])
    outmatrix = [outmatrix; result']
end

println("PROGRESS: Calculating distance matrix (Euclidean).")

R = pairwise(Euclidean(), outmatrix')

println("PROGRESS: Clustering by distance.")

kmatrix = Matrix{Int64}(0, 2)
for i in 5:5:100
    Rk = kmeans(outmatrix', i; maxiter=200)
    println(Rk.totalcost)
    karray = [i,Rk.totalcost]
    println(karray)
    kmatrix = [kmatrix; karray']
end

println("PROGRESS: Drawing plot.")

draw(PNG(fileout, 6inch, 3inch), plot(x = kmatrix[:, 1], y = kmatrix[:, 2], Geom.point, Geom.line))

Rfinal = kmeans(outmatrix', 40; maxiter=200)

println(Rfinal.counts)

println("PROGRESS: Complete.")
