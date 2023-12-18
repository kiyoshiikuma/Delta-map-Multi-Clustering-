using Pkg

# include("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/Jump_install_file.jl")
# Jump_install()

function Jump_install()
    
    Pkg.add("JuMP")
    Pkg.add("GLPKMathProgInterface")
    Pkg.add("Clp")
    Pkg.add("Cbc")
    Pkg.add("Culp")
    Pkg.add("ECOS")
    Pkg.add("Ipopt")
    Pkg.add("NLopt")

    print("complete")
    
    return 0

end