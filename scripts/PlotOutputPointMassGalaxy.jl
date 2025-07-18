using DelimitedFiles
using ArgParse
using Glob
using Plots
using LaTeXStrings
using Plots.PlotMeasures
using LinearAlgebra

tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--framerate"
    help = "Number of frames per second."
    arg_type = Int64
    default = 60
    "--galactic_frame"
    help = "Convert output to galacto-centric frame."
    arg_type = Bool
    default = true
end

parsed_args = parse_args(tabargs)

const framepersec = parsed_args["framerate"]
const GALACTIC_FRAME = parsed_args["galactic_frame"]

################################################################################################################
# Read the cluster's position and velocity
# Hack to read the position in the OUT file (since RG, VG are not saved anywhere else a priori ?)
################################################################################################################

# In kpc
function get_cluster_position()

    filename = "path/to/data/OUT"

    # To extract numbers in scientific notation
    # https://docs.julialang.org/en/v1/base/strings/#Base.Regex
    # https://www.regular-expressions.info/floatingpoint.html
    # [+-]? → [...] means “any character in the set”, ? means “zero or one”.
    # \d → digit character.
    # + → one or more.
    # \. → escaped literal dot.
    # [Ee] → matches E or e.
    pattern = r"[+-]?\d+\.\d+[Ee][+-]?\d+"

    nbt = 0
    open(filename, "r") do io
        while !eof(io)
            line = readline(io)
            # Look for the header line
            if occursin("CLUSTER ORBIT", line)
                nbt += 1
            end
        end
    end


    tab_Rc_Vc = zeros(Float64, nbt, 6)

    it = 1
    open(filename, "r") do io
        while !eof(io)
            line = readline(io)
            # Look for the header line
            if occursin("CLUSTER ORBIT", line)

                matches = collect(eachmatch(pattern, line))
                numbers = parse.(Float64, [m.match for m in matches])

                tab_Rc_Vc[it, 1] = numbers[2]
                tab_Rc_Vc[it, 2] = numbers[3]
                tab_Rc_Vc[it, 3] = numbers[4]
                tab_Rc_Vc[it, 4] = numbers[5]
                tab_Rc_Vc[it, 5] = numbers[6]
                tab_Rc_Vc[it, 6] = numbers[7]

                it += 1
            end
        end
    end

    return tab_Rc_Vc

end

################################################################################################################
# Read the run's information
# Switch back to galactocentric frame if needed
################################################################################################################

function get_list_files()

    list_files = glob("path/to/data/output/out_*.txt")
    nbt = length(list_files)
    list_time = zeros(Int64, nbt)

    for i=1:nbt 
        time = parse(Int64, split(split(list_files[i], "_")[end],".")[1])
        list_time[i] = time
    end
    p = sortperm(list_time)

    return list_files[p], nbt
end

# temp 

function plot_data()

    list_files, nbt = get_list_files()
    tab_Rc_Vc = get_cluster_position()

    rmax = 25.0 # kpc

    # Plot the orbit of cluster
    p = scatter(tab_Rc_Vc[:,1], tab_Rc_Vc[:,2], tab_Rc_Vc[:,3],
            xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", zlabel=L"z"*" [kpc]",
            xlims=(-rmax, rmax), ylims=(-rmax, rmax),  zlims=(-rmax, rmax),
            title="Cluster's center", 
            framestyle=:box, label=false,
            aspect_ratio=1, size=(800,800), 
            camera=(35, 10),
            markerstrokewidth=0,
            markersize=3)

    scatter!(p, [0], [0], [0], label=false, color=:red)

    display(p)
    readline()

    mkpath("path/to/plot/fig/")
    namefile_pdf = "path/to/plot/fig/ClusterOrbit.pdf"
    savefig(p, namefile_pdf)

    anim = @animate for i=1:1:nbt
        namefile = list_files[i]
        data = readdlm(namefile)

        println("Progress : ", i, "/", nbt)

        data_info = Float64.(data[1,:])
        data_star = Float64.(data[2:end,1:8])

        # Read general information 
        time_in_HU = data_info[1]
        pc_per_HU = data_info[2]
        Msun_per_HU = data_info[3]
        r_tidal_in_HU = data_info[4]
        r_dens_in_HU = data_info[5:7] # 3D vector
        Myr_per_HU = data_info[8]
        kms_per_HU = data_info[9]
        r_core_in_HU = data_info[10]
        n_core_in_HU = data_info[11]
        v_core_in_HU = data_info[12]

        # Read cluster information
        tab_m_in_HU = data_star[:,1]
        tab_r_in_HU = data_star[:,2:4]
        tab_v_in_HU = data_star[:,5:7]
        tab_pot_in_HU = data_star[:,8]

        Npart = length(tab_m_in_HU)

        # Convert to kpc
        tab_r_in_kpc = tab_r_in_HU .* pc_per_HU ./ 1000.0

        if (GALACTIC_FRAME)

            tab_r_in_kpc[:, 1] = tab_r_in_kpc[:, 1] .+ tab_Rc_Vc[i, 1]
            tab_r_in_kpc[:, 2] = tab_r_in_kpc[:, 2] .+ tab_Rc_Vc[i, 2]
            tab_r_in_kpc[:, 3] = tab_r_in_kpc[:, 3] .+ tab_Rc_Vc[i, 3]

        end

        # Convert to Myr 
        time = time_in_HU * Myr_per_HU
        time = round(time, digits=1) # Cutoff digits

        # Markersize
        s = 1.0

        # Plot the cluster
        scatter(tab_r_in_kpc[:,1], tab_r_in_kpc[:,2], tab_r_in_kpc[:,3], 
                xlims=(-rmax, rmax), ylims=(-rmax, rmax), zlims=(-rmax, rmax),
                xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", 
                framestyle=:box, labels=:false,
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                camera=(35, 10),
                title="t = "*string(time)*" Myr")

        scatter!([0], [0], [0], label=false, color=:red)


    end

    mkpath("path/to/plot/gif_output/")
    namefile_gif = "path/to/plot/gif_output/StreamEvolution.gif"
    gif(anim, namefile_gif, fps = framepersec)

end

plot_data()