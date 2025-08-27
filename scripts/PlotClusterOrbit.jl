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
    "--name"
    help = "Name of the model. Default: eccentric"
    arg_type = String
    default = "eccentric"
    "--rmax"
    help = "Position bounds (in kpc)."
    arg_type = Float64
    default = 25.0
end

parsed_args = parse_args(tabargs)

const framepersec = parsed_args["framerate"]
const GALACTIC_FRAME = parsed_args["galactic_frame"]
const name_model = parsed_args["name"]
const rmax = parsed_args["rmax"]

const path_to_data = "path/to/data/" # Location of the folder containing the data
const path_to_plot = "path/to/plot/" # Location of the folder containing the plots

################################################################################################################
# Read the run's information
# Switch back to galactocentric frame if needed
################################################################################################################

function get_list_files()

    list_files = glob(path_to_data * name_model * "/output/out_*.txt")
    nbt = length(list_files)
    list_time = zeros(Int64, nbt)

    for i=1:nbt 
        time = parse(Int64, split(split(list_files[i], "_")[end],".")[1])
        list_time[i] = time
    end
    p = sortperm(list_time)

    return list_files[p], nbt
end


function plot_data()

    list_files, nbt = get_list_files()
    tab_Rc_Vc = zeros(Float64, nbt, 6)

    Threads.@threads for i=1:nbt
        namefile = list_files[i]
        data = readdlm(namefile)
        data_info = Float64.(data[1,:])

        pc_per_HU = data_info[2]
        kpc_per_HU = pc_per_HU/1000.0
        kms_per_HU = data_info[9]

        # Cluster's position (in HU)
        xcl = data_info[13]
        ycl = data_info[14]
        zcl = data_info[15]
        vxcl = data_info[16]
        vycl = data_info[17]
        vzcl = data_info[18] 

        tab_Rc_Vc[i, 1] = xcl * kpc_per_HU
        tab_Rc_Vc[i, 2] = ycl * kpc_per_HU
        tab_Rc_Vc[i, 3] = zcl * kpc_per_HU

        tab_Rc_Vc[i, 4] = vxcl * kpc_per_HU
        tab_Rc_Vc[i, 5] = vycl * kpc_per_HU
        tab_Rc_Vc[i, 6] = vzcl * kpc_per_HU

    end


    # Plot the orbit of cluster in the (x, y) place
    p = plot(tab_Rc_Vc[:,1] , tab_Rc_Vc[:,2],
            xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]",
            xlims=(-rmax, rmax), ylims=(-rmax, rmax), 
            title="Cluster's center", 
            framestyle=:box, label=false,
            aspect_ratio=1, size=(600,600), 
            marker=true,
            markerstrokewidth=0,
            markersize=2)

    display(p)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/ClusterOrbit_xy.pdf"
    savefig(p, namefile_pdf)


    # Plot the orbit of cluster in the (x, z) place
    p = plot(tab_Rc_Vc[:,1] , tab_Rc_Vc[:,3],
            xlabel=L"x"*" [kpc]", ylabel=L"z"*" [kpc]",
            xlims=(-rmax, rmax), ylims=(-rmax, rmax), 
            title="Cluster's center", 
            framestyle=:box, label=false,
            aspect_ratio=1, size=(600,600), 
            marker=true,
            markerstrokewidth=0,
            markersize=2)

    display(p)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/ClusterOrbit_xz.pdf"
    savefig(p, namefile_pdf)


    # Plot the orbit of cluster in the (y, z) place
    p = plot(tab_Rc_Vc[:,2] , tab_Rc_Vc[:,3],
            xlabel=L"y"*" [kpc]", ylabel=L"z"*" [kpc]",
            xlims=(-rmax, rmax), ylims=(-rmax, rmax), 
            title="Cluster's center", 
            framestyle=:box, label=false,
            aspect_ratio=1, size=(600,600), 
            marker=true,
            markerstrokewidth=0,
            markersize=2)

    display(p)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/ClusterOrbit_yz.pdf"
    savefig(p, namefile_pdf)


    # Plot the orbit of cluster in the (R, z) place
    p = plot(sqrt.(tab_Rc_Vc[:,1] .^2 .+ tab_Rc_Vc[:,2] .^2), tab_Rc_Vc[:,3],
            xlabel=L"R"*" [kpc]", ylabel=L"z"*" [kpc]",
            xlims=(0.0, rmax), ylims=(-rmax, rmax), 
            title="Cluster's center", 
            framestyle=:box, label=false,
            aspect_ratio=1, size=(300,600), 
            left_margin = [5mm 0mm], 
            marker=true,
            markerstrokewidth=0,
            markersize=2)

    display(p)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/ClusterOrbit_Rz.pdf"
    savefig(p, namefile_pdf)


    # Plot the orbit of cluster in 3D space
    p = scatter(tab_Rc_Vc[:,1], tab_Rc_Vc[:,2], tab_Rc_Vc[:,3],
            xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", zlabel=L"z"*" [kpc]",
            xlims=(-rmax, rmax), ylims=(-rmax, rmax),  zlims=(-rmax, rmax),
            title="Cluster's center", 
            framestyle=:box, label=false,
            aspect_ratio=1, size=(600,600), 
            camera=(35, 10),
            markerstrokewidth=0,
            markersize=2)

    scatter!(p, [0], [0], [0], label=false, color=:red)

    display(p)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/ClusterOrbit.pdf"
    savefig(p, namefile_pdf)

    
    # Gif of the orbit
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
                xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", zlabel=L"z"*" [kpc]",  
                framestyle=:box, labels=:false,
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                camera=(35, 10),
                title="t = "*string(time)*" Myr")

        scatter!([0], [0], [0], label=false, color=:red)


    end

    mkpath(path_to_plot * "gif_output/"*name_model*"/")
    namefile_gif = path_to_plot * "gif_output/"*name_model*"/StreamEvolution.gif"
    gif(anim, namefile_gif, fps = framepersec)

end

plot_data()