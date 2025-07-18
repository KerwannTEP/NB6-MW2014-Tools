using DelimitedFiles
using ArgParse
using Glob
using Plots
using LaTeXStrings
using Plots.PlotMeasures
using HDF5
using SpecialFunctions
using Statistics

tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--isnap"
    help = "Snapshot index for IOM analysis"
    arg_type = Int64
    default = 1
end

parsed_args = parse_args(tabargs)

const isnap = parsed_args["isnap"]


const G_in_kpc_MSun_Myr = 4.49851e-12
const Mtot = 1.0 # HU

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

    tab_out_time = zeros(Float64, nbt)
    tab_Rc_Vc = zeros(Float64, nbt, 6)

    it = 1
    open(filename, "r") do io
        while !eof(io)
            line = readline(io)
            # Look for the header line
            if occursin("CLUSTER ORBIT", line)

                matches = collect(eachmatch(pattern, line))
                numbers = parse.(Float64, [m.match for m in matches])

                tab_out_time[it] = numbers[1]
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

    return tab_Rc_Vc, tab_out_time

end

################################################################################################################
# MWPotential2014 (Bovy 2015)
# https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract
################################################################################################################



# Parameters (in physical units)

const Mb = 4.501370 * 10^9 # Msun
const alpha = -1.8
const rc = 1.9 # kpc

const Md = 6.819390 * 10^10 # Msun
const ad = 3.0 # kpc
const bd = 280.0 # pc

const Mh = 8.092401 * 10^11 # Msun
const rs = 16.0 # kpc
const ch = 15.3

const gc = log(1+ch) - ch/(1.0+ch)
const sb = 0.5*(3.0-alpha)
const gamma_s = gamma(sb)


# Potentials (in HU)

const _G = 1.0 # HU

function gamma_low_sx(xr)

    return gamma_s - gamma(sb, xr)

end

function psi_bulge(r::Float64, rc_HU::Float64, Mb_HU::Float64)

    xr = (r/rc_HU)^2

    # Inner potential
    gamma_lower = gamma_s - gamma(sb, xr)
    psi_inner = -_G*Mb_HU/r * gamma_lower/gamma_s

    # Outer potential 
    gamma_upper = gamma(1.0-alpha/2.0, xr)
    psi_outer = -_G*Mb_HU/rc_HU * gamma_upper/gamma_s
     
    return psi_inner + psi_outer

end

function psi_disk(R::Float64, z::Float64, ad_HU::Float64, bd_HU::Float64, Md_HU::Float64)

    return -_G * Md_HU/sqrt(R^2 + (sqrt(z^2 + bd_HU^2) + ad_HU)^2)

end

function psi_halo(r::Float64, rs_HU::Float64, Mh_HU::Float64)

    x = r/rs_HU
    cst = _G*Mh_HU/rs_HU * 1.0/gc

    return -cst/x * log(1+x)

end

################################################################################################################
# Save important quantities and plots
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

function treat_data()

    ###############################################
    # Read data
    ###############################################

    data_global = readdlm("path/to/data/global.30")
    data_lag = readdlm("path/to/data/lagr.7")

    list_files, nbt = get_list_files()
    tab_Rc_Vc, tab_out_time_in_Myr = get_cluster_position() # IN KPC AND KM/S !!!!

    tab_dist_cluster_in_kpc = sqrt.(tab_Rc_Vc[:, 1] .^ 2 + tab_Rc_Vc[:, 2] .^ 2 + tab_Rc_Vc[:, 3] .^ 2)

    tab_time_in_HU = zeros(Float64, nbt)
    tab_E_in_HU = zeros(Float64, nbt)
    tab_L_in_HU = zeros(Float64, nbt, 3)

    tab_rcore_in_HU = zeros(Float64, nbt)
    tab_lag_in_HU = zeros(Float64, nbt, 5)
    tab_frac_unbound = zeros(Float64, nbt)

    # Get units
    data0 = readdlm(list_files[1])
    data0_info = Float64.(data0[1, :])
    data0_star = Float64.(data0[2:end, 1:8])

    pc_per_HU = data0_info[2]
    kpc_per_HU = pc_per_HU/1000.0
    Msun_per_HU = data0_info[3]
    Myr_per_HU = data0_info[8]
    kms_per_HU = data0_info[9]
    Npart = length(data0_star[:, 1])

    Nub = 0
    tab_E_wrt_cluster_in_HU = zeros(Float64, Npart)
    tab_Lz_wrt_cluster_in_HU = zeros(Float64, Npart)

    Threads.@threads for i=1:nbt
        namefile = list_files[i]
        data = readdlm(namefile)

        data_info = Float64.(data[1, :])
        data_star = Float64.(data[2:end, 1:8])

        # Cluster's position (convert to HU)
        xcl = tab_Rc_Vc[1] / kpc_per_HU
        ycl = tab_Rc_Vc[2] / kpc_per_HU
        zcl = tab_Rc_Vc[3] / kpc_per_HU
        vxcl = tab_Rc_Vc[4] / kms_per_HU
        vycl = tab_Rc_Vc[5] / kms_per_HU
        vzcl = tab_Rc_Vc[6] / kms_per_HU

        # Read stellar data
        tab_time_in_HU[i] = data_info[1]
        r_dens_in_HU = data_info[5:7]
        tab_rcore_in_HU[i] = data_info[10]

        # Read cluster information
        tab_m_in_HU = data_star[:, 1]
        tab_r_in_HU = data_star[:, 2:4]
        tab_v_in_HU = data_star[:, 5:7]
        tab_pot_in_HU = data_star[:, 8]

        Utot = 0.0
        Ktot = 0.0

        Lxtot = 0.0
        Lytot = 0.0
        Lztot = 0.0

        for k=1:Npart 
            xc, yc, zc = tab_r_in_HU[k, :]
            x = xc + xcl
            y = yc + ycl
            z = zc + zcl
            R = sqrt(x^2 + y^2)
            r = sqrt(R^2 + z^2)

            vxc, vyc, vzc = tab_v_in_HU[k, :]
            vx = vxc + vxcl
            vy = vyc + vycl
            vz = vzc + vzcl

            m = tab_m_in_HU[k]
            phin = tab_pot_in_HU[k]

            # Bulge
            rc_HU = rc / kpc_per_HU
            Mb_HU = Mb / Msun_per_HU
            psib =  psi_bulge(r, rc_HU, Mb_HU)

            # Disk
            ad_HU = ad / kpc_per_HU
            bd_HU = bd / pc_per_HU
            Md_HU = Md / Msun_per_HU
            psid = psi_disk(R, z, ad_HU, bd_HU, Md_HU)

            # Halo
            rs_HU = rs / kpc_per_HU
            Mh_HU = Mh / Msun_per_HU
            psih = psi_halo(r, rs_HU, Mh_HU)

            # Total galactic potential 
            psi_gal = psib + psid + psih

            # Potential energy of particle
            Ui = m * (phin + psi_gal)
            Utot += Ui 

            # Kinetic energy of particle
            Ki = 0.5 * m * (vx^2 + vy^2 + vz^2)
            Ktot += Ki 

            # Angular momentum
            Lxi = m*(y*vz-z*vy)
            Lyi = m*(z*vx-x*vz)
            Lzi = m*(x*vy-y*vx) 

            Lxtot += Lxi
            Lytot += Lyi
            Lztot += Lzi
            
            if (i == isnap)
                Ui_c = m * phin
                Ki_c = 0.5 * m * (vxc^2 + vyc^2 + vzc^2)
                Ei_c = Ki_c + Ui_c 
                if (Ei_c >= 0.0)
                    Nub += 1
                end
            end


        end

        Utot *= 0.5
        
        Etot = Ktot + Utot
        tab_E_in_HU[i] = Etot 

        tab_L_in_HU[i, 1] = Lxtot
        tab_L_in_HU[i, 2] = Lytot
        tab_L_in_HU[i, 3] = Lztot


        # Snapshot IOM
        if (i == isnap)
            index = 1

            for k=1:Npart 
                xc, yc, zc = tab_r_in_HU[k, :]
                x = xc + xcl
                y = yc + ycl
                z = zc + zcl
                R = sqrt(x^2 + y^2)
                r = sqrt(R^2 + z^2)

                vxc, vyc, vzc = tab_v_in_HU[k, :]
                vx = vxc + vxcl
                vy = vyc + vycl
                vz = vzc + vzcl

                m = tab_m_in_HU[k]
                phin = tab_pot_in_HU[k]

                # Bulge
                rc_HU = rc / kpc_per_HU
                Mb_HU = Mb / Msun_per_HU
                psib =  psi_bulge(r, rc_HU, Mb_HU)

                # Disk
                ad_HU = ad / kpc_per_HU
                bd_HU = bd / pc_per_HU
                Md_HU = Md / Msun_per_HU
                psid = psi_disk(R, z, ad_HU, bd_HU, Md_HU)

                # Halo
                rs_HU = rs / kpc_per_HU
                Mh_HU = Mh / Msun_per_HU
                psih = psi_halo(r, rs_HU, Mh_HU)

                # Total galactic potential 
                psi_gal = psib + psid + psih

                # Potential energy of particle
                Ui = m * (phin + psi_gal)
                Utot += Ui 

                # Kinetic energy of particle
                Ki = 0.5 * m * (vx^2 + vy^2 + vz^2)
                Ktot += Ki 

                # Angular momentum
                Lxi = m*(y*vz-z*vy)
                Lyi = m*(z*vx-x*vz)
                Lzi = m*(x*vy-y*vx) 
                
                
                Ui_c = m * phin
                Ki_c = 0.5 * m * (vxc^2 + vyc^2 + vzc^2)
                Ei_c = Ki_c + Ui_c 
                    
                if (Ei_c >= 0.0)
                    tab_E_wrt_cluster_in_HU[index] = Ki + Ui
                    tab_Lz_wrt_cluster_in_HU[index] = Lzi
                    index += 1
                end

            end

        end

    end

    ###############################################
    # Save in HDF5 file
    ###############################################

    mkpath("path/to/data/post_treatment")
    namefile_hf5 = "path/to/data/post_treatment/iom_cluster.hf5"
    file = h5open(namefile_hf5, "w")

    write(file, "data_Etot_wrt_host", tab_E_in_HU)
    write(file, "data_Rc", tab_Rc_Vc[:, 1:3])
    write(file, "data_Vc", tab_Rc_Vc[:, 4:6])

    write(file, "data_time_HU", tab_time_in_HU)
    write(file, "data_time_Myr", tab_time_in_HU .* Myr_per_HU)
    
    write(file, "data_Lx_wrt_host", tab_L_in_HU[:, 1])
    write(file, "data_Ly_wrt_host", tab_L_in_HU[:, 2])
    write(file, "data_Lz_wrt_host", tab_L_in_HU[:, 3])

    write(file, "data_unbound_frac", data_global[2:end,31] ./ Npart)

    write(file, "data_lagrange_rad_01_pc", data_lag[3:end, 5] .* pc_per_HU)
    write(file, "data_lagrange_rad_10_pc", data_lag[3:end, 8] .* pc_per_HU)
    write(file, "data_lagrange_rad_20_pc", data_lag[3:end, 9] .* pc_per_HU)
    write(file, "data_lagrange_rad_50_pc", data_lag[3:end, 12] .* pc_per_HU)
    write(file, "data_lagrange_rad_90_pc", data_lag[3:end, 16] .* pc_per_HU)
    write(file, "data_r_core_pc", data_lag[3:end, 20] .* pc_per_HU)

    write(file, "Npart", Npart)
    write(file, "pc_per_HU", pc_per_HU)
    write(file, "kpc_per_HU", kpc_per_HU)
    write(file, "G_in_kpc_MSun_Myr", G_in_kpc_MSun_Myr)
    write(file, "Msun_per_HU", Msun_per_HU)

    write(file, "Myr_per_HU", Myr_per_HU)
    write(file, "kms_per_HU", kms_per_HU)

    close(file)


    ###############################################
    # Plot data
    ###############################################

    # Total energy

    plt = plot(tab_time_in_HU .* Myr_per_HU, tab_E_in_HU, 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Energy [HU]", 
        xticks=0:200:5000,
        xminorticks=4,
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        frame=:box)

    display(plt)
    readline()

    # Fractional (total) energy
    dataFracE = abs.(1.0 .- tab_E_in_HU[2:end] ./ tab_E_in_HU[1])
    minE = 10.0^floor(Int64, max(log10(minimum(dataFracE)), -16))
    maxE = 10.0^(floor(Int64, log10(maximum(dataFracE)))+1)

    plt = plot(tab_time_in_HU[2:end] .* Myr_per_HU, dataFracE, 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Fractional energy", 
        yaxis=:log10,
        xticks=0:200:5000,
        xminorticks=4,
        yticks=10.0 .^ (-20:1:2),
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        ylims=(minE, maxE),
        frame=:box)

    display(plt)
    readline()

    # Total angular momentum Lz
    plt = plot(tab_time_in_HU .* Myr_per_HU, tab_L_in_HU[:, 3], 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel=L"L_z"*" [HU]", 
        xticks=0:200:5000,
        xminorticks=4,
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        frame=:box)

    display(plt)
    readline()

    # Distance to cluster
    dmax = maximum(tab_dist_cluster_in_kpc) * 1.1
    plt = plot(tab_out_time_in_Myr, tab_dist_cluster_in_kpc,
                labels=false,
                xlabel="Time [Myr]", ylabel="Distance to cluster [kpc]",
                xlims=(tab_out_time_in_Myr[1], tab_out_time_in_Myr[end]),
                xticks=0:200:5000,
                ylims=(0.0, dmax),
                xminorticks=4,
                marker=true, markersize=2,
                frame=:box)

    display(plt)
    readline()

    # Unbound particles
    plt = plot(data_global[2:end,1] .* Myr_per_HU, data_global[2:end,31] ./Npart .* 100.0,
                labels=false,
                xlabel="Time [Myr]", ylabel="Fraction of unbound stars [%]",
                xlims=(data_global[2,1], data_global[end,1]) .* Myr_per_HU,
                xticks=0:200:5000,
                xminorticks=4,
                yticks=0:10:100,
                yminorticks=5,
                frame=:box)

    display(plt)
    readline()

    # Lagrange radii and core radius

    plt = plot(data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,20] .* pc_per_HU, yaxis=:log10, label=L"R_{c}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,5] .* pc_per_HU, yaxis=:log10, label=L"r_{01}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,8] .* pc_per_HU, yaxis=:log10, label=L"r_{10}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,9] .* pc_per_HU, yaxis=:log10, label=L"r_{20}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,12] .* pc_per_HU, yaxis=:log10, label=L"r_{50}", legend=:bottomleft)
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,16] .* pc_per_HU, yaxis=:log10, label=L"r_{90}", legend=:topleft)
    plot!(plt, legend=:topleft, xlabel="Time [Myr]", ylabel="Radii [pc]", frame=:box)
    plot!(plt, yticks=10.0 .^ (-5:1:10), yminorticks=10)
    plot!(plt, xticks=0:200:5000, xminorticks=4)
    plot!(plt, xlims=(data_lag[3,1], data_lag[end,1]) .* Myr_per_HU)

    display(plt)
    readline()


    # Snapshot IOM
    meanE = mean(tab_E_wrt_cluster_in_HU[1:Nub])
    varE = var(tab_E_wrt_cluster_in_HU[1:Nub], corrected=true, mean=meanE)
    sigmaE = sqrt(varE)
    tabDeltaE =  (tab_E_wrt_cluster_in_HU[1:Nub] .- meanE) ./ sigmaE

    meanLz = mean(tab_Lz_wrt_cluster_in_HU[1:Nub])
    varLz = var(tab_Lz_wrt_cluster_in_HU[1:Nub], corrected=true, mean=meanLz)
    sigmaLz = sqrt(varLz)
    tabDeltaLz =  (tab_Lz_wrt_cluster_in_HU[1:Nub] .- meanLz) ./ sigmaLz

    println("<E>  [unbound, HU] = ", meanE)
    println("<Lz> [unbound, HU] = ", meanLz)

    time_snapshot_in_Myr = tab_time_in_HU[isnap] .* Myr_per_HU
    time_snapshot_in_Myr = round(time_snapshot_in_Myr, digits=1) # Cutoff digits

    s = 2.0

    plt = scatter(tabDeltaLz, tabDeltaE,
                xlabel=L"\Delta L_z", ylabel=L"\Delta E", 
                framestyle=:box, labels=:false, 
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                xticks=-3:0.5:3, yticks=-3:0.5:3,
                title="t = "*string(time_snapshot_in_Myr)*" Myr")

    display(plt)
    readline()

end

treat_data()