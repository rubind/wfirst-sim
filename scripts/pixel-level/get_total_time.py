from numpy import *
import commands
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from pixel_level_ETC2 import *
from copy import deepcopy


slew_time = 70.
dither_time = 25.425
generate_data = 0
#SN_counts = [69, 208, 403, 223, 327, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136]
SN_counts = [100] + [400]*4 + [100]*11


#SN_counts = [30, 70, 137, 224, 280, 300, 406, 87, 83, 90, 95, 85, 114, 112, 99]

# for pixel_scale, slice_in_pixels, PSF_source, read_noise_floor, white_noise, gal_flamb, dark_current, live_dithers, ref_dithers, effareafl, IFURfl, run_color in runs:

runs = [#(0.15, 1, "test", 4.0, 0),
    #(0.15, 1, "AiryObstruct", 4.0, 0),
    #(0.15, 1, "WebbPSF", 4.0, 0),
    #(0.075, 2, "WebbPSF", 4.0, 0),
    #(0.075, 2, "WebbPSF", 3, 0),
    #(0.15, 1, "WebbPSF", 4, 6e-19),
    #(0.15, 1, "WebbPSF", 4, 6e-19, 0.01, 4),
    #(0.1, 1.5, "WebbPSF", 4, 6e-19, 0.01, 4),
    #(0.075, 2, "WebbPSF", 4, 6e-19, 0.01, 4, "IFU_effective_area_150518_blue_lower.txt", "IFU_R_Content.txt"),
    #(0.075, 2, "WebbPSF", 4, 6e-19, 0.01, 4, "IFU_effective_area_160513.txt", "IFU_R_160121.txt"),
    #(0.075, 2, "WebbPSF", 4, 6e-19, 0.01, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.1, 1.5, "WebbPSF", 4, 6e-19, 0.01, 5, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 2, 15, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 4, 15, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 2, 12, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 4., 15, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 4., 15, 0, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_500.txt"),
    (0.05, 3, "WebbPSF", 4., 15, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    (0.2, 1, "Gauss_PSFs_1.2m", 4., 15, 6e-19, 0.003, 1, 4, "1.2m_area.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 4., 15, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 4., 15, 0, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 4., 15, 6e-19, 0.005, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 2., 15, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 3., 20, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.045, 3.333, "WebbPSF", 3., 20, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 4, 15, 6e-19, 0.01, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 2, 10, 6e-19, 0.002, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 2, 15, 6e-19, 0.01, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 2, 15, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 4, 15, 6e-19, 0.001, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    
    #(0.04, 3.75, "WebbPSF", 0, 30, 6e-19, 0.003, 1, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.04, 3.75, "WebbPSF", 0, 30, 6e-19, 0.003, 2, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    
    #(0.04, 3.75, "WebbPSF", 2, 30, 6e-19, 0.003, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    ##(0.04, 3.75, "WebbPSF", 2, 30, 6e-19, 0.003, 1, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.04, 3.75, "WebbPSF", 2, 30, 6e-19, 0.003, 2, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 0, 6e-19, 0.01, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 0, 6e-19, 0.01, 1, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.05, 3, "WebbPSF", 0, 6e-19, 0.01, 2, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 0, 15, 6e-19, 0.01, 1, 4, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 0, 15, 6e-19, 0.01, 1, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 0, 15, 6e-19, 0.01, 2, 8, "IFU_effective_area_160720.txt", "IFU_R_160720.txt"),
    #(0.075, 2, "WebbPSF", 2, 6e-19, 0.01, 4),
    #(0.075, 2, "WebbPSF", 3, 6e-19, 0.01, 4),
    #(0.075, 2, "WebbPSF", 4, 6e-19, 0.01, 4),
    #(0.075, 2, "WebbPSF", 5, 6e-19, 0.01, 4),
    #(0.075, 2, "WebbPSF", 6, 6e-19, 0.01, 4),
    #(0.075, 2, "WebbPSF", 7, 6e-19, 0.01, 4),
    
    #(0.15, 1, "WebbPSF", 3.),
    #(0.075, 2, "WebbPSF", 3.)
    #(0.15, 1, "NIC"),
    #(0.075, 2, "NIC")]:
]

#for pix_scale in arange(0.05, 0.1501, 0.005):
#    runs.append((pix_scale, 2, "WebbPSF", 4, 6e-19, 0.01, 4*(pix_scale/0.075)**2., "IFU_effective_area_160720.txt", "IFU_R_160720.txt"))

print runs

offset_par = 4
offset_perp = 0.

read_noise = None

run_colors = ['r', 'g', 'b', 'k', 'orange', 'cyan', 'm', 'y']
for i in range(len(runs)):
    try:
        runs[i] = runs[i] + (run_colors[i],)
    except:
        runs[i] = runs[i] + ('k',)

plt.figure(figsize = (12, 10))

for pixel_scale, slice_in_pixels, PSF_source, read_noise_floor, white_noise, gal_flamb, dark_current, live_dithers, ref_dithers, effareafl, IFURfl, run_color in runs:

    PSFs = initialize_PSFs(pixel_scales = [int(round(pixel_scale/0.005))],
                           slice_scales = [int(round(pixel_scale*slice_in_pixels/0.005))], PSF_source = PSF_source)
    redshifts = arange(0.15, 1.66, 0.1)[:len(SN_counts)]
    args = {"gal_flamb": lambda x:gal_flamb, "pixel_scale": pixel_scale, "slice_scale": slice_in_pixels*pixel_scale, "dark_current": dark_current, "mdl": 'hsiao', "PSFs": PSFs, "TTel": 260.,
            "IFURfl": IFURfl, "effareafl": effareafl, "min_wave": 4200., "max_wave": 21000., "read_noise": read_noise, "white_noise": white_noise, "read_noise_floor": read_noise_floor, "offset_par": offset_par, "offset_perp": offset_perp}

    if 0:
        args["key1"] = "obs_frame"
        args["key2"] = (10200, 12850)
        StoNs = [3.5, 6.0, 10]
    elif 0:
        args["key1"] = "rest_frame_mean_S/N"
        args["key2"] = (5000, 6000)
        StoNs = [3.5, 6.0, 10]
    elif 1:
        args["key1"] = "rest_frame_band_S/N"
        args["key2"] = (5000, 6000)
        StoNs = [14, 24, 40, 46.6/sqrt(ref_dithers)]
    elif 0:
        args["key1"] = "rest_frame_band_S/N"
        args["key2"] = (6000, 7000)
        StoNs = [3.5*25/sqrt(136.*live_dithers), 25*6/sqrt(136.*live_dithers), 25*10/sqrt(136.*live_dithers), 25./sqrt(ref_dithers)]
    elif 0:
        args["key1"] = "rest_frame_band_S/N"
        args["key2"] = (5000, 6000)
        StoNs = [3.5*45/sqrt(136.*live_dithers), 45*6/sqrt(136.*live_dithers), 45*10/sqrt(136.*live_dithers), 45./sqrt(ref_dithers)]
    else:
        args["key1"] = "rest_frame_band_S/N"
        args["key2"] = (5000, 6000)
        StoNs = [3.5*45/sqrt(136.*live_dithers), sqrt(68.)*sqrt(15.)/sqrt(live_dithers), sqrt(68.)*sqrt(15.)/sqrt(live_dithers), 45./sqrt(ref_dithers)]

    args["restframe_bins"] = args["key2"]

    total_time = 0.

    for SN_count, redshift in zip(SN_counts, redshifts):
        # eight at S/N 3.5 + three at S/N 6 + one at S/N 10 = 8 + 3*2 + 4 = 18
        exptime35 = solve_for_exptime(S_to_N = StoNs[0], redshift = redshift, **args)
        exptime6 = solve_for_exptime(S_to_N = StoNs[1], redshift = redshift, **args)
        exptime10 = solve_for_exptime(S_to_N = StoNs[2], redshift = redshift, **args)#StoNs[2], redshift = redshift, **args)
        exptimeref = solve_for_exptime(S_to_N = StoNs[3], redshift = redshift, **args)
        
        time_at_this_redshift_per_SN = (
            (exptime35*live_dithers + slew_time + dither_time*(live_dithers - 1))*0 +
            (exptime6*live_dithers + slew_time + dither_time*(live_dithers - 1))*1.25 +
            (exptime10*live_dithers + slew_time + dither_time*(live_dithers - 1))*1.25 +
            (exptimeref*ref_dithers + slew_time + dither_time*(ref_dithers - 1))*1)

        total_time += time_at_this_redshift_per_SN*SN_count


        toprint = [redshift, "%.1f: %.1f" % (StoNs[0], exptime35), "%.1f: %.1f" % (StoNs[1], exptime6), "%.1f: %.1f" % (StoNs[2], exptime10), "live=%0.1f" % live_dithers, "ref=%0.1f" % ref_dithers, SN_count, "%.1f" % time_at_this_redshift_per_SN,
                   "Days: %.1f" % (time_at_this_redshift_per_SN*SN_count/86400.), "%.1f" % (total_time/86400.), "offset %i %i" % (args["offset_par"], args["offset_perp"]),
                   "PSF: " + PSF_source, str(pixel_scale) + '"', str(slice_in_pixels) + " pixel/slice", "RN white/floor: " + str(white_noise) + "/" + str(read_noise_floor), "DC: " + str(dark_current), "gal_flamb: %.1g" % gal_flamb, args["key1"] + " " + str(args["key2"])]

        toprint = [str(item) for item in toprint]

        print '\t'.join(toprint) + " Total"*int(redshift == redshifts[-1])

        label = ["Resolution:" + IFURfl, "Eff Area: " + effareafl, "Ref Dithers: %.1f" % ref_dithers, "PSF: " + PSF_source, str(pixel_scale) + '" pixels', str(slice_in_pixels) + " pixel/slice", "Read Noise Floor: " + str(read_noise_floor), "Dark: " + str(dark_current)]
        label = [str(item) for item in label]

        plt.plot(redshift, time_at_this_redshift_per_SN*SN_count/86400., 'o', color = run_color, label = int(redshift == redshifts[0])*(', '.join(label)))

    
    print
    print

plt.xlabel("Redshift")
plt.ylabel("Days of Time at This Redshift")

plt.legend(loc = 'best', fontsize = 7)
plt.savefig("time_per_redshift.pdf", bbox_inches = 'tight')
plt.close()
