from numpy import *
import numpy as np
import pickle as pickle
import sys
import subprocess
import multiprocessing as mp
from scipy.stats import scoreatpercentile

colors = {"R062": (1, 0.5, 1), "Z087": 'm', "Y106": 'b', "J129": 'g', "H158": 'orange', "F184": 'r', "K213": 'r', "Ground_g": 'm', "Ground_r": 'b', "Ground_i": 'c', "Ground_z": 'g', "Ground_Y": 'orange', "g": 'm', "r": 'b', "i": 'c', "z": 'g', "Y": 'orange', "W146": 'orange', "Euclid_Y": 'b', "Euclid_J": 'g', "Euclid_H": 'orange', "P100": 'k',
          "REX_1.03": "m", "REX_1.09": "b", "REX_1.16": "c", "REX_1.22": "g", "REX_1.28": "orange", "REX_1.34": "red", "REX_1.40": "brown", "REX_1.46": "k"
          }

tier_labels = {"DeepImaging": "Deep Imaging",
               "WideImaging": "Wide Imaging"}


def plot_a_SN(lc_data, daymax, plot_to_make, phase_not_date, redshift, plt, stacked_SNRs, SN_ind, flc = None):


    filts_used = unique(lc_data["filts"])
    filts_to_sort = [[item[1:], item] for item in filts_used]
    filts_to_sort.sort()
    filts_used = [item[1] for item in filts_to_sort]
    
    for filt in filts_used:
        inds = where(array(lc_data["filts"]) == filt)

        if phase_not_date:
            xvals = (array(lc_data["dates"])[inds] - daymax)/(1. + redshift)
        else:
            xvals = array(lc_data["dates"])[inds]
        if plot_to_make == "LC":
            plt.errorbar(xvals, array(lc_data["fluxes"])[inds] + (len(filt) == 1)*3, yerr = array(lc_data["dfluxes"])[inds], fmt = '.', capsize = 0, color = colors[filt], label = "$" + filt.replace("Ground_", "").replace("_", "") + "$: %.1f" % stacked_SNRs[filt[0]][SN_ind])
            if flc != None:
                for i in range(len(xvals)):
                    flc.write("%.3f  %.3f  %.3f  %s  %.3f\n" % (xvals[i], array(lc_data["fluxes"])[inds][i], array(lc_data["dfluxes"])[inds][i], filt, redshift))
        else:
            plt.plot(xvals, (2.5/log(10.))*array(lc_data["dfluxes"])[inds]/array(lc_data["fluxes"])[inds], '.', color = colors[filt], label = "$" + filt.replace("Ground_", "") + "$")
    if plot_to_make == "dmag":
        plt.ylim(0., 0.3)
    #ylim = plt.ylim()
    #plt.plot([daymax*(1 - phase_not_date)]*2, ylim, color = (0.8, 0.8, 0.8))
    #plt.ylim(ylim)
    plt.axvline(daymax*(1 - phase_not_date), color = (0.8, 0.8, 0.8))
    ylim = plt.ylim()
    

    for date in lc_data["IFS_dates"]:
        phase = (date - daymax)/(1. + redshift)
        if phase_not_date:
            xval = phase
        else:
            xval = date
        if phase < 75:
            print("xval", xval, ylim)
            plt.annotate('', xy=(xval, 0), xytext=(xval, ylim[1]*0.1),
                         arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8, headlength=10), zorder = -1)

            #plt.axvline(xval, color = 'k', linestyle = ':')

    if phase_not_date == 0:
        plt.xlim(lc_data["dates"].min() - 1, lc_data["dates"].max() + 1)
    plt.legend(loc = 'best', fontsize = 8)


def get_useful_redshifts(z_set, tier_set, redshifts, stacked_SNRs, survey_fields, suffix, working_dir, plt):
    useful_redshift_mask = {}

    for SNR_thresh, plot_color in zip((0., 20., 30., 40., 80., 120., 160), ['r', "orange", (0.5, 1, 0), 'g', 'c',  'b', 'm']):
        useful_redshifts = {}

        for tier in tier_set:
            for z in z_set:
                inds = where((survey_fields == tier)*(redshifts == z))
                SNRs = stacked_SNRs[inds]
                total_SNe = float(len(SNRs))
                if total_SNe == 0:
                    useful_redshifts[(tier, z)] = 0
                else:
                    inds = where((survey_fields == tier)*(redshifts == z)*(stacked_SNRs >= SNR_thresh))
                    thresh_SNe = float(len(stacked_SNRs[inds]))
                    plt.plot(z, thresh_SNe/total_SNe, 'o', label = str(SNR_thresh), color = plot_color)

                    useful_redshifts[(tier, z)] = (thresh_SNe/total_SNe > 0.7)

        print("useful_redshifts", useful_redshifts)
        print("survey_fields", survey_fields)
        print("HACK on malmquist cut!!!!!")
        useful_redshift_mask[SNR_thresh] = ones(len(redshifts)) #array([useful_redshifts[(survey_fields[i], redshifts[i])] for i in range(len(redshifts))])

    plt.savefig(working_dir + "/efficiency_by_redshift_" + suffix + ".pdf", bbox_inches = 'tight')
    plt.close()

    return useful_redshift_mask


def make_IFS_time_plot(SN_data, working_dir, nsne, outputname, plt):
    all_IFS_times = []
    IFS_times_by_SN = []
    color_by_SN = []
    for i in range(nsne):
        all_IFS_times.extend(SN_data["SN_observations"][i]["IFS_exptimes"])
        IFS_times_by_SN.append(sum(SN_data["SN_observations"][i]["IFS_exptimes"]))
        color_by_SN.append(SN_data["SN_observations"][i]["c"])

    color_by_SN = array(color_by_SN)
    IFS_times_by_SN = array(IFS_times_by_SN)
    cbins = linspace(min(color_by_SN) - 0.001, max(color_by_SN) + 0.001, 20)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    for i in range(len(cbins) - 1):
        inds = where((color_by_SN >= cbins[i])*(color_by_SN < cbins[i+1]))
        ax1.plot(mean(color_by_SN[inds]), sum(IFS_times_by_SN[inds]), 'o', color = 'b', label = "Time (s)"*(i==5))
        ax2.plot(mean(color_by_SN[inds]), len(IFS_times_by_SN[inds]), 'o', color = 'g', label = "Number"*(i==5))
    plt.legend(loc = 'best')
    plt.xlabel("Color")
    plt.savefig(working_dir + "/IFS_exptime_by_c_%s.pdf" % outputname, bbox_inches = 'tight')
    plt.close()

    print("Total IFS time", sum(all_IFS_times), "for N_spectra", len(all_IFS_times))

def make_IFS_date_plot(SN_data, working_dir, nsne, outputname, plt):
    all_IFS_dates = []
    for i in range(nsne):
        all_IFS_dates.extend(SN_data["SN_observations"][i]["IFS_dates"])

    if len(all_IFS_dates) > 0:
        plt.hist(all_IFS_dates, bins = arange(min(all_IFS_dates) - 0.5,
                                              max(all_IFS_dates) + SN_data["survey_parameters"]["tier_parameters"]["cadences"][0][0],
                                              SN_data["survey_parameters"]["tier_parameters"]["cadences"][0][0]), linewidth = 0)
    plt.savefig(working_dir + "/IFS_dates_%s.pdf" % outputname, bbox_inches = 'tight')
    plt.close()


def make_IFS_phase_plot(SN_data, working_dir, nsne, outputname, plt):
    plt.figure(1, figsize = (12, 200))
    plt.figure(2, figsize = (15, 10))

    max_IFC_redshift = 0
    for k in range(nsne):
        if SN_data["SN_observations"][k]["IFS_dates"] != []:
            max_IFC_redshift = max(SN_data["SN_table"]["redshifts"][k], max_IFC_redshift)
    

    redshift_set = sort(unique(SN_data["SN_table"]["redshifts"]))
    redshift_scales = (redshift_set - redshift_set.min())/(max_IFC_redshift - redshift_set.min())

    date_stacks = [[], [], []]
    phase_stacks = [[], [], []]

    for i in range(len(redshift_set)):
        for j in range(3):
            phases = []
            for k in range(nsne):
                if len(SN_data["SN_observations"][k]["IFS_dates"]) >=3  and SN_data["SN_table"]["redshifts"][k] == redshift_set[i]:
                    phases.append(
                        (SN_data["SN_observations"][k]["IFS_dates"][j] - SN_data["SN_table"]["daymaxes"][k])/(1 + SN_data["SN_table"]["redshifts"][k])
                        )
            
            plt.figure(1)
            plt.subplot(len(redshift_set), 3,1 + i*3 + j)
            plt.hist(phases, bins = arange(-20., 20., 1))
            plt.title("Phase of Spectrum %i Redshift %.2f" % (j+1, redshift_set[i]))

            if len(phases) > 0:
                plt.figure(2)
                plt.subplot(2, 3, j + 1)
                date_stacks[j].extend(array(phases)*(1. + redshift_set[i]))
                the_label = "Redshift %.1f" % redshift_set[i]
                the_label = the_label * int((i == 0) + (redshift_set[i] == max_IFC_redshift))

                plt.hist(date_stacks[j], color = (redshift_scales[i], 0, 1 - redshift_scales[i]), zorder = 3 - redshift_scales[i], linewidth = 0, bins = arange(-40., 16, 2.5), label = the_label)
                if redshift_set[i] == max_IFC_redshift:
                    plt.legend(loc = 'best')
                plt.title(["1st", "2nd", "3rd"][j] + " Spectrum")
                plt.xlabel("Observer Date - Maximum")

                plt.subplot(2, 3, j + 1 + 3)
                phase_stacks[j].extend(phases)
                plt.hist(phase_stacks[j], color = (redshift_scales[i], 0, 1 - redshift_scales[i]), zorder = 3 - redshift_scales[i], linewidth = 0, bins = arange(-15., 11, 1.5), label = the_label)
                if redshift_set[i] == max_IFC_redshift:
                    plt.legend(loc = 'best')
                plt.title(["1st", "2nd", "3rd"][j] + " Spectrum")
                plt.xlabel("Rest Date - Maximum (Phase)")


    plt.figure(1)
    plt.savefig(working_dir + "/IFS_phases_%s.pdf" % outputname, bbox_inches = 'tight')
    plt.close()

    plt.figure(2)
    plt.savefig(working_dir + "/IFS_date_stack_%s.pdf" % outputname, bbox_inches = 'tight')
    plt.close()


def plot_field(SN_data, working_dir, nsne, outputname, plt):
    plt.figure(figsize=(20,20))
    #survey_fields_set = list(unique(SN_data["SN_table"]["survey_fields"]))
    #print "Survey Fields", survey_fields_set

    print(SN_data["observation_table"])

    for i in range(nsne):
        if SN_data["SN_observations"][i]["found_date"] != None:
            plt.subplot(3,3, 1 + int(around(SN_data["SN_table"]["RAs"][i]/20.)))#survey_fields_set.index(SN_data["SN_table"]["survey_fields"][i]))

            plt.plot(SN_data["SN_table"]["RAs"][i], SN_data["SN_table"]["Decs"][i], '.', color = 'b', markersize = 1)

    
    for i in range(len(SN_data["observation_table"]["RA"])):
        plt.subplot(3,3, 1 + int(around(SN_data["observation_table"]["RA"][i]/20.))
                    )
        if SN_data["observation_table"]["SNind"][i] == -1:
            color = 'r'
        else:
            color = 'g'

        plt.plot(SN_data["observation_table"]["RA"][i], SN_data["observation_table"]["dec"][i], '.', color = color, markersize = 1)

    
    plt.savefig(working_dir + "/survey_pointings_%s.pdf" % outputname, bbox_inches = 'tight')

    plt.close()

def plot_time_used(SN_data, working_dir, outputname, plt):
    

    dates = sort(unique(SN_data["observation_table"]["date"]))
    filts = sort(unique(SN_data["observation_table"]["filt"]))[::-1]
    labeled = []

    fields = around(SN_data["observation_table"]["RA"]/20.)
    nfields = len(unique(fields))

    plt.figure(figsize = (8, 6*nfields))
    for field in range(nfields):
        plt.subplot(nfields, 1, field+1)
        for date in dates:
            tot_time = 0
            for filt in filts:
                inds = where((SN_data["observation_table"]["date"] == date)*(SN_data["observation_table"]["filt"] == filt)*(fields == field))
                
                this_time = sum(SN_data["observation_table"]["exptime"][inds])
                if this_time > 0:
                    tot_time += this_time
                    plt.plot(date, tot_time, '.', color = colors[filt], label = filt*int(labeled.count(filt) == 0))
                    labeled.append(filt)
            
                assert all(SN_data["observation_table"]["date"][inds] == date)
        plt.legend(loc = 'best', fontsize = 9)
    plt.savefig(working_dir + "/exposure_time_%s.pdf" % outputname, bbox_inches = 'tight')
    plt.close()

def phases_of_discovery(working_dir, SN_data, outputname, plt, phase_not_observer = 1):
    ntiers = len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])

    plt.figure(figsize=(6,3*ntiers))
    for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]):

        redshift_set = sort(unique(SN_data["SN_table"]["redshifts"]))
        
        for redshift in redshift_set:
            phases = []
            not_found = 0.

            inds = where((SN_data["SN_table"]["redshifts"] == redshift)*(SN_data["SN_table"]["survey_fields"] == tier_name))[0]
            
            for ind in inds:
                daymax = SN_data["SN_table"]["daymaxes"][ind]
                found_date = SN_data["SN_observations"][ind]["found_date"]
                
                if found_date == None:
                    not_found += 1
                else:
                    if phase_not_observer:
                        phases.append((found_date - daymax)/(1. + redshift))
                    else:
                        phases.append(found_date - daymax)

            plt.subplot(ntiers, 2, 2*i+1)
            plt.plot(redshift, len(phases)/float(len(phases) + not_found + 1e-20), 'o', color = 'b')
            plt.ylim(0,1)
            plt.title("Fraction Found")
            
            plt.subplot(ntiers, 2, 2*i+2)
            plt.plot(redshift, scoreatpercentile(phases, 10), 'o', color = 'b')
            plt.plot(redshift, scoreatpercentile(phases, 50), 'o', color = 'g')
            plt.plot(redshift, scoreatpercentile(phases, 90), 'o', color = 'r')
            plt.title("10/50/90")

    plt.savefig(working_dir + "/" + "phase"*phase_not_observer + "date"*(1 - phase_not_observer) + "_of_discovery_" + outputname + ".pdf", bbox_inches = 'tight')
    plt.close()

            

def plot_time_remaining(SN_data, working_dir, outputname, plt):
    print(SN_data["time_remaining_values"])

    """
    plt.figure()
    plt.plot([item[0] for item in SN_data["time_remaining_values"]],
             [item[1] for item in SN_data["time_remaining_values"]], '.', color = 'b')
    plt.savefig(working_dir + "/timeremaining_" + outputname + ".pdf", bbox_inches = 'tight')
    plt.close()
    """


def light_curve_cuts(SN_data, nsne, crit, include_z = True):
    this_useful_redshift_mask = []
    all_reasonable_SNe = []
    
    for ind in range(nsne):
        lc_data = SN_data["SN_observations"][ind]

        good_SN = len(lc_data["fluxes"]) > 5
        reasonable_SN = len(lc_data["fluxes"]) > 5

        good_SN *= any([len(item) == 4 for item in lc_data["filts"]]) # Any WFIRST data?
        reasonable_SN *= any([len(item) == 4 for item in lc_data["filts"]])  # Any WFIRST data?

        if crit == "any point in any filter with S/N > 8":
            SNRs_all = lc_data["fluxes"]/lc_data["dfluxes"]
            if any(SNRs_all > 8):
                pass
            else:
                good_SN = 0

        for filt in list(unique(lc_data["filts"])):
            if len(filt) == 4 and ((filt != "Z087") or include_z): # WFIRST Filter 

                filtinds = where(array(lc_data["filts"]) == filt)
                SNRs = lc_data["fluxes"][filtinds]/lc_data["dfluxes"][filtinds]
                SNRs_cut = SNRs[where(SNRs > 2.)]
                crit_found = 0

                if crit == "any point in any filter with S/N > 8":
                    crit_found = 1


                if crit == ">= 1 point per filter with S/N > 10": # All WFIRST filters must have S/N > 10
                    crit_found = 1
                    if SNRs.max() < 10:
                        good_SN = 0

                if crit == "stacked S/N per filter > 5": # All WFIRST filters must have stacked S/N > 15
                    crit_found = 1
                    if sum(SNRs_cut**2.) < 5**2.:
                        good_SN = 0

                if crit == "stacked S/N per filter > 10": # All WFIRST filters must have stacked S/N > 15
                    crit_found = 1
                    if sum(SNRs_cut**2.) < 10**2.:
                        good_SN = 0

                if crit == "stacked S/N per filter > 15": # All WFIRST filters must have stacked S/N > 15
                    crit_found = 1
                    if sum(SNRs_cut**2.) < 15**2.:
                        good_SN = 0

                if crit == "stacked S/N per filter > 20": # All WFIRST filters must have stacked S/N > 20
                    crit_found = 1
                    if sum(SNRs_cut**2.) < 20**2.:
                        good_SN = 0

                for filt_to_find in ["R062", "Z087", "Y106", "J129", "H158", "F184"]:
                    if crit == ("stacked %s S/N > 20" % filt_to_find):
                        crit_found = 1

                        if filt == filt_to_find:
                            if sum(SNRs_cut**2.) < 20**2.:
                                good_SN = 0
                        if 1 - any(array(lc_data["filts"]) == filt_to_find):
                            good_SN = 0
                        
                assert crit_found == 1, "Couldn't find " + crit

        this_useful_redshift_mask.append(good_SN)
        all_reasonable_SNe.append(reasonable_SN)

    return array(this_useful_redshift_mask), array(all_reasonable_SNe)
                

def make_selection_figure(SN_data, working_dir, nsne, plt):
    survey_fields = array(SN_data["SN_table"]["survey_fields"])
    n_tiers = len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])
    extra_tier = n_tiers > 1

    crit_list = [">= 1 point per filter with S/N > 10",
                 "any point in any filter with S/N > 8",
                 "stacked S/N per filter > 5",
                 "stacked S/N per filter > 10",
                 "stacked S/N per filter > 15",
                 "stacked S/N per filter > 20"]
    for filt_to_find in ["R062", "Z087", "Y106", "J129", "H158", "F184"]:
        crit_list.append("stacked %s S/N > 20" % filt_to_find)
        
    n_crit = len(crit_list)

    plt.figure(figsize = (1+5*n_tiers + 5*extra_tier, 3*n_crit))
    f_crit = open(working_dir + "/redshifts_selection_crit.txt", 'w')
    

    for j in range(n_crit):
        this_useful_redshift_mask, all_reasonable_SNe = light_curve_cuts(SN_data, nsne, crit_list[j], include_z = False)
        f_crit.write(crit_list[j] + '\n')

        for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"] + ["All"]*extra_tier):
            plt.subplot(n_crit, n_tiers+extra_tier, i+1 + (n_tiers + extra_tier)*j)
            if tier_name != "All":
                inds = where((survey_fields == tier_name)*this_useful_redshift_mask)
                all_inds = where(survey_fields == tier_name)
            else:
                inds = where(this_useful_redshift_mask)
                all_inds = where(survey_fields != "AAAAAAAA")
            f_crit.write(tier_name + '\n')


            if len(inds[0]) > 0:
                plt.hist(SN_data["SN_table"]["redshifts"][all_inds], bins = arange(0., 3.01, 0.1), color = 'r')
                n_list, bins_list, NA = plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 3.01, 0.1), color = 'b', label = crit_list[j])
                plt.legend(loc = 'best', fontsize = 8)

                if crit_list[j] == "stacked S/N per filter > 15" and tier_name == "All":
                    f = open(working_dir + "/StoN_gt_15.txt", 'w')
                    f.write(str(list(n_list)) + '\n')
                    f.write(str(list(bins_list)) + '\n')
                    f.close()
                f_crit.write(str(list(n_list)) + '\n')
                f_crit.write(str(list(bins_list)) + '\n')

            plt.axvline(0.8, color = 'gray')
            plt.axvline(1.7, color = 'gray')

            plt.title(working_dir + " " + tier_name + ", NSNe: " + str(len(inds[0])), fontsize = 8)

    plt.savefig(working_dir + "/redshifts_selection_crit.pdf", bbox_inches = 'tight')
    plt.close()
    f_crit.close()



def make_SNR_vs_z(SN_data, working_dir, nsne, plt):
    print("Making SNR vs z...")

    survey_fields = array(SN_data["SN_table"]["survey_fields"])
    n_tiers = len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])
    extra_tier = n_tiers > 1

    filts = ["R062", "Z087", "Y106", "J129", "H158", "F184"]

    plt.figure(figsize = (1+7*n_tiers + 7*extra_tier, 5*len(filts)))
    fSNR = open(working_dir + "/max_SNR_vs_redshift.txt", 'w')

    for j in range(len(filts)):
        for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"] + ["All"]*extra_tier):
            print(tier_name, filts[j])

            plt.subplot(len(filts), n_tiers + extra_tier, i+1 + (n_tiers + extra_tier)*j)
            if tier_name != "All":
                inds = where(survey_fields == tier_name)[0]
            else:
                inds = where(survey_fields != "AAAAAAAA")[0]

            print("inds", inds)

            plot_zs = []
            plot_SNRs = []

            for ind in inds:
                lc_data = SN_data["SN_observations"][ind]

                if any(array(lc_data["filts"]) == filts[j]):
                    filtinds = where(array(lc_data["filts"]) == filts[j])
                    SNRs = lc_data["fluxes"][filtinds]/lc_data["dfluxes"][filtinds]
                    
                    
                    plot_zs.append(SN_data["SN_table"]["redshifts"][ind])
                    plot_SNRs.append(SNRs.max())
                    

            plt.plot(plot_zs, plot_SNRs, '.', color = 'b')
            for k in range(len(plot_zs)):
                fSNR.write("%s  %s  %f  %f\n" % (filts[j], tier_name, plot_zs[k], plot_SNRs[k]))

            plt.title(tier_name + ", " + filts[j], fontsize = 10)
            plt.yscale('log')
            plt.ylim(1, 100)

    plt.savefig(working_dir + "/max_SNR_vs_redshift.png", bbox_inches = 'tight')
    plt.close()
    fSNR.close()


def make_SNR_vs_z(SN_data, working_dir, nsne, plt):
    print("Making SNR vs host for F184...")

    survey_fields = array(SN_data["SN_table"]["survey_fields"])

    plt.figure(figsize = (8, 6))

    inds = where(survey_fields == "Deep")[0]

    print("inds", inds)

    plot_gals = []
    plot_SNRs = []

    for ind in inds:
        lc_data = SN_data["SN_observations"][ind]

        if any(array(lc_data["filts"]) == "J129") and abs(SN_data["SN_table"]["redshifts"][ind] - 1.3) < 0.1:
            filtinds = where(array(lc_data["filts"]) == "J129")
            SNRs = lc_data["true_fluxes"][filtinds]/lc_data["dfluxes"][filtinds]

            waveinds = where((SN_data["IFC_waves"] > 12500.)*(SN_data["IFC_waves"] > 13300.))
            flambmean = mean(SN_data["SN_observations"][ind]["gal_background"][waveinds])
            abref = 0.10884806248/12900**2.

            plot_gals.append(-2.5*log10(flambmean/abref))
            plot_SNRs.append(SNRs.max())


    plt.plot(plot_gals, plot_SNRs, '.', color = 'b')

    #plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1, 100)
    plt.xlabel("1.29 $\mu$m AB mag")

    plt.savefig(working_dir + "/max_SNR_vs_host.png", bbox_inches = 'tight')
    plt.close()

def get_cadence_stops(SN_data):
    print(SN_data["observation_table"])
    
    survey_fields = array(SN_data["observation_table"]["tier"])
    
    
    
    cadence_stops = {}

    for tier_name in SN_data["survey_parameters"]["tier_parameters"]["tier_name"]:
        last_date = 0.0

        #print(SN_data["observation_table"])
    
        
        inds = where((SN_data["observation_table"]["tier"] == tier_name))
        print("inds", inds)
        
        if len(inds[0]) > 0:
            dmax = SN_data["observation_table"]["date"][inds].max()
            print("tier ", tier_name, dmax)
            last_date = max(last_date, dmax)
        cadence_stops[tier_name] = last_date
    print("cadence_stops", cadence_stops)
    
    return cadence_stops


def make_lc_sampling(SN_data, working_dir, nsne, n_tiers, cadence_stops, plt):
    return 0

    plt.figure(figsize = (8, (n_tiers+1)*8))

    extra_tier = n_tiers > 1

    for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"] + ["All"]*extra_tier):
        plt.subplot(n_tiers+1, 1, i+1)
        this_useful_redshift_mask = (SN_data["SN_table"]["daymaxes"] < cadence_stops[tier_name] - 20*(1. + SN_data["SN_table"]["redshifts"]))
        inds = where((survey_fields == tier_name)*(stacked_SNRs[SNR_key] >= SNR_thresh)*this_useful_redshift_mask)
        
    


def write_pointings(SN_data, working_dir):
    try:
        SN_data["ordered_sol"]
    except:
        return 0

    f = open(working_dir + "/pointings.txt", 'w')
    f.write("date\tfilter\tRA1_from_cent\tdec1_from_cent\tRA2_from_cent\tdec2_from_cent\tdegrees\tseconds\n")

    for i in range(len(SN_data["ordered_sol"]["date"])):
        towrite = []
        for key in ["date", "filt", "RA1", "dec1", "RA2", "dec2", "deg", "time"]:
            towrite.append(str(SN_data["ordered_sol"][key][i]))
        f.write("\t".join(towrite) + '\n')
    f.close()


def make_SNRMAX_plots(SN_data, working_dir, plt, pickle_to_read):
    nsne = SN_data["nsne"]
    SNR_pass = []



    for i in range(nsne):
        SNRs = array(SN_data["SN_observations"][i]["fluxes"])/array(SN_data["SN_observations"][i]["dfluxes"])
        SNR_max_by_filt = []
        for filt in np.unique(SN_data["SN_observations"][i]["filts"]):
            inds = where(SN_data["SN_observations"][i]["filts"] == filt)            

            SNR_max_by_filt.append(np.max(SNRs[inds]))

        SNR_max_by_filt.sort()
        SNR_max_by_filt = SNR_max_by_filt[::-1]
        #print("SNR_max_by_filt", SNR_max_by_filt)

        if len(SNR_max_by_filt) > 2:
            if SNR_max_by_filt[0] >= 10 and SNR_max_by_filt[2] >= 5:
                SNR_pass.append(1)
            else:
                SNR_pass.append(0)
        else:
            SNR_pass.append(0)
            

    assert len(SN_data["SN_table"]["redshifts"]) == len(SNR_pass)

    SNR_pass = np.array(SNR_pass)
    
    plt.hist(SN_data["SN_table"]["redshifts"][np.where(SNR_pass == 1)], bins = np.arange(0, 3.01, 0.1))
    plt.title(pickle_to_read)
    plt.savefig(working_dir + "/redshift_SNRPEAK.pdf", bbox_inches = 'tight')


def plot_median_LCs(plt, SN_data, working_dir, stacked_SNRs, phase_not_date = 1):
    tiers_to_plot = SN_data["survey_parameters"]["tier_parameters"]["tier_name"]
    tiers_to_plot = [item for item in tiers_to_plot if item.count("Nearby") == 0]
    
    n_tiers = len(tiers_to_plot)
    survey_fields = array([str(item) for item in SN_data["SN_table"]["survey_fields"]])
    
    cadence_stops = get_cadence_stops(SN_data)

    print("cadence_stops", cadence_stops)
    print("survey_fields", survey_fields)

    z_to_plot = [0.225, 0.475, 0.725, 1.025, 1.475, 2.025]
    horizontal_layout = 0

    if horizontal_layout:
        plt.figure(figsize = (3.2*len(z_to_plot), 2.4*n_tiers))
    else:
        plt.figure(figsize = (3.2*n_tiers, 2.4*len(z_to_plot)))

    subplots_made = []
    for j, this_z in enumerate(z_to_plot):
        xlim = [1e7, -1e7]
        for i, tier_name in enumerate(tiers_to_plot):
            survey_mask = [item == tier_name for item in survey_fields]


            inds = where(isclose(SN_data["SN_table"]["redshifts"], this_z)*(survey_mask)*(SN_data["SN_table"]["daymaxes"] < cadence_stops[tier_name] - 20*(1. + this_z)))[0]
            print(tier_name, inds)

            S_to_N_inds = []
            for ind in inds:
                SNRs = SN_data["SN_observations"][ind]["true_fluxes"]/SN_data["SN_observations"][ind]["dfluxes"]
                SNRs = SNRs[where(SNRs > 1)]
                this_SNR = sqrt(sum(SNRs**2.))

                if this_SNR > 0:
                    S_to_N_inds.append((sqrt(sum(SNRs**2.)), ind))
            S_to_N_inds.sort()
            print("S_to_N_inds", S_to_N_inds)

            if len(S_to_N_inds) > 0:
                ind = S_to_N_inds[int(len(S_to_N_inds)/2.)][1]

                if horizontal_layout:
                    plt.subplot(n_tiers, len(z_to_plot), len(z_to_plot)*i+1 + j)
                    subplots_made.append(len(z_to_plot)*i+1 + j)
                else:
                    plt.subplot(len(z_to_plot), n_tiers, n_tiers*j+1 + i)
                    subplots_made.append(n_tiers*j+1 + i)
                    
                plot_a_SN(SN_data["SN_observations"][ind], SN_data["SN_table"]["daymaxes"][ind], plot_to_make = "LC",
                          phase_not_date = phase_not_date, redshift = this_z, plt = plt, stacked_SNRs = stacked_SNRs, SN_ind = ind)
                plt.xticks(fontsize = 8)
                plt.yticks(fontsize = 8)
                plt.ylim(0, plt.ylim()[1])
                
                xlim[0] = min(xlim[0], plt.xlim()[0])
                xlim[1] = max(xlim[1], plt.xlim()[1])
                #plt.text(tier_name + " z=%.2f" % this_z)

                if horizontal_layout:
                    if j == 0:
                        plt.ylabel("Flux, Zeropoint=25 AB")
                    if i == n_tiers - 1:
                        plt.xlabel("Relative Date " + "(Observer-Frame)"*(1 - phase_not_date) + "(Rest-Frame)"*phase_not_date)
                else:
                    if i == 0:
                        plt.ylabel("Flux, Zeropoint=25 AB")
                    if j == len(z_to_plot) - 1:
                        plt.xlabel("Relative Date " + "(Observer-Frame)"*(1 - phase_not_date) + "(Rest-Frame)"*phase_not_date)
            else:
                if horizontal_layout:
                    plt.subplot(n_tiers, len(z_to_plot), len(z_to_plot)*i+1 + j)
                else:
                    plt.subplot(len(z_to_plot), n_tiers, n_tiers*j+1 + i)
                plt.axis('off')



        for i in range(len(tiers_to_plot)):
            if horizontal_layout:
                plt.subplot(n_tiers, len(z_to_plot), len(z_to_plot)*i+1 + j)
                if subplots_made.count(len(z_to_plot)*i+1 + j):
                    plt.text(0.95*xlim[0] + 0.05*xlim[1], plt.ylim()[1]*0.95, tier_labels.get(tiers_to_plot[i], tiers_to_plot[i])
                             + " z=%.2f" % this_z, ha = 'center', va = 'top', rotation = 90)
                    plt.xlim(xlim)

            else:
                plt.subplot(len(z_to_plot), n_tiers, n_tiers*j+1 + i)
                if subplots_made.count(n_tiers*j+1 + i):
                    plt.text(0.95*xlim[0] + 0.05*xlim[1], plt.ylim()[1]*0.95, tier_labels.get(tiers_to_plot[i], tiers_to_plot[i])
                             + " z=%.2f" % this_z, ha = 'center', va = 'top', rotation = 90)
                    plt.xlim(xlim)
                
    plt.tight_layout()
    plt.savefig(working_dir + "/median_SN_LCs.pdf", bbox_inches = 'tight')
    plt.close()

    



def collection_of_plots(pickle_to_read):
    from matplotlib import use
    use("PDF")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages


    SN_data = pickle.load(open(pickle_to_read, 'rb'))
    
    if pickle_to_read.count("/"):
        working_dir = pickle_to_read[:pickle_to_read.rfind("/")] + "/output"
    else:
        working_dir = "output"

    outputname = pickle_to_read.split("/")[-1].replace("_pickle", "").replace("pickle", "").replace(".txt", "")
    print("outputname ", outputname)
    subprocess.getoutput("rm -fr " + working_dir)
    subprocess.getoutput("mkdir " + working_dir)

    for i in range(10):
        for j in range(len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])):
            SN_data["survey_parameters"]["tier_parameters"]["tier_name"][j] = SN_data["survey_parameters"]["tier_parameters"]["tier_name"][j].replace(str(i), "")
        for j in range(len(SN_data["SN_table"]["survey_fields"])):
            SN_data["SN_table"]["survey_fields"][j] = SN_data["SN_table"]["survey_fields"][j].replace(str(i), "")

    for key in SN_data:
        print("SN_data:%s" % key)
    print(SN_data["survey_parameters"])

    for key in SN_data:
        print(key)
    print() 
    for key in SN_data["SN_observations"][0]:
        print(key)

    #print SN_data["SN_observations"]
    print(SN_data["SN_table"])

    #write_pointings(SN_data, working_dir)


    plt.figure()
    plt.plot(SN_data["SN_table"]["redshifts"] + random.random(size = len(SN_data["SN_table"]["redshifts"]))*0.05 - 0.025, SN_data["SN_table"]["daymaxes"], '.')
    plt.xlabel("Redshift")
    plt.ylabel("Date of Maximum")
    plt.savefig(working_dir + "/" + outputname + "_" + "daymax_vs_redshift.pdf")
    plt.close()

    make_SNRMAX_plots(SN_data, working_dir, plt, pickle_to_read = pickle_to_read)

    z_set = list(set(list(SN_data["SN_table"]["redshifts"])))
    z_set.sort()


    nsne = SN_data["nsne"]
    assert len(SN_data["SN_table"]["redshifts"]) == nsne
    stacked_SNRs = {"All": [], "RZYJHFK": [], "9000": [], "10000": []}
    eff_lamb_dict = dict(Y106 = 10600, J129 = 12900, H158 = 15800, F184 = 18400, K213 = 21300) 

    for i in range(nsne):
        SNRs = array(SN_data["SN_observations"][i]["true_fluxes"])/array(SN_data["SN_observations"][i]["dfluxes"])
        inds = where(SNRs > 0)
        total_SNR = sqrt(dot(SNRs[inds], SNRs[inds]))
        stacked_SNRs["All"].append(total_SNR)

        for rest_lim in [9000, 10000]:
            total_SNR_squared = 0.
            for IR_filt in eff_lamb_dict:
                if eff_lamb_dict[IR_filt] / (1 + SN_data["SN_table"]["redshifts"][i]) >= rest_lim:
                    inds = np.where((SNRs > 0)*(SN_data["SN_observations"][i]["filts"] == IR_filt))
                    total_SNR_squared += np.dot(SNRs[inds], SNRs[inds])
            stacked_SNRs[str(rest_lim)].append(np.sqrt(total_SNR_squared))
            


        #IRinds = where((SN_data["SN_observations"][i]["filts"] == "K193") + (SN_data["SN_observations"][i]["filts"] == "F184") + (SN_data["SN_observations"][i]["filts"] == "H158"))
        #SNRs = array(SN_data["SN_observations"][i]["fluxes"][IRinds])/array(SN_data["SN_observations"][i]["dfluxes"][IRinds])
        #inds = where(SNRs > 0)
        #total_SNR = sqrt(dot(SNRs[inds], SNRs[inds]))
        #stacked_SNRs["HFK"].append(total_SNR)

        

        IRinds = where((SN_data["SN_observations"][i]["filts"] == "K213") + (SN_data["SN_observations"][i]["filts"] == "F184") + (SN_data["SN_observations"][i]["filts"] == "H158") + (SN_data["SN_observations"][i]["filts"] == "J129") + (SN_data["SN_observations"][i]["filts"] == "Y106") + (SN_data["SN_observations"][i]["filts"] == "Z087") + (SN_data["SN_observations"][i]["filts"] == "R062"))
        SNRs = array(SN_data["SN_observations"][i]["true_fluxes"][IRinds])/array(SN_data["SN_observations"][i]["dfluxes"][IRinds])
        inds = where(SNRs > 0)
        total_SNR = sqrt(dot(SNRs[inds], SNRs[inds]))
        stacked_SNRs["RZYJHFK"].append(total_SNR)
        
        for filt in set(SN_data["SN_observations"][i]["filts"]):
            if not filt[0] in stacked_SNRs:
                stacked_SNRs[filt[0]] = []

    for i in range(nsne):
        for filt in stacked_SNRs:
            if len(filt) == 1:
                IRinds = where(np.array([item[0] for item in SN_data["SN_observations"][i]["filts"]]) == filt)
                SNRs = np.array(SN_data["SN_observations"][i]["true_fluxes"][IRinds])/np.array(SN_data["SN_observations"][i]["dfluxes"][IRinds])
                inds = np.where(SNRs > 0)
                total_SNR = np.sqrt(np.dot(SNRs[inds], SNRs[inds]))

                stacked_SNRs[filt[0]].append(total_SNR)
            
        

    for key in stacked_SNRs:
        stacked_SNRs[key] = array(stacked_SNRs[key])

    print("stacked_SNRs", stacked_SNRs)
    plot_median_LCs(plt, SN_data, working_dir = working_dir, stacked_SNRs = stacked_SNRs)


    survey_fields = array(SN_data["SN_table"]["survey_fields"])
    print(survey_fields)


    n_tiers = len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])
    extra_tier = n_tiers > 1

    plot_time_used(SN_data, working_dir, outputname, plt)
    phases_of_discovery(working_dir, SN_data, outputname, plt, phase_not_observer = 1)
    phases_of_discovery(working_dir, SN_data, outputname, plt, phase_not_observer = 0)

    plot_time_remaining(SN_data, working_dir, outputname, plt)
    cadence_stops = get_cadence_stops(SN_data)

    useful_redshift_mask = {}
    for key in stacked_SNRs:
        useful_redshift_mask[key] = get_useful_redshifts(z_set, SN_data["survey_parameters"]["tier_parameters"]["tier_name"], SN_data["SN_table"]["redshifts"], stacked_SNRs[key], survey_fields, suffix = outputname, working_dir = working_dir, plt = plt)
    has_IFS_mask = array([len(SN_data["SN_observations"][i]["IFS_dates"]) > 0 for i in range(nsne)])

    #make_IFS_time_plot(SN_data, working_dir, nsne, outputname, plt = plt)
    #make_IFS_date_plot(SN_data, working_dir, nsne, outputname, plt = plt)
    #make_IFS_phase_plot(SN_data, working_dir, nsne, outputname, plt = plt)

    make_lc_sampling(SN_data, working_dir, nsne, n_tiers, cadence_stops = cadence_stops, plt = plt)
    #make_selection_figure(SN_data, working_dir, nsne, plt = plt)
    make_SNR_vs_z(SN_data, working_dir, nsne, plt)
    survey_fields = array([str(item).split("'")[1] for item in survey_fields])
    print("survey_fields", survey_fields)
    print("unique survey_fields", np.unique(survey_fields))
    

    for use_malm in [1]:
        for cumulative in [0]:
            for has_IFS in [0]:
                for SNR_key in stacked_SNRs:
                    plt.figure(figsize=(6,1+3.5*n_tiers + 3.5*extra_tier))

                    for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"] + ["All"]*extra_tier):
                        print("tier_name", tier_name)
                        plt.subplot(n_tiers+1, 1, i+1)
                        for SNR_thresh, SNR_color in zip((0., 10, 20., 30., 40., 80., 120., 160), ['k', 'r', "orange", (0.5, 1, 0), 'g', 'c',  'b', 'm']):

                            if use_malm:
                                this_useful_redshift_mask = 1
                            else:
                                this_useful_redshift_mask = useful_redshift_mask[SNR_key][SNR_thresh]


                            if has_IFS:
                                this_useful_redshift_mask *= has_IFS_mask

                            if tier_name != "All":
                                this_useful_redshift_mask *= (SN_data["SN_table"]["daymaxes"] < cadence_stops[tier_name] - 20*(1. + SN_data["SN_table"]["redshifts"]))
                                inds = where((survey_fields == tier_name)*(stacked_SNRs[SNR_key] > SNR_thresh)*this_useful_redshift_mask)
                                print(tier_name, "survey_fields == tier_name", sum(survey_fields == tier_name))
                                print(tier_name, "stacked_SNRs[SNR_key] > SNR_thresh", sum(stacked_SNRs[SNR_key] > SNR_thresh))
                                print(tier_name, "this_useful_redshift_mask", sum(this_useful_redshift_mask))
                                print(tier_name, "inds", inds)
                            else:
                                good_date = SN_data["SN_table"]["daymaxes"]*0
                                for j, tmp_tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]):
                                    good_date += (survey_fields == tmp_tier_name)*(SN_data["SN_table"]["daymaxes"] < cadence_stops[tmp_tier_name] - 20*(1. + SN_data["SN_table"]["redshifts"]))

                                inds = where((stacked_SNRs[SNR_key] > SNR_thresh)*good_date)

                            if len(inds[0]) > 0:
                                plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 3.01, 0.1), color = SNR_color, label = "SNR Sum>%.0f: %i"% (SNR_thresh, len(inds[0])), cumulative=cumulative, histtype="step"*(SNR_thresh == 0) + "stepfilled"*(SNR_thresh > 0))

                        plt.title("Tier: " + tier_labels.get(tier_name, tier_name))
                        plt.legend(loc = 'best', fontsize = 8)
                        plt.ylabel("Number of SNe Ia per 0.1")
                    plt.xlabel("Redshift")
                    plt.savefig(working_dir + "/redshifts_" + outputname + "_cumulative"*cumulative + "_nomalm"*(1 - use_malm) + "_hasIFS"*has_IFS + "_SNR_key=" + SNR_key + ".pdf", bbox_inches = 'tight')
                    plt.close()

    if SN_data["survey_parameters"]["tier_parameters"]["tier_name"].count("Wide"):
        plt.figure(figsize=(6,4))
        


        inds = where((stacked_SNRs["HFK"] >= 20.)*(survey_fields == "Wide"))#*(SN_data["SN_table"]["redshifts"] < 1.01))
        plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 2.01, 0.1), color = 'r', label = "Parallel-Observed, $F184$ S/N > 20, NSNe = %i" % len(SN_data["SN_table"]["redshifts"][inds]))

        inds = where((stacked_SNRs["RZYJHFK"] >= 40.)*(survey_fields == "Deep"))#*(SN_data["SN_table"]["redshifts"] < 1.01))
        plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 2.01, 0.1), color = 'g', label = "WFIRST Deep Tier, NSNe = %i" % len(SN_data["SN_table"]["redshifts"][inds]))
        
        inds = where(has_IFS_mask)#*(survey_fields == "Wide"))
        plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 2.01, 0.1), color = 'b', label = "IFC Observed, NSNe = %i" % len(SN_data["SN_table"]["redshifts"][inds]))
        
        plt.legend(loc = 'best')
        plt.xlabel("Redshift")
        plt.ylabel("Number of SNe")
        plt.xlim(0, 2)
        plt.savefig(working_dir + "/LSST.pdf", bbox_inches = 'tight')
        plt.close()


        

    plt.figure(figsize=(6,8))
    for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]):
        plt.subplot(len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]), 1, i+1)
        inds = where(survey_fields == tier_name)
        if len(inds[0]) > 0:
            plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 3.01, 0.1), color = 'r', label = str(len(inds[0])))

        mask = zeros(len(SN_data["SN_table"]["redshifts"]))
        for j in range(len(SN_data["SN_table"]["redshifts"])):
            if SN_data["SN_table"]["survey_fields"][j] == tier_name:
                daymax = SN_data["SN_table"]["daymaxes"][j]
                found_date = SN_data["SN_observations"][j]["found_date"]
                if found_date != None:
                    if (found_date - daymax) <= -10.:
                        mask[j] = 1


        inds = where(mask)
        plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 3.01, 0.1), color = 'b', label = str(len(inds[0])))
        plt.title(tier_name)
        plt.legend(loc = 'best', fontsize = 8)
    plt.savefig(working_dir + "/redshifts_for_trigger_" + outputname + ".pdf", bbox_inches = 'tight')
    plt.close()


    flc = open(working_dir + "/LC_summary.txt", 'w')

    for plot_to_make in ["LC", "dmag"]:
        print("Plotting LCs")
        pdf = PdfPages(working_dir + "/LC_samples_" + plot_to_make + "_" + outputname + ".pdf")



        for j in range(len(z_set)):
            plt.figure(figsize=(3*n_tiers, 9))
            for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]):
                #survey_mask = [item.decode('UTF-8') == tier_name for item in survey_fields]
                survey_mask = [item == tier_name for item in survey_fields]
                print(survey_fields)
                print(tier_name)
                #print("survey_mask", survey_mask)
                
                inds = where((SN_data["SN_table"]["redshifts"] == z_set[j])*(survey_mask)*(stacked_SNRs["All"] > 1)*(SN_data["SN_table"]["daymaxes"] < cadence_stops[tier_name] - 20*(1. + z_set[j])))[0]
                #print((SN_data["SN_table"]["redshifts"] == z_set[j]))
                #print((SN_data["SN_table"]["daymaxes"] < cadence_stops[i] - 20*(1. + z_set[j])))
                #print((SN_data["SN_table"]["redshifts"] == z_set[j])*(survey_mask)*(SN_data["SN_table"]["daymaxes"] < cadence_stops[i] - 20*(1. + z_set[j])))
                #print(sum((SN_data["SN_table"]["redshifts"] == z_set[j])*(survey_mask)*(SN_data["SN_table"]["daymaxes"] < cadence_stops[i] - 20*(1. + z_set[j]))))
                
                
                try:
                    sne_chosen = random.choice(inds, size = 3, replace = False)
                except:
                    sne_chosen = inds

                print("sne_chosen", sne_chosen, z_set[j])

                for k, ind in enumerate(sne_chosen):
                    plt.subplot(3, n_tiers, n_tiers*k+1 + i)

                    label_items = tier_name + " z=%.3f SNRs" %  z_set[j]
                    for key in ["All"]:#stacked_SNRs:
                        label_items += " " + key + "=%.1f" % stacked_SNRs[key][ind]

                    label_items += '\n'
                
                    label_items += "DayMax=%.1f RA=%.1f Dec=%.1f" % (SN_data["SN_table"]["daymaxes"][ind], SN_data["SN_table"]["RAs"][ind], SN_data["SN_table"]["Decs"][ind])

                    plt.title(label_items, size = 7)
                    if plot_to_make == "LC":
                        flc.write("_"*42 + '\n')
                    plot_a_SN(SN_data["SN_observations"][ind], SN_data["SN_table"]["daymaxes"][ind], plot_to_make = plot_to_make, phase_not_date = 1, redshift = z_set[j], plt = plt, flc = flc, stacked_SNRs = stacked_SNRs, SN_ind = ind)
                    plt.xticks(fontsize = 6)
                    plt.yticks(fontsize = 6)
                    plt.ylim(0, plt.ylim()[1])
            #plt.close()

            pdf.savefig(plt.gcf())
        pdf.close()
    #plot_field(SN_data, working_dir, nsne, outputname, plt = plt)
    flc.close()


if __name__ == "__main__":

    if len(sys.argv[1:]) > 1:
        processors_to_use = clip(len(sys.argv[1:]), 1, 32)
        
        pool = mp.Pool(processes = processors_to_use)
        pool.map(collection_of_plots, sys.argv[1:])
    else:
        collection_of_plots(sys.argv[1])

    print("Done!")
