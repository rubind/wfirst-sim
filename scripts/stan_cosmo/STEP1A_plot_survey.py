from numpy import *
import cPickle as pickle
import sys
import commands
import multiprocessing as mp
from scipy.stats import scoreatpercentile


def plot_a_SN(lc_data, daymax, plot_to_make, phase_not_date, redshift, plt):

    colors = {"Z087": 'm', "Y106": 'b', "J129": 'g', "H158": 'orange', "F184": 'r', "K193": 'r', "Ground_g": 'm', "Ground_r": 'b', "Ground_i": 'c', "Ground_z": 'g', "Ground_Y": 'orange', "g": 'm', "r": 'b', "i": 'c', "z": 'g', "Y": 'orange'}

    

    for filt in list(set(list(lc_data["filts"]))):
        inds = where(array(lc_data["filts"]) == filt)

        if phase_not_date:
            xvals = (array(lc_data["dates"])[inds] - daymax)/(1. + redshift)
        else:
            xvals = array(lc_data["dates"])[inds]
        if plot_to_make == "LC":
            plt.errorbar(xvals, array(lc_data["fluxes"])[inds] + (len(filt) == 1)*3, yerr = array(lc_data["dfluxes"])[inds], fmt = '.', capsize = 0, color = colors[filt], label = "$" + filt.replace("Ground_", "") + "$")
        else:
            plt.plot(xvals, (2.5/log(10.))*array(lc_data["dfluxes"])[inds]/array(lc_data["fluxes"])[inds], '.', color = colors[filt], label = "$" + filt.replace("Ground_", "") + "$")
    if plot_to_make == "dmag":
        plt.ylim(0., 0.3)
    #ylim = plt.ylim()
    #plt.plot([daymax*(1 - phase_not_date)]*2, ylim, color = (0.8, 0.8, 0.8))
    #plt.ylim(ylim)
    plt.axvline(daymax*(1 - phase_not_date), color = (0.8, 0.8, 0.8))

    for date in lc_data["IFS_dates"]:
        phase = (date - daymax)/(1. + redshift)
        if phase_not_date:
            xval = phase
        else:
            xval = date
        if phase < 75:
            plt.axvline(xval, color = 'k', linestyle = ':')
    plt.legend(loc = 'best', fontsize = 8)


def get_useful_redshifts(z_set, tier_set, redshifts, stacked_SNRs, survey_fields, suffix, working_dir, plt):

    useful_redshift_mask = {}

    for SNR_thresh, plot_color in zip((0., 20., 40., 80., 120.), ['r', "orange", 'g', 'c',  'b']):
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

        useful_redshift_mask[SNR_thresh] = array([useful_redshifts[(survey_fields[i], redshifts[i])] for i in range(len(redshifts))])

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

    print "Total IFS time", sum(all_IFS_times), "for N_spectra", len(all_IFS_times)

def make_IFS_date_plot(SN_data, working_dir, nsne, outputname, plt):
    all_IFS_dates = []
    for i in range(nsne):
        all_IFS_dates.extend(SN_data["SN_observations"][i]["IFS_dates"])

    if len(all_IFS_dates) > 0:
        plt.hist(all_IFS_dates, bins = arange(min(all_IFS_dates) - 0.5,
                                              max(all_IFS_dates) + SN_data["survey_parameters"]["tier_parameters"]["cadence"][0],
                                              SN_data["survey_parameters"]["tier_parameters"]["cadence"][0]), linewidth = 0)
    plt.savefig(working_dir + "/IFS_dates_%s.pdf" % outputname, bbox_inches = 'tight')
    plt.close()


def make_IFS_phase_plot(SN_data, working_dir, nsne, outputname, plt):
    plt.figure(figsize = (12, 200))
    redshift_set = sort(unique(SN_data["SN_table"]["redshifts"]))
    for i in range(len(redshift_set)):
        for j in range(3):
            phases = []
            for k in range(nsne):
                if SN_data["SN_observations"][k]["IFS_dates"] != [] and SN_data["SN_table"]["redshifts"][k] == redshift_set[i]:
                    phases.append(
                        (SN_data["SN_observations"][k]["IFS_dates"][j] - SN_data["SN_table"]["daymaxes"][k])/(1 + SN_data["SN_table"]["redshifts"][k])
                        )
            
            plt.subplot(len(redshift_set), 3,1 + i*3 + j)
            plt.hist(phases, bins = arange(-20., 20., 1))
            plt.title("Phase of Spectrum %i Redshift %.2f" % (j+1, redshift_set[i]))
                
    plt.savefig(working_dir + "/IFS_phases_%s.pdf" % outputname, bbox_inches = 'tight')
    plt.close()


def plot_field(SN_data, working_dir, nsne, outputname, plt):
    plt.figure(figsize=(20,20))
    #survey_fields_set = list(unique(SN_data["SN_table"]["survey_fields"]))
    #print "Survey Fields", survey_fields_set

    print SN_data["observation_table"]

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
    print SN_data["time_remaining_values"]

    """
    plt.figure()
    plt.plot([item[0] for item in SN_data["time_remaining_values"]],
             [item[1] for item in SN_data["time_remaining_values"]], '.', color = 'b')
    plt.savefig(working_dir + "/timeremaining_" + outputname + ".pdf", bbox_inches = 'tight')
    plt.close()
    """

def light_curve_cuts(SN_data, nsne, j):
    this_useful_redshift_mask = []
    all_reasonable_SNe = []
    
    for ind in range(nsne):
        lc_data = SN_data["SN_observations"][ind]

        good_SN = len(lc_data["filts"]) > 5
        reasonable_SN = len(lc_data["filts"]) > 5

        for filt in list(unique(lc_data["filts"])):
            filtinds = where(array(lc_data["filts"]) == filt)
            SNRs = lc_data["fluxes"][filtinds]/lc_data["dfluxes"][filtinds]
            SNRs_cut = SNRs[where(SNRs > 2.)]

            if j == 0: # All filters must have S/N > 10
                the_name = ">= 1 point per filter with S/N > 10"
                if SNRs.max() < 10:
                    good_SN = 0
            
            if j == 1: # All filters must have stacked S/N > 15
                the_name = "stacked S/N per filter > 15"
                if sum(SNRs_cut**2.) < 15**2.:
                    good_SN = 0

            if j == 2: # All filters must have stacked S/N > 20
                the_name = "stacked S/N per filter > 20"
                if sum(SNRs_cut**2.) < 20**2.:
                    good_SN = 0

            if j > 2 and j <= 7:
                filt_to_find = ["Z087", "Y106", "J129", "H158", "F184"][j-3]
                the_name = "stacked %s S/N > 20" % filt_to_find
 
                if filt == filt_to_find:
                    if sum(SNRs_cut**2.) < 20**2.:
                        good_SN = 0
                if 1 - any(array(lc_data["filts"]) == filt_to_find):
                    good_SN = 0
                        

        this_useful_redshift_mask.append(good_SN)
        all_reasonable_SNe.append(reasonable_SN)

    return array(this_useful_redshift_mask), array(all_reasonable_SNe), the_name
                

def make_selection_figure(SN_data, working_dir, nsne, plt):
    survey_fields = array(SN_data["SN_table"]["survey_fields"])
    n_tiers = len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])
    extra_tier = n_tiers > 1
    n_crit = 8

    plt.figure(figsize = (1+5*n_tiers + 5*extra_tier, 3*n_crit))

    for j in range(n_crit):
        this_useful_redshift_mask, all_reasonable_SNe, the_name = light_curve_cuts(SN_data, nsne, j)

        for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"] + ["All"]*extra_tier):
            plt.subplot(n_crit, n_tiers+extra_tier, i+1 + (n_tiers + extra_tier)*j)
            if tier_name != "All":
                inds = where((survey_fields == tier_name)*this_useful_redshift_mask)
                all_inds = where(survey_fields == tier_name)
            else:
                inds = where(this_useful_redshift_mask)
                all_inds = where(survey_fields != "AAAAAAAA")

            if len(inds[0]) > 0:
                plt.hist(SN_data["SN_table"]["redshifts"][all_inds], bins = arange(0., 2.6, 0.1), color = 'r')
                plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 2.6, 0.1), color = 'b', label = the_name)
                plt.legend(loc = 'best', fontsize = 8)

            plt.axvline(0.8, color = 'gray')
            plt.axvline(1.7, color = 'gray')

            plt.title(tier_name + ", NSNe: " + str(len(inds[0])))

    plt.savefig(working_dir + "/redshifts_selection_crit.pdf", bbox_inches = 'tight')
    plt.close()



def collection_of_plots(pickle_to_read):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages


    SN_data = pickle.load(open(pickle_to_read, 'rb'))
    if pickle_to_read.count("/"):
        working_dir = pickle_to_read[:pickle_to_read.rfind("/")] + "/output"
    else:
        working_dir = "output"

    outputname = pickle_to_read.split("/")[-1].replace("_pickle", "").replace("pickle", "").replace(".txt", "")
    print "outputname ", outputname
    commands.getoutput("mkdir " + working_dir)

    for i in range(10):
        for j in range(len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])):
            SN_data["survey_parameters"]["tier_parameters"]["tier_name"][j] = SN_data["survey_parameters"]["tier_parameters"]["tier_name"][j].replace(str(i), "")
        for j in range(len(SN_data["SN_table"]["survey_fields"])):
            SN_data["SN_table"]["survey_fields"][j] = SN_data["SN_table"]["survey_fields"][j].replace(str(i), "")

    for key in SN_data:
        print "SN_data:%s" % key
    print SN_data["survey_parameters"]

    for key in SN_data:
        print key
    print 
    for key in SN_data["SN_observations"][0]:
        print key

    #print SN_data["SN_observations"]
    print SN_data["SN_table"]

    z_set = list(set(list(SN_data["SN_table"]["redshifts"])))
    z_set.sort()


    nsne = SN_data["nsne"]
    stacked_SNRs = {"All": [], "ZYJHFK": [], "HFK": []}

    for i in range(nsne):
        SNRs = array(SN_data["SN_observations"][i]["fluxes"])/array(SN_data["SN_observations"][i]["dfluxes"])
        inds = where(SNRs > 0)
        total_SNR = sqrt(dot(SNRs[inds], SNRs[inds]))
        stacked_SNRs["All"].append(total_SNR)

        IRinds = where((SN_data["SN_observations"][i]["filts"] == "K193") + (SN_data["SN_observations"][i]["filts"] == "F184") + (SN_data["SN_observations"][i]["filts"] == "H158"))
        SNRs = array(SN_data["SN_observations"][i]["fluxes"][IRinds])/array(SN_data["SN_observations"][i]["dfluxes"][IRinds])
        inds = where(SNRs > 0)
        total_SNR = sqrt(dot(SNRs[inds], SNRs[inds]))
        stacked_SNRs["HFK"].append(total_SNR)

        IRinds = where((SN_data["SN_observations"][i]["filts"] == "K193") + (SN_data["SN_observations"][i]["filts"] == "F184") + (SN_data["SN_observations"][i]["filts"] == "H158") + (SN_data["SN_observations"][i]["filts"] == "J129") + (SN_data["SN_observations"][i]["filts"] == "Y106") + (SN_data["SN_observations"][i]["filts"] == "Z087"))
        SNRs = array(SN_data["SN_observations"][i]["fluxes"][IRinds])/array(SN_data["SN_observations"][i]["dfluxes"][IRinds])
        inds = where(SNRs > 0)
        total_SNR = sqrt(dot(SNRs[inds], SNRs[inds]))
        stacked_SNRs["ZYJHFK"].append(total_SNR)



    for key in stacked_SNRs:
        stacked_SNRs[key] = array(stacked_SNRs[key])


    survey_fields = array(SN_data["SN_table"]["survey_fields"])
    print survey_fields


    n_tiers = len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"])
    extra_tier = n_tiers > 1

    phases_of_discovery(working_dir, SN_data, outputname, plt, phase_not_observer = 1)
    phases_of_discovery(working_dir, SN_data, outputname, plt, phase_not_observer = 0)

    plot_time_remaining(SN_data, working_dir, outputname, plt)
    

    useful_redshift_mask = {}
    for key in stacked_SNRs:
        useful_redshift_mask[key] = get_useful_redshifts(z_set, SN_data["survey_parameters"]["tier_parameters"]["tier_name"], SN_data["SN_table"]["redshifts"], stacked_SNRs[key], survey_fields, suffix = outputname, working_dir = working_dir, plt = plt)
    has_IFS_mask = array([len(SN_data["SN_observations"][i]["IFS_dates"]) > 0 for i in range(nsne)])

    make_IFS_time_plot(SN_data, working_dir, nsne, outputname, plt = plt)
    make_IFS_date_plot(SN_data, working_dir, nsne, outputname, plt = plt)
    make_IFS_phase_plot(SN_data, working_dir, nsne, outputname, plt = plt)

    for use_malm in (0,1):
        for cumulative in (0,1):
            for has_IFS in (0,1):
                for SNR_key in stacked_SNRs:
                    plt.figure(figsize=(6,1+2.5*n_tiers + 2.5*extra_tier))

                    for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"] + ["All"]*extra_tier):
                        plt.subplot(n_tiers+1, 1, i+1)
                        for SNR_thresh, SNR_color in zip((0., 20., 40., 80., 120.), ['r', "orange", 'g', 'c',  'b']):

                            if use_malm:
                                this_useful_redshift_mask = 1
                            else:
                                this_useful_redshift_mask = useful_redshift_mask[SNR_key][SNR_thresh]

                            if has_IFS:
                                this_useful_redshift_mask *= has_IFS_mask

                            if tier_name != "All":
                                inds = where((survey_fields == tier_name)*(stacked_SNRs[SNR_key] >= SNR_thresh)*this_useful_redshift_mask)
                            else:
                                inds = where((stacked_SNRs[SNR_key] >= SNR_thresh)*this_useful_redshift_mask)

                            if len(inds[0]) > 0:
                                plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 2.6, 0.1), color = SNR_color, label = "SNR>%.0f: %i"% (SNR_thresh, len(inds[0])), cumulative=cumulative)

                        plt.title(tier_name)
                        plt.legend(loc = 'best', fontsize = 8)
                    plt.savefig(working_dir + "/redshifts_" + outputname + "_cumulative"*cumulative + "_nomalm"*(1 - use_malm) + "_hasIFS"*has_IFS + "_SNR_key=" + SNR_key + ".pdf", bbox_inches = 'tight')
                    plt.close()


    make_selection_figure(SN_data, working_dir, nsne, plt = plt)

    plt.figure(figsize=(6,8))
    for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]):
        plt.subplot(len(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]), 1, i+1)
        inds = where(survey_fields == tier_name)
        if len(inds[0]) > 0:
            plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 2.6, 0.1), color = 'r', label = str(len(inds[0])))

        mask = zeros(len(SN_data["SN_table"]["redshifts"]))
        for j in range(len(SN_data["SN_table"]["redshifts"])):
            if SN_data["SN_table"]["survey_fields"][j] == tier_name:
                daymax = SN_data["SN_table"]["daymaxes"][j]
                found_date = SN_data["SN_observations"][j]["found_date"]
                if found_date != None:
                    if (found_date - daymax) <= -10.:
                        mask[j] = 1


        inds = where(mask)
        plt.hist(SN_data["SN_table"]["redshifts"][inds], bins = arange(0., 2.6, 0.1), color = 'b', label = str(len(inds[0])))
        plt.title(tier_name)
        plt.legend(loc = 'best', fontsize = 8)
    plt.savefig(working_dir + "/redshifts_for_trigger_" + outputname + ".pdf", bbox_inches = 'tight')
    plt.close()



    for plot_to_make in ["LC", "dmag"]:
        print "Plotting LCs"
        pdf = PdfPages(working_dir + "/LC_samples_" + plot_to_make + "_" + outputname + ".pdf")



        for j in range(len(z_set)):
            plt.figure(figsize=(3*n_tiers, 9))
            for i, tier_name in enumerate(SN_data["survey_parameters"]["tier_parameters"]["tier_name"]):
                inds = where((SN_data["SN_table"]["redshifts"] == z_set[j])*(survey_fields == tier_name))[0]
                try:
                    sne_chosen = random.choice(inds, size = 3, replace = False)
                except:
                    sne_chosen = inds

                print "sne_chosen", sne_chosen, z_set[j]

                for k, ind in enumerate(sne_chosen):
                    plt.subplot(3, n_tiers, n_tiers*k+1 + i)

                    label_items = tier_name + " z=%.3f SNRs" %  z_set[j]
                    for key in stacked_SNRs:
                        label_items += " " + key + "=%.1f" % stacked_SNRs[key][ind]

                    label_items += '\n'
                
                    label_items += "DayMax=%.1f RA=%.1f Dec=%.1f" % (SN_data["SN_table"]["daymaxes"][ind], SN_data["SN_table"]["RAs"][ind], SN_data["SN_table"]["Decs"][ind])

                    plt.title(label_items, size = 7)
                    plot_a_SN(SN_data["SN_observations"][ind], SN_data["SN_table"]["daymaxes"][ind], plot_to_make = plot_to_make, phase_not_date = 1, redshift = z_set[j], plt = plt)
                    plt.xticks(fontsize = 6)
                    plt.yticks(fontsize = 6)
            #plt.close()

            pdf.savefig(plt.gcf())
        pdf.close()
    plot_field(SN_data, working_dir, nsne, outputname, plt = plt)

processors_to_use = 18


pool = mp.Pool(processes = processors_to_use)
pool.map(collection_of_plots, sys.argv[1:])


print "Done!"
