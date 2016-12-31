import glob
from numpy import *
import cPickle as pickle
import sys
import os
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/cosmo/")
from astro_functions import get_FoM
import pyfits
import multiprocessing as mp
from STEP1A_plot_survey import light_curve_cuts

def get_SNe_with_imaging_LCs(SN_data):
    nsne = SN_data["nsne"]

    SNe_at_15 = sum(light_curve_cuts(SN_data, nsne, crit = "stacked S/N per filter > 15")[0])
    SNe_at_20 = sum(light_curve_cuts(SN_data, nsne, crit = "stacked S/N per filter > 20")[0])

    return SNe_at_15, SNe_at_20

def get_survey_efficiency(SN_data):
    #for key in SN_data["SN_observations"]:
    #    print key


    observation_table = SN_data["observation_table"]

    inds = where((observation_table["SNind"] == -2)*(observation_table["instr"] == "WFI"))
    CCSNe = sum(observation_table["exptime"][inds])

    inds = where((observation_table["SNind"] == -1)*(observation_table["instr"] == "WFI"))
    imaging_only = sum(observation_table["exptime"][inds])

    inds = where((observation_table["SNind"] >= 0)*(observation_table["instr"] == "WFI"))
    SNeIa = sum(observation_table["exptime"][inds])

    total_exp_time = CCSNe*2. + SNeIa*2. + imaging_only
    print "total_exp_time ", total_exp_time
    total_survey_time = SN_data["survey_parameters"]["total_survey_time"]*365.24*86400
    return total_exp_time/total_survey_time

    



def get_FoM_from_file(item):
    f = pyfits.open(item)
    cmat = f[0].data
    f.close()

    z_list = concatenate(([0.05], arange(len(cmat) - 1)*0.05 + 0.125))
    print z_list

    FoM = get_FoM(cmat, z_list, adddesi = 0)[0]

    J = zeros([len(cmat), 2], dtype=float64) + 1.
    J[0,0] = 0.

    agg_prec_w1 = sum(linalg.inv(cmat[1:,1:]))
    agg_prec_w2 = dot(transpose(J), dot(linalg.inv(cmat), J))

    agg_prec_1 = 1./sqrt(agg_prec_w1)
    agg_prec_2 = sqrt(linalg.inv(agg_prec_w2)[0,0])
    
    print "Aggregate Precision Mag: ", agg_prec_1, agg_prec_2
    agg_prec_1 *= log(10.)/5.
    agg_prec_2 *= log(10.)/5.
    print "Aggregate Precision Dist: ", agg_prec_1, agg_prec_2

    print "FoM", item, FoM
    return FoM, agg_prec_1, agg_prec_2



def get_line_to_write(pic_cmat):
    pic, cmat = pic_cmat
    print pic, cmat

    SN_data = pickle.load(open(pic, 'rb'))

    #print "%s  %.3f  SN_data["total_time_used"]/(SN_data["survey_parameters"]["total_survey_time"]*365.24*86400)

    SNe_at_15, SNe_at_20 = get_SNe_with_imaging_LCs(SN_data)

    line_to_write = ""
    line_to_write += pic.split("/")[-2] + ","

    eff_area_fl = SN_data["survey_parameters"]["IFU_effective_area"]
    line_to_write += "6:"*(eff_area_fl.count("_16")) + "5:"*(eff_area_fl.count("_15")) + eff_area_fl.split("_")[3] + ","


    IFC_pixel_scale = SN_data["survey_parameters"]["IFU_pixel_scale"]
    line_to_write += "%.3f," % IFC_pixel_scale


    line_to_write += "%.3f," % (SN_data["total_time_used"]/86400.)
    line_to_write += "%s," % (str(SN_data["survey_parameters"]["tier_parameters"]["tier_fraction_time"]).replace('[', '').replace(']', '').replace(' ', '').replace(',', '+'))
    line_to_write += "+".join([item[0] for item in SN_data["survey_parameters"]["spectra_depth_and_phase"]]) + ","
    line_to_write += "%.2f+%.2f+%.2f+%.2f," % (SN_data["survey_parameters"]["shallow_SNR"], SN_data["survey_parameters"]["medium_SNR"], SN_data["survey_parameters"]["deep_SNR"], SN_data["survey_parameters"]["reference_SNR"])

    current_place = 0
    while len(SN_data["SN_observations"][current_place]["filts"]) == 0:
        current_place += 1

    lens = [len(item) for item in SN_data["SN_observations"][current_place]["filts"]]
    print lens

    line_to_write += str(min(lens) == 1) + ","

    wide_square_degrees = 0.
    deep_square_degrees = 0.
    wide_exp_time_per_filter = 0
    deep_exp_time_per_filter = 0
    wide_dithers_per_filter = 0
    deep_dithers_per_filter = 0

    for i in range(len(SN_data["survey_parameters"]["tier_parameters"]["square_degrees"])):
        trigger_redshifts = SN_data["survey_parameters"]["tier_parameters"]["trigger_redshifts"]
        trigger_fraction = SN_data["survey_parameters"]["tier_parameters"]["trigger_fraction"]

        highest_nonzero = where(array(trigger_fraction[i]) > 0)[0][-1]
        print trigger_redshifts[i], trigger_fraction[i], highest_nonzero

        tmp_exp_time = 0
        tmp_dithers = 0

        for j in range(len(SN_data["survey_parameters"]["tier_parameters"]["filters"][i])):
            if len(SN_data["survey_parameters"]["tier_parameters"]["filters"][i][j]) > 1:
                # If WFC filter:
                tmp_exp_time = SN_data["survey_parameters"]["tier_parameters"]["exp_times_per_dither"][i][j]
                tmp_dithers = SN_data["survey_parameters"]["tier_parameters"]["dithers_per_filter"][i][j]

        if trigger_redshifts[i][highest_nonzero] < 1:
            wide_square_degrees += SN_data["survey_parameters"]["tier_parameters"]["square_degrees"][i]
            wide_exp_time_per_filter = tmp_exp_time
            wide_dithers_per_filter = tmp_dithers
        else:
            deep_square_degrees += SN_data["survey_parameters"]["tier_parameters"]["square_degrees"][i]
            deep_exp_time_per_filter = tmp_exp_time
            deep_dithers_per_filter = tmp_dithers

    wide_exp_time_per_filter = int(around(wide_exp_time_per_filter))
    deep_exp_time_per_filter = int(around(deep_exp_time_per_filter))


    line_to_write += str(wide_square_degrees) + "," + str(wide_exp_time_per_filter) + "x" + str(wide_dithers_per_filter) + "," + str(deep_square_degrees) + "," + str(deep_exp_time_per_filter) + "x" + str(deep_dithers_per_filter) + ","

    has_IFS_mask = array([len(SN_data["SN_observations"][i]["IFS_dates"]) > 0 for i in range(SN_data["nsne"])])


    for i in range(len(zbins) - 1):
        counts = sum(has_IFS_mask*(SN_data["SN_table"]["redshifts"] > zbins[i])*(SN_data["SN_table"]["redshifts"] <= zbins[i+1]))
        line_to_write += "%i," % counts

    line_to_write += "%i,%i," % (SNe_at_15, SNe_at_20)


    if FoM_table:
        FoM, agg_prec_1, agg_prec_2 = get_FoM_from_file(cmat)
        line_to_write += "%s,%.1f,%.2g,%.2g" % (cmat.split("/")[-2], FoM, agg_prec_1, agg_prec_2)
    survey_eff = get_survey_efficiency(SN_data)
    line_to_write += ",%.2f" % survey_eff
    
    return line_to_write



FoM_table = int(sys.argv[2])
cores_to_use = 20

if FoM_table:
    zbins = arange(0.1, 2.01, 0.1)
else:
    zbins = [0, 0.4, 0.8, 1.1, 1.4, 1.7, 2.0]
 
zbins_txt = ""
for i in range(len(zbins) - 1):
    zbins_txt += ",%.1f<z<%.1f" % (zbins[i], zbins[i+1])

if FoM_table:
    f = open("FoM_table.csv", 'w')
    f.write("Survey,Cycle,Pixel Scale,Time Used,Tier Fraction,Spectra Type,Spectra S/N,Has Ground,Wide Square Degrees,Wide Exp Time per Filter,Deep Square Degrees,Deep Exp Time per Filter" + zbins_txt + ",SNe at S/N 15,SNe at S/N 20,FoM Params,FoM,Aggregate Precision 1,Aggregate Precision 2,Survey Efficiency\n")
else:
    f = open("summary.csv", 'w')
    f.write("Survey,Cycle,Pixel Scale,Time Used,Tier Fraction,Spectra Type,Spectra S/N,Has Ground,Wide Square Degrees,Wide Exp Time per Filter,Deep Square Degrees,Deep Exp Time per Filter" + zbins_txt + ",SNe at S/N 15,SNe at S/N 20,Survey Efficiency\n")

pic_cmats = []
for pic in glob.glob(sys.argv[1] + "/*/*pickle*"):
    print pic

    if FoM_table:
        for cmat in glob.glob(pic[:pic.rfind("/")] + "/*/cmat.fits")*FoM_table + [None]*(1 - FoM_table):
            print cmat
            pic_cmats.append((pic, cmat))
    else:
        pic_cmats.append((pic, None))

if len(pic_cmats) == 0:
    print "Couldn't match ", sys.argv[1], ", no cmats found!"

pool = mp.Pool(processes = cores_to_use)
lines_to_write = pool.map(get_line_to_write, pic_cmats)

f.write('\n'.join(lines_to_write))

f.close()



