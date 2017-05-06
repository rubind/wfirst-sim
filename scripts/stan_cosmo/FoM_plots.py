from numpy import *
import glob
import matplotlib.pyplot as plt
import pyfits
import os
import sys
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/cosmo/")
from astro_functions import get_FoM

def get_FoM_from_file(item):
    f = pyfits.open(item)
    cmat = f[0].data
    f.close()

    z_list = concatenate(([0.05], arange(len(cmat) - 1)*0.05 + 0.125))
    print z_list

    FoM = get_FoM(cmat, z_list, adddesi = 0)[0]

    print "FoM", item, FoM
    return FoM

def get_survey_type(survey_name):
    #survey_00001,5:150518,194.652,shallow+medium+deep,12.00+20.60+34.30+40.00,False,9.0,50x2,3.0,330x2,207,1167,191,214,169,0,
    f = open("summary.csv")
    lines = f.read().split('\n')
    f.close()

    for line in lines:
        parsed = line.split(",")
        if len(parsed) > 3:
            if parsed[0] == survey_name:
                print "parsed ", parsed, parsed[12:17]
                has_ground = eval(parsed[7])
                

                if sum([float(item) for item in parsed[14:18]]) == 0:
                    max_z = "$z=0.8$"
                    color = 'r'
                elif sum([float(item) for item in parsed[15:17]]) == 0:
                    max_z = "$z=1.1$"
                    color = 'orange'
                elif sum([float(item) for item in parsed[16:18]]) == 0:
                    max_z = "$z=1.4$"
                    color = 'green'
                elif sum([float(item) for item in parsed[17:18]]) == 0:
                    max_z = "$z=1.7$"
                    if has_ground:
                        color = 'blue'
                    else:
                        color = 'cyan'
                else:# float(parsed[16]) == 0:
                    max_z = "$z=2.0$"
                    color = (0.5, 0, 1.)
                
                

                cycle = parsed[1][0]
                nsne = "+".join(parsed[10:16])

                if has_ground and max_z == "$z=0.8$":
                    return "Ground Discovery to " + max_z, color, cycle, nsne, has_ground
                elif has_ground:
                    return "Ground and Space Discovery to " + max_z, color, cycle, nsne, has_ground
                else:
                    return "Space Discovery to " + max_z, color, cycle, nsne, has_ground
                
    print "Couldn't find ", survey_name

def cmat_to_vals(cmat):
    parsed = cmat.split("/")[-2].split("_")

    #['FoM', 'IPC=0.02', 'nredcoeff=2', 'fundcal=0.005', 'crnl=0.005', 'include', 'sys=1', 'read', 'noise', 'floor=4', 'graydisp=0.08', 'dark', 'current=0.01', 'TTel=282']
    #assert len(parsed) == 14, parsed
    #assert parsed[0] == "FoM"
    #assert parsed[11] == "dark"
    
    vals = {}
    for item in parsed:
        if item.count("="):
            vals[item.split("=")[0]] = float(item.split("=")[1])
    return vals

def how_diff(vals1, vals2):
    for key in vals1:
        if not isclose(vals1[key], vals2[key]):
            return [key]
    return vals1.keys()
        

plot_type = "Saul"#"paper" # paper, Saul, or ""
#cmat_list = glob.glob("survey_tests_19/*/*/cmat.fits")
cmat_list = glob.glob("survey_00155/*/cmat.fits") + glob.glob("survey_space_low/*/cmat.fits") + glob.glob("survey_lsst_low/*/cmat.fits") + glob.glob("survey_lsst_high/*/cmat.fits") + glob.glob("survey_tests_new/survey_00181/*/cmat.fits")
cmat_list.sort()


survey_list = [item.split("/")[-3] for item in cmat_list]
#survey_colors = ['b', 'm', 'c', 'g', 'r', 'orange', 'k', 'brown', 'blueviolet', 'gold', (0.8, 1, 0), 'gray', (0.6, 0.3, 0.3)]
survey_colors = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"] # http://stackoverflow.com/questions/2328339/how-to-generate-n-different-colors-for-any-natural-number-n

key_names = {"crnl": "Count-Rate Nonlinearity Uncertainty (mag)",
             "IPC": "Inter-Pixel Capacitance",
             "TTel": "Telescope Temperature (K)",
             "current": "Dark Current (e-/s/pix)",
             "floor": "Read-Noise Floor (e-)",
             "fundcal": "Fundamental Calibration Uncertainty (mag)",
             "graydisp": "Gray Dispersion (mag)",
             "MWZP": "Milky-Way Extinction ZP Uncertainty (mag)",
             "redwave": "Reddest Rest-Frame Wavelength (A)",
             "MWnorm": "Milky-Way Extinction Normalization Uncertainty",
             "nredcoeff": "Redshift Coefficients (1 = const, 2 = slope, 3 = curve)",
             "sys": "Include Systematics (calibration+mw extinction+intergalactic extinction)"}

survey_set = list(unique(survey_list))
assert len(survey_colors) >= len(survey_set)

print len(cmat_list), "found"
print "survey_set ", survey_set


attributes = {}

vals = cmat_to_vals(cmat_list[0])
print "vals ", vals

for key in vals:
    attributes[key] = []

key_list = vals.keys()

for cmat in cmat_list:
    vals = cmat_to_vals(cmat)
    for key in key_list:
        attributes[key].append(vals[key])

print "attributes ", attributes
base_case = {}
for key in key_list:
    base_case[key] = median(attributes[key])

print "base_case ", base_case




for key in key_list:
    labeled_items = []


    for survey in survey_set:
        xvals = []
        yvals = []

        survey_ind = survey_set.index(survey)
        for i in range(len(cmat_list)):
            
            vals = cmat_to_vals(cmat_list[i])
            diff_keys = how_diff(base_case, vals)
            #print cmat_list[i], diff_key
            
            if diff_keys.count(key):
                if survey_list[i] == survey:
                    xvals.append(vals[key])
                    yvals.append(get_FoM_from_file(cmat_list[i]))

        label, color, cycle, nsne, has_ground = get_survey_type(survey)
        if plot_type == "Saul":
            if len(xvals) > 0:
                plt.plot(xvals, yvals, 'o-', color = survey_colors[survey_ind])
                plt.text(xvals[0], yvals[0], survey + "_ground"*has_ground + "_cycle=" + cycle + "_nsne=" + nsne  + "  ", size = 7, color = survey_colors[survey_ind], horizontalalignment = random.choice(['right', 'left']))
        else:
            plt.plot(xvals, yvals, 'o-', color = color, label = label*(labeled_items.count(label) == 0))
            labeled_items.append(label)

    plt.legend(loc = 'best', fontsize = 12 - 6*(plot_type == ""))
    if not (plot_type == "Saul"):
        plt.ylim(0, plt.ylim()[1])
    if plot_type == "Saul":
        xlim = plt.xlim()
        plt.xlim(xlim[0]*2.5 - 1.5*xlim[1], -1.5*xlim[0] + 2.5*xlim[1])
    plt.xlabel(key_names[key])
    plt.ylabel("FoM")
    #if plot_type == "paper":
    #    plt.yticks([])
    plt.savefig("FoM_key=%s%s.%s" % (key, "_small"*(len(survey_set) < 5), "pdf"*(plot_type != "paper") + "eps"*(plot_type == "paper")), bbox_inches = 'tight')
    plt.close()


