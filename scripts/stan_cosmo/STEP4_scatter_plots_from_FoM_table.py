from numpy import *
import matplotlib.pyplot as plt
import colorsys
from scipy.stats import gaussian_kde

def read_csv(csv):
    f = open(csv)
    lines = f.read().split('\n')
    f.close()

    data = {}
    headers = lines[0].split(",")
    for item in headers:
        data[item] = []

    for line in lines[1:]:
        parsed = line.split(",")

        if len(parsed) > 2:
            assert len(parsed) == len(headers), line
            
            for i in range(len(headers)):
                item = parsed[i]
                
                try:
                    item = float(item)
                except:
                    pass

                data[headers[i]].append(item)
    return data

def scatter_by_param(x, y, labels, key, shape_mask, FoM_run, dolabel = None):
    plt.figure()

    label_set = unique(labels)

    colors = [colorsys.hsv_to_rgb(item, 1, 1) for item in linspace(0, 2./3., len(label_set))]
    print colors


    for i in range(len(label_set)):
        inds = where(labels == label_set[i])
        for ind in inds[0]:
            label = ""
            if ind == inds[0][0]:
                label = label_set[i]
            plt.plot(x[ind], y[ind], '*'*int(shape_mask[ind]) + '.'*int(1 - shape_mask[ind]), label = label, color = colors[i], markeredgewidth = 0, markersize = 8)
            if dolabel != None:
                plt.text(x[ind], y[ind], dolabel[ind].replace("survey_", ""), color = colors[i], fontsize = 6)
                    
    plt.legend(loc = 'best', fontsize = 8)
    plt.savefig("FoM_scatter/FoM_%s_%s.pdf" % (key.replace("/", "-").replace(" ", "_"), FoM_run))
    plt.close()


def KDE_by_param(x, y, labels, key, shape_mask, FoM_run, xlim, ylim):
    plt.figure()

    label_set = unique(labels)

    colors = [colorsys.hsv_to_rgb(item, 1, 1) for item in linspace(0, 2./3., len(label_set))]
    lightcolors = [colorsys.hsv_to_rgb(item, 0.5, 1) for item in linspace(0, 2./3., len(label_set))]
    print "colors ", colors

    labels_used = []
    proxy = []

    for i in range(len(label_set)):
        inds = where(labels == label_set[i])

        if len(inds[0]) > 1:
            xs = x[inds]
            ys = y[inds]
            
            print "xs, ys", xs, ys
            
            X, Y = meshgrid(linspace(xs.min()*0.8, xs.max()*1.2, 100),
                            linspace(ys.min()*0.8, ys.max()*1.2, 100))
            
            xys = transpose(array(zip(xs, ys)))
            
            print "xys", xys
            
            kde_fn = gaussian_kde(xys)
            
            XYs = transpose(array(zip(reshape(X, 10000), reshape(Y, 10000))))
            
            Z = kde_fn(XYs)
            Z /= max(Z)
            
            plt.contourf(X, Y, reshape(Z, (100, 100)), levels = [Z.max()/4., Z.max()/2., Z.max()], colors = [lightcolors[i], colors[i], 'k'], alpha = 0.2)
            for ind in inds[0]:
                label = ""
                if ind == inds[0][0]:
                    label = label_set[i]
                plt.plot(x[ind], y[ind], '*'*shape_mask[ind] + '.'*(1 - shape_mask[ind]), color = colors[i], markeredgewidth = 0, markersize = 8, zorder = 3, label = label)

            print Z
            plt.xlabel("Number of SNe")
            plt.ylabel("DETF Figure of Merit")
            plt.title(key)
            
            #labels_used.append(label_set[i])
            #proxy.append( plt.Rectangle((0,0),1,1, fc = colors[i], ec = None, alpha = 0.3))
    plt.xlim(xlim)
    plt.ylim(ylim)


    #plt.legend(proxy, labels_used)
    plt.legend(loc = 'best')

    plt.savefig("FoM_scatter/FoM_KDE_%s_%s.pdf" % (key.replace("/", "-").replace(" ", "_"), FoM_run))
    plt.close()


def key_to_zmax(key):
    if key.count("<z<"):
        zmax = float(key.split("<")[-1])
        return zmax
    else:
        return None



def make_scatter(csv):
    data = read_csv(csv)
    data["nsne"] = zeros(len(data["1.8<z<1.9"]))
    for key in data:
        if key.count("<z<"):
            data["nsne"] += data[key]

    data["Max z"] = zeros(len(data["1.8<z<1.9"]), dtype=float32)
    for key in data:
        zmax = key_to_zmax(key)
        if zmax != None:
            for i in range(len(data["Max z"])):
                if data[key][i] > 0:
                    data["Max z"][i] = max(data["Max z"][i], zmax)

    data["1.0 < z < 1.4"] = zeros(len(data["1.8<z<1.9"]), dtype=float32)
    for key in ["1.0<z<1.1","1.1<z<1.2","1.2<z<1.3","1.3<z<1.4"]:
        data["1.0 < z < 1.4"] += data[key]
    data["1.0 < z < 1.4"] = [str(around(item/50.)*50.) for item in data["1.0 < z < 1.4"]]


    NSNe_ultimate = zeros(len(data["1.8<z<1.9"]), dtype=float32)
    NSNe_penultimate = zeros(len(data["1.8<z<1.9"]), dtype=float32)


    print data["Max z"]
    
    for key in data:
        zmax = key_to_zmax(key)
        if zmax != None:
            for i in range(len(data["Max z"])):
                if isclose(zmax, data["Max z"][i]):
                    NSNe_ultimate[i] = float(data[key][i])
                elif isclose(zmax, data["Max z"][i] - 0.2):
                    NSNe_penultimate[i] = float(data[key][i])

    assert all(NSNe_ultimate > 0), NSNe_ultimate
    assert all(NSNe_penultimate > 0), NSNe_penultimate
    
    data["Taper"] = around(NSNe_ultimate/NSNe_penultimate, 1)
    

    for key in data:
        data[key] = array(data[key])

    data["Has Ground"] = (data["Wide Square Degrees"] == 28.8) + (data["Wide Square Degrees"] > 99)
    print data["Has Ground"]

    FoM_runs = unique(data["FoM Params"])


    for FoM_run in FoM_runs:
        inds = where(data["FoM Params"] == FoM_run)
        ylim = (50*floor(min(data["FoM"][inds])/50.), 50*ceil(max(data["FoM"][inds])/50.))
        print "ylim", ylim

        for key in ["Wide Square Degrees", "Deep Square Degrees", "Spectra S/N", "Pixel Scale", "Max z", "Taper", "Tier Fraction", "1.0 < z < 1.4"]:
            inds = where(data["FoM Params"] == FoM_run)

            scatter_by_param(data["nsne"][inds], data["FoM"][inds], data[key][inds], key,
                             data["Has Ground"][inds], FoM_run)
            
            scatter_by_param(data["nsne"][inds], data["FoM"][inds], data[key][inds], key,
                             data["Has Ground"][inds], FoM_run + "_label", dolabel = data["Survey"])


            KDE_by_param(data["nsne"][inds], data["FoM"][inds], data[key][inds], key,
                             data["Has Ground"][inds], FoM_run + "_IFCSNe", xlim = [0, 11000], ylim = ylim)


            """
            KDE_by_param(data["SNe at S/N 40"][inds], data["FoM"][inds], data[key][inds], key,
                             data["Has Ground"][inds], FoM_run  + "_WFC40", xlim = [0, 40000], ylim = ylim)

            KDE_by_param(data["SNe at S/N 60"][inds], data["FoM"][inds], data[key][inds], key,
                             data["Has Ground"][inds], FoM_run  + "_WFC60", xlim = [0, 40000], ylim = ylim)
            """

make_scatter("FoM_table.csv")

