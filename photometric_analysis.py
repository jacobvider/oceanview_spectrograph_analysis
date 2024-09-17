import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simpson, trapezoid, trapz
import glob
from scipy import ndimage, datasets
from IPython.display import Image, display
import xml.etree.ElementTree as ET
import astropy.units as u

# # Import reference transmission curves for chroma filters (B-Bessell, V-Bessell, R-Sloan)
def calc_transmission(filter, ref, dark):
    subtracted_filters = []
    subtracted_refs = []
    avg_filter = np.zeros_like(np.loadtxt(filter[0], skiprows=14)[:, 1])
    for i in range(len(filter)):
        subtracted_filter = np.loadtxt(filter[i], skiprows=14)[:, 1] - dark
        subtracted_filters.append(subtracted_filter)
        avg_filter += subtracted_filter
    avg_filter /= len(subtracted_filters)
    avg_fmed = ndimage.median_filter(avg_filter, size=40)

    avg_ref = np.zeros_like(np.loadtxt(ref[0], skiprows=14)[:, 1])
    for i in range(len(ref)):
        subtracted_ref = np.loadtxt(ref[i], skiprows=14)[:, 1] - dark
        subtracted_refs.append(subtracted_ref)
        avg_ref += subtracted_ref
    avg_ref /= len(subtracted_refs)
    avg_rmed = ndimage.median_filter(avg_ref, size=120)

    trans_med = [avg_fmed[i] / avg_rmed[i] for i in range(len(wv))]

    return trans_med


def calc_transmission2(filter, ref, dark):
    subtracted_filters = []
    subtracted_refs = []
    avg_filter = np.zeros_like(np.loadtxt(filter[0], skiprows=14)[:, 1])
    for i in range(len(filter)):
        subtracted_filter = np.loadtxt(filter[i], skiprows=14)[:, 1] - dark
        subtracted_filters.append(subtracted_filter)
        avg_filter += subtracted_filter
    avg_filter /= len(subtracted_filters)
    avg_fmed = ndimage.median_filter(avg_filter, size=60)
    avg_ref = np.zeros_like(np.loadtxt(ref[0], skiprows=14)[:, 1])
    for i in range(len(ref)):
        subtracted_ref = np.loadtxt(ref[i], skiprows=14)[:, 1] - dark
        subtracted_refs.append(subtracted_ref)
        avg_ref += subtracted_ref
    avg_ref /= len(subtracted_refs)
    avg_rmed = ndimage.median_filter(avg_ref, size=200)

    avg_rmed[avg_rmed == 0] = 1e-10 #new line

    trans_med = [avg_fmed[i] / avg_rmed[i] for i in range(len(avg_fmed))]
    return trans_med, wv_AA

# Import reference transmission curves for chroma filters (B-Bessell, V-Bessell, R-Sloan)
def calc_chroma(file):
    data = np.loadtxt(file)
    wavelengths = data[:, 0] * u.nm
    transmissions = data[:, 1]
    return wavelengths, transmissions

# Import reference transmission curves for chroma filters (B-Bessell, V-Bessell, R-Sloan)
def calc_svo(input):
    tree_u = ET.parse(input)  
    root_u = tree_u.getroot()
    wavelengths_u = [] 
    transmissions_u = []

    for tr in root_u.findall('.//TR'):
        tds = tr.findall('TD')
        wavelength = float(tds[0].text) 
        transmission = float(tds[1].text)
        wavelengths_u.append(wavelength)
        transmissions_u.append(transmission)

    wavelengths_u = np.array(wavelengths_u)
    transmissions_u = np.array(transmissions_u)
    return transmissions_u, wavelengths_u

#calculate average transmission curve for two filters tested in the lab
def calc_avg_transmission(filter1, ref1, dark1, filter2, ref2, dark2):
    trans_med1 = calc_transmission(filter1, ref1, dark1)
    trans_med2 = calc_transmission(filter2, ref2, dark2)
    avg_trans = [(trans_med1[i] + trans_med2[i])/2 for i in range(len(wv))]
    return avg_trans

#plot transmission curve for lab data
def plot_transmission(filter, ref, dark, wv1, wv2, label, color, ls):
    wv_AA_trans = wv.to(u.AA)
    trans_med = calc_transmission(filter, ref, dark)
    index1 = np.abs(wv_AA_trans.value - wv1).argmin()
    index2 = np.abs(wv_AA_trans.value - wv2).argmin()
    plt.plot(wv_AA_trans[index1:index2], trans_med[index1:index2], label=label, color=color, linestyle=ls)

#plot reference transmission curves for chroma filters
def plot_chroma_curves(file, wv1, wv2, label, color, ls):
    wavelengths, transmissions = calc_chroma(file)
    wv_AA_chroma = wavelengths.to(u.AA)
    index1 = np.abs(wv_AA_chroma.value - wv1).argmin()
    index2 = np.abs(wv_AA_chroma.value - wv2).argmin()
    plt.plot(wv_AA_chroma[index1:index2], transmissions[index1:index2], color = color, label = label, linestyle = ls)

#plot reference transmission curve for svo data
def plot_svo(input, wv1, wv2, color, label, ls):
    transmissions_u, wavelengths_u = calc_svo(input)
    wavelengths_u = np.array(wavelengths_u)
    transmissions_u = np.array(transmissions_u)
    index1 = np.abs(wavelengths_u - wv1).argmin()
    index2 = np.abs(wavelengths_u - wv2).argmin()
    plt.plot(wavelengths_u[index1:index2], transmissions_u[index1:index2], color = color, label= label, linestyle = ls)

#plot average transmission curve for two filters tested in the lab
def plot_avg_trans(filter1, ref1, dark1, filter2, ref2, dark2, wv1, wv2, label, color, ls):
    wv_AA_avg = wv.to(u.AA)
    avg = calc_avg_transmission(filter1, ref1, dark1, filter2, ref2, dark2)
    index1 = np.abs(wv_AA_avg.value - wv1).argmin()
    index2 = np.abs(wv_AA_avg.value - wv2).argmin()
    plt.plot(wv_AA_avg[index1:index2], avg[index1:index2], label=label, color=color, linestyle=ls)

#calculate central wavelength values for lab data
def central_wavelength_weighted(filter, ref, dark, wv1, wv2):
    trans = calc_transmission(filter, ref, dark)
    index1 = np.abs(wv_AA.value - wv1).argmin()
    index2 = np.abs(wv_AA.value - wv2).argmin()
    transmissions = trans[index1:index2]
    wavelengths = wv_AA[index1:index2]
    wavelengths_arr = np.array(wavelengths)
    transmissions_arr = np.array(transmissions)
    transmissions_normalized = transmissions_arr / np.sum(transmissions_arr)
    weighted_avg = np.sum(wavelengths_arr * transmissions_normalized)
    return weighted_avg

#calculate central wavelength value for chroma data
def central_wavelength_chroma(input, wv1, wv2):
    wavelengths, transmissions = calc_chroma(input)
    wavelengths_interp = np.linspace(wv1, wv2, num=1000)  # Increase num for finer resolution
    transmissions_interp = np.interp(wavelengths_interp, wavelengths.value, transmissions)
    transmissions_normalized = transmissions_interp / np.sum(transmissions_interp)
    weighted_avg = np.sum(wavelengths_interp * transmissions_normalized)
    return weighted_avg

#read in atmospheric data
def read_atmo_data(data):
    with open(data, 'r') as file:
        file.readline()
        file.readline()
        atmo = np.loadtxt(file, dtype=float)
    atmo[:, 1] = 10 ** (-0.4 * atmo[:, 1])
    wv_atmo = atmo[:, 0][:59]
    trans_atmo = atmo[:, 1][:59]
    return wv_atmo, trans_atmo


def read_pixis_qe(wv_pixis, qe):
    pixis_wv_AA = wv_pixis.to(u.AA)  # Convert to Angstroms
    return pixis_wv_AA, qe

#get spectral data for SN 1992a
def get_spectra_92a(input):
        transmissions_u, wavelengths_u = calc_svo('Swift.UVOT.U_trn.xml')
        spec_92A = np.loadtxt(input, skiprows=1)
        normresp = transmissions_u / np.max(transmissions_u)

        spec_92A_b = np.interp(wavelengths_u, spec_92A[:, 0], spec_92A[:, 1])

        spec_92A_suvot = spec_92A_b * normresp * 1e14
        return spec_92A_suvot

#combine spectral data for SN 1992a, pixis qe data, atmospheric data, and data for filter tested in the lab
def plot_interp_lab(filter, ref, dark, color, label, observation_site, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        trans_uvot = calc_transmission(filter, ref, dark)
        # transmissions_u, wavelengths_u = calc_svo(input)
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp = trans_uvot / np.max(trans_uvot)
        normresp_interp = np.interp(wv_val, wv_val, normresp)
        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        wv_atmo, trans_atmo = read_atmo_data(observation_site)
        apo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
        product1 = (pixis_interp * apo_interp) * (spec_92A_b)*1e14*normresp_interp
        plt.plot(wv_AA[773:1543], product1[773:1543], color = color, label = label)
        # plt.legend()

def plot_interp_lab2(filter, ref, dark, color, label, ls, observation_site, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        trans_uvot, wv_uvot = calc_transmission2(filter, ref, dark)
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp = trans_uvot / np.max(trans_uvot)
        normresp_interp = np.interp(wv_val, wv_uvot.value, normresp)
        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        wv_atmo, trans_atmo = read_atmo_data(observation_site)
        apo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
        product1 = (pixis_interp * apo_interp) * (spec_92A_b)*1e14*normresp_interp

        wavelength_range = wv_AA[773:1543]
        product1_range = product1[773:1543]
        product2_range = product1[773:1079]

        area = trapz(product1_range, wavelength_range)
        area_3400 = trapz(product2_range, wv_AA[773:1079])
        percent_3400 = (area_3400/area)*100
        percent_3400_rounded = round(percent_3400.value, 4)
        print(f"{percent_3400_rounded}%")

        plt.plot(wv_AA[773:1543], product1[773:1543], color = color, label = label, linestyle = ls)

def plot_interp_chroma(filter, ref, dark, color, label, observation_site, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        trans_uvot = calc_transmission(filter, ref, dark)
        # transmissions_u, wavelengths_u = calc_svo(input)
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp = trans_uvot / np.max(trans_uvot)
        normresp_interp = np.interp(wv_val, wv_val, normresp)
        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        wv_atmo, trans_atmo = read_atmo_data(observation_site)
        apo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
        product1 = (pixis_interp * apo_interp) * (spec_92A_b)*1e14*normresp_interp
        plt.plot(wv_AA[773:1543], product1[773:1543], color = color, label = label)
        
#combine spectral data for SN 1992a, pixis qe data, atmospheric data, and reference data for chroma filter
def plot_interp_chroma(input, color, label, observation_site, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        transmissions_u, wavelengths_u  = calc_chroma(input) 
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp_svo = transmissions_u / np.max(transmissions_u)
        normresp_svo_interp = np.interp(wv_val, wavelengths_u, normresp_svo)
        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        wv_atmo, trans_atmo = read_atmo_data(observation_site)
        atmo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
        product = (pixis_interp * atmo_interp) * (spec_92A_b)*1e14*normresp_svo_interp
        plt.plot(wv_AA[773:1420], product[773:1420], color = color, label = label)
        # plt.legend()

#combine spectral data for SN 1992a, pixis qe data, atmospheric data, and reference data for chroma filter
def plot_interp_svo(input, color, label, observation_site, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        transmissions_u, wavelengths_u = calc_svo(input)
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp_svo = transmissions_u / np.max(transmissions_u)
        normresp_svo_interp = np.interp(wv_val, wavelengths_u, normresp_svo)

        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        svo_interp = np.interp(wv_val, wavelengths_u, transmissions_u)
        wv_atmo, trans_atmo = read_atmo_data(observation_site)
        apo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
        # flux_interp = np.interp(wv_val, wv_92a, spec_92a)
        product1 = (pixis_interp * apo_interp * svo_interp) * (spec_92A_b)*1e14*normresp_svo_interp

        plt.plot(wv_AA[773:1543], product1[773:1543], color = color, label = label)
        # plt.legend()

def plot_interp_svo_analyze(input, color, label, observation_site, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        transmissions_u, wavelengths_u = calc_svo(input)
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp_svo = transmissions_u / np.max(transmissions_u)
        normresp_svo_interp = np.interp(wv_val, wavelengths_u, normresp_svo)
        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        svo_interp = np.interp(wv_val, wavelengths_u, transmissions_u)
        wv_atmo, trans_atmo = read_atmo_data(observation_site)
        apo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
        # flux_interp = np.interp(wv_val, wv_92a, spec_92a)
        product1 = (pixis_interp * apo_interp * svo_interp) * (spec_92A_b)*1e14*normresp_svo_interp
        wavelength_range = wv_AA[773:1543]
        product1_range = product1[773:1543]
        product2_range = product1[773:1079]
        area = trapz(product1_range, wavelength_range)
        area_3400 = trapz(product2_range, wv_AA[773:1079])
        # area_3200 = trapz(product3_range, wv_AA[773:926])
        percent_3400 = (area_3400/area)*100
        percent_3400_rounded = round(percent_3400.value, 4)
        # percent_3200 = (area_3200/area)*100
        print(f"{percent_3400_rounded}%")
        plt.plot(wv_AA[773:1543], product1[773:1543], color=color, label=label)
        # plt.legend()

def plot_ground_instr_resp(color, label, observation_site, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        # spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        wv_atmo, trans_atmo = read_atmo_data(observation_site)
        apo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
        # flux_interp = np.interp(wv_val, wv_92a, spec_92a)
        product1 = (pixis_interp * apo_interp)
        # percent_3200 = (area_3200/area)*100
        plt.plot(wv_AA[773:1543], product1[773:1543], color=color, label=label)
        # plt.legend()

def plot_interp_lab_analyze(filter, ref, dark, color, label, ls, observation_site, coating, coating_qe):
    # Read data
    pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
    wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
    pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
    trans_uvot = calc_transmission(filter, ref, dark)
    spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
    normresp = trans_uvot / np.max(trans_uvot)
    normresp_interp = np.interp(wv_val, wv_val, normresp)
    spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
    pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
    wv_atmo, trans_atmo = read_atmo_data(observation_site)
    apo_interp = np.interp(wv_val, wv_atmo, trans_atmo)
    product1 = (pixis_interp * apo_interp) * (spec_92A_b) * 1e14 * normresp_interp
    
    # Compute area under the curve
    wavelength_range = wv_AA[773:1543]
    product1_range = product1[773:1543]
    product2_range = product1[773:1079]

    area = trapz(product1_range, wavelength_range)
    area_3400 = trapz(product2_range, wv_AA[773:1079])
    percent_3400 = (area_3400/area)*100
    percent_3400_rounded = round(percent_3400.value, 4)

    plt.plot(wv_AA[773:1543], product1[773:1543], color=color, label=label, linestyle=ls)
    print(f"{percent_3400_rounded}%")

    # plt.legend()


def plot_interp_lab_analyze_noatmo(filter, ref, dark, color, label, coating, coating_qe):
    pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
    wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
    pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
    trans_uvot = calc_transmission(filter, ref, dark)
    spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
    normresp = trans_uvot / np.max(trans_uvot)
    normresp_interp = np.interp(wv_val, wv_val, normresp)
    spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
    pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
    product1 = (pixis_interp) * (spec_92A_b) * 1e14 * normresp_interp
    
    plt.plot(wv_AA[773:1543], product1[773:1543], color=color, label=label)

    # plt.legend()


def plot_interp_svo_analyze_noatmo(input, color, label, coating, coating_qe):
        pixis_wv_AA, pixis_qe = read_pixis_qe(coating, coating_qe)
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        transmissions_u, wavelengths_u = calc_svo(input)
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp_svo = transmissions_u / np.max(transmissions_u)
        normresp_svo_interp = np.interp(wv_val, wavelengths_u, normresp_svo)
        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        pixis_interp = np.interp(wv_val, pixis_wv_val[:6], pixis_qe[:6])
        svo_interp = np.interp(wv_val, wavelengths_u, transmissions_u)
        # flux_interp = np.interp(wv_val, wv_92a, spec_92a)
        product1 = (pixis_interp * svo_interp) * (spec_92A_b)*1e14*normresp_svo_interp

        plt.plot(wv_AA[773:1543], product1[773:1543], color=color, label=label)
        # plt.legend()

def plot_svo_interp_trn(input, color, label):
        wv_val = [wv_AA[i].value for i in range(len(wv_AA))]
        # pixis_wv_val = [pixis_wv_AA[i].value for i in range(len(pixis_wv_AA))]
        transmissions_u, wavelengths_u = calc_svo(input)
        spec_92A = np.loadtxt('sn1992a-19920124.214-hst.flm', skiprows=1)
        normresp_svo = transmissions_u / np.max(transmissions_u)
        normresp_svo_interp = np.interp(wv_val, wavelengths_u, normresp_svo)

        spec_92A_b = np.interp(wv_val, spec_92A[:, 0], spec_92A[:, 1])
        svo_interp = np.interp(wv_val, wavelengths_u, transmissions_u)
        product1 = (svo_interp) * (spec_92A_b)*1e14*normresp_svo_interp

        wavelength_range = wv_AA[773:1543]
        product1_range = product1[773:1543]
        product2_range = product1[773:1079]

        area = trapz(product1_range, wavelength_range)
        area_3400 = trapz(product2_range, wv_AA[773:1079])
        percent_3400 = (area_3400/area)*100
        percent_3400_rounded = round(percent_3400.value, 4)

        plt.plot(wv_AA[773:1543], product1[773:1543], color = color, label = label)
        print(f"{percent_3400_rounded}%")


def subtract_trans(svo_input, filter, ref, dark):
    trans_med, wv_med = calc_transmission2(filter, ref, dark)
    trans_u, wavelengths_u = calc_svo(svo_input) #'Swift.UVOT.U_fil.xml'
    trans_svo_interp = np.interp(wv_med.value, wavelengths_u, trans_u)
    trans_subtracted = np.abs(np.array(trans_med) - trans_svo_interp)
    return trans_subtracted, wv_med


def calc_fwhm_cwl_chroma(input, wv1, wv2):
    wavelengths, transmissions = calc_chroma(input)
    wavelengths_AA = wavelengths.to(u.AA)
    wavelengths_interp = np.linspace(wv1, wv2, num=1000)
    transmissions_interp_cwl = np.interp(wavelengths_interp, wavelengths.value, transmissions) #this line
    transmissions_normalized_cwl = transmissions_interp_cwl / np.sum(transmissions_interp_cwl)
    weighted_avg = np.sum(wavelengths_interp * transmissions_normalized_cwl)
    transmissions_interp = np.interp(wavelengths_interp, wavelengths_AA.value, transmissions) #this line
    amp_y = np.max(transmissions_interp)
    halfamp_y = amp_y * 0.5
    intersect_indices = np.where(np.abs(transmissions_interp - halfamp_y) < 0.05)[0]
    intersect_wavelengths = wavelengths_interp[intersect_indices]
    # intersect_transmissions = transmissions_interp[intersect_indices]
    wv_max = np.max(intersect_wavelengths)
    wv_min = np.min(intersect_wavelengths)
    fwhm = wv_max - wv_min
    # fwhm = intersect_wavelengths[1] - intersect_wavelengths[0]
    print(fwhm, weighted_avg)


def calc_fwhm_cwl(filter, ref, dark, wv1, wv2):
    wv_AA_fwhm = wv.to(u.AA)
    trans = calc_transmission(filter, ref, dark)
    index1 = np.abs(wv_AA_fwhm.value - wv1).argmin()
    index2 = np.abs(wv_AA_fwhm.value - wv2).argmin()
    wv_arr = wv_AA_fwhm[index1: index2]
    trans_arr = trans[index1:index2]

    amp_y = np.max(trans_arr)
    amp_x = np.argmax(trans_arr)
    halfamp_y = amp_y * 0.5
    intersections = []
    for i in range(1, len(trans_arr)):
        if (trans_arr[i-1] < halfamp_y and trans_arr[i] > halfamp_y) or (trans_arr[i-1] > halfamp_y and trans_arr[i] < halfamp_y):
            x1, x2 = wv_arr[i-1], wv_arr[i]
            y1, y2 = trans_arr[i-1], trans_arr[i]
            
            x_intersect = x1 + (halfamp_y - y1) * (x2 - x1) / (y2 - y1)
            intersections.append(x_intersect)
    if len(intersections) == 1:
        intersections.append((wv2) * u.AA)
        fwhm = intersections[1] - intersections[0]
        cwl = fwhm/2 + intersections[0]
    else:
        fwhm = intersections[1] - intersections[0]
        cwl = fwhm/2 + intersections[0]
    return fwhm, cwl