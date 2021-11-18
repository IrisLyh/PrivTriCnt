import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy.core.fromnumeric import prod
from numpy.core.records import array
import xlrd
from scipy.interpolate import make_interp_spline

table = xlrd.open_workbook("./results.xlsx")
#plot result
'''
methods = ["2Phase-RLDP","DDP","2Rounds-LDP","Basic-RLDP","Sample Error"]
colors = ["red","blue","black","green","deeppink"]
markers = ["o","s","x","*","+"]
datasets = ["result-wiki","result-cit","result-Enron"]

name = datasets[2]
result_sheet = table.sheet_by_name(name)
result_list = []
for i in range(1, result_sheet.ncols-1):              
    col_value = result_sheet.col_values(i)
    temp_list = []
    for j in range(1, len(col_value)):
        temp_list.append(col_value[j])
    result_list.append(temp_list)
result_array_wiki = np.array(result_list)

x=np.array([1,2,3,4,5,6,7,8,9,10])
for method in range(0,len(methods)):
    y=result_array_wiki[method]
    plt.plot(x,y,c=colors[method],marker=markers[method],label=methods[method])
    plt.yscale("log")
    plt.xticks(range(1,11,1))
    plt.xlabel("Privacy budget $\epsilon$")
    plt.ylabel("MRE")
    plt.legend(loc="upper right")
    fig_name = name+".png"
plt.savefig(fig_name)
'''

#plot sample rate
'''
colors = ["green","blue","red","black","deeppink"]
markers = ["o","x","s","*","+"]
datasets=["sample-rate-wiki","sample-rate-cit","sample-rate-enron"]
sample_rates = ["p=0.001","p=0.008","p=0.010","p=0.012","p=0.100"]
name = datasets[0]
result_sheet = table.sheet_by_name(name)
result_list = []
for i in range(1, result_sheet.ncols):              
    col_value = result_sheet.col_values(i)
    temp_list = []
    for j in range(1, len(col_value)):
        temp_list.append(col_value[j])
    result_list.append(temp_list)
result_array = np.array(result_list)

x=np.array([1,2,3,4,5,6,7,8,9,10])
for rate in range(0,len(sample_rates)):
    y=result_array[rate]
    plt.plot(x,y,c=colors[rate],marker=markers[rate],label=sample_rates[rate])
    plt.yscale("log")
    plt.xticks(range(1,11,1))
    plt.xlabel("Privacy budget $\epsilon$")
    plt.ylabel("MRE")
    plt.legend(loc="upper right")
    fig_name = name+".png"
plt.savefig(fig_name)
'''
#plot noisy local sensitivity
'''
colors = ["red","green","blue"]
markers = ["s","o","*"]
datasets = ["wiki-Vote","cit-HepTh","email-Enron"]
name = "upper_LS"
result_sheet = table.sheet_by_name(name)
result_list = []
for i in range(1, result_sheet.ncols):              
    col_value = result_sheet.col_values(i)
    temp_list = []
    for j in range(1, len(col_value)):
        temp_list.append(col_value[j])
    result_list.append(temp_list)
result_array = np.array(result_list)

x=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
for dataset in range(0,len(datasets)):
    y=result_array[dataset]
    plt.plot(x,y,c=colors[dataset],marker=markers[dataset],label=datasets[dataset])
    # plt.yscale("log")
    plt.xticks(np.arange(0.1,1.1,0.1))
    plt.xlabel("Privacy budget $\epsilon_1/2$")
    plt.ylabel("Average LS after subsampling")
    plt.legend(loc="upper right")
    # fig_name = name+".png"
plt.savefig("local-sensitivity.png")
'''

#plot k
print(mpl.get_configdir())
colors = ["red","green","blue"]
markers = ["s","o","*"]
datasets = ["$\epsilon_2=2.4$","$\epsilon_2=4.0$","$\epsilon_2=8.0$"]
linestyles = ["-","-.",":"]
name = "k"
result_sheet = table.sheet_by_name(name)
result_list = []
for i in range(1, result_sheet.ncols):              
    col_value = result_sheet.col_values(i)
    temp_list = []
    for j in range(1, len(col_value)):
        temp_list.append(col_value[j])
    result_list.append(temp_list)
result_array = np.array(result_list)

x=np.arange(1,5001,1)
# x=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
for dataset in range(0,len(datasets)):
    y=result_array[dataset]
    # x_new = np.linspace(1,5000,10000)
    # y_smooth = make_interp_spline(x,y,x_new)\
    fig=plt.figure(figsize=(8,6))
    plt.plot(x,y,c=colors[dataset],linestyle=linestyles[dataset],label=datasets[dataset])
    # plt.yscale("log")
    x_mid = x[int(len(x)/2)]
    font_legend = {'size':16,'weight':'normal'}
    plt.tick_params(labelsize=14)
    plt.xticks(np.arange(0,5001,1000))
    plt.xlabel("The maximum edge using frequency $k^*$",size=14,labelpad = 3)
    plt.ylabel("Privacy budget $\epsilon_2'$",size=14)
    plt.legend(loc="upper right",prop=font_legend)
    # fig_name = name+".png"
plt.savefig("k.png",pad_inches = 2)

