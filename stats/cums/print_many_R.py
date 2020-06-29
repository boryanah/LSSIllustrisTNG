import numpy as np
import matplotlib.pyplot as plt
import plotparams
plotparams.buba()

proxy = 'm200m'
opts = np.array(['shuff','partial_env_cw','partial_vani','partial_s2r','partial_tot_pot'])

bin_centers = np.load("data_many/rs.npy")
chosen = 2
print(bin_centers[chosen])

# if not ratios
second_fid = np.load("data_many/second_true_mean.npy")[chosen]
second_error_fid = np.load("data_many/second_true_err.npy")[chosen]
third_fid = np.load("data_many/third_true_mean.npy")[chosen]
third_error_fid = np.load("data_many/third_true_err.npy")[chosen]

#second_fid = np.load("data_many/second_rat_"+proxy+"_shuff_mean.npy")
#second_error_fid = np.load("data_many/second_rat_"+proxy+"_shuff_err.npy")
#third_fid = np.load("data_many/third_rat_"+proxy+"_shuff_mean.npy")
#third_error_fid = np.load("data_many/third_rat_"+proxy+"_shuff_err.npy")

nprops = len(opts)

for i in range(nprops):
    opt = opts[i]
    print(opt)

    second = np.load("data_many/second_"+proxy+"_"+opt+"_mean.npy")[chosen]
    second_err = np.load("data_many/second_"+proxy+"_"+opt+"_err.npy")[chosen]
    third = np.load("data_many/third_"+proxy+"_"+opt+"_mean.npy")[chosen]
    third_err = np.load("data_many/third_"+proxy+"_"+opt+"_err.npy")[chosen]


    print("Option: true",file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("<D_g^2> = ",second_fid,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("error D_g^2 = ",second_error_fid,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("<D_g^3> = ",third_fid,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("error D_g^3 = ",third_error_fid,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("Option:", opt,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("<D_g_opt^2> = ",second,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("error D_g_opt^2 = ",second_err,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("<D_g_opt^3> = ",third,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
    print("error D_g_opt^3 = ",third_err,file=open("txt/cums_"+proxy+"_"+opt+".txt","a"))
