
#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import numpy as np

def plot_benchmarks(max_pointings, gpu=False):
    """
    #Galaxy MPB benchmarks serial, cal once upgrade
    galaxy_mpb_times = np.array([24.946077569999996, 16.736973759999998, 13.74829847,
                        12.11295253, 11.083652268, 10.641167811666666,
                        10.711637881428569, 9.822233005000001, 9.584228442222221,
                        9.22881431, 9.166687645454546, 9.4694958725,
                        9.33200447923077, 8.810422717142856, 8.674555905999998])
    galaxy_mpb_t_std = np.array([1.4452587252036433, 0.8276700183039727, 0.7197822452137351,
                        0.38995991243788464, 0.1323128179640961, 0.38698380790094133,
                        0.23876854414689958, 0.2991264454866008, 0.28994737767526846,
                        0.2038852330143943, 0.17738044480151316, 0.20314511164648275,
                        0.19631739333354165, 0.47328326372738955, 0.19456733116853275])
    galaxy_orig_times = np.array([1.23*24]*15)
    galaxy_orig_t_std = np.array([0.03*24]*15)

    #Ozstar MPB benchmarks serial, cal once upgrade
    ozstar_mpb_times = np.array([10.8919, 6.6485, 5.138566666666667,
                        4.611625, 4.37054, 4.0846833333333326,
                        3.665057142857143, 3.8473875, 3.688588888888889,
                        3.65932, 3.472590909090909, 3.642475,
                        3.6434384615384605, 3.1627214285714285, 3.0993533333333327])
    ozstar_mpb_t_std = np.array([0.5828762390079045, 0.38231331391935586, 0.2598150089754033,
                        0.3046873551281707, 0.3518930786474778, 0.33209539752239203,
                        0.32578035796016463, 0.3030123241614934, 0.300121486307145,
                        0.3495969816803344, 0.2543067229525066, 0.2715681916357903,
                        0.251464319230689, 0.24366376957200842, 0.22391966882989287])
    ozstar_orig_times = np.array([0.73*24]*15)
    ozstar_orig_t_std = np.array([0.02*24]*15)
 
    #Shangia ARM MPB benchmarks serial, cal once upgrade
    #sugon gpu
    arm_mpb_times = np.array([10.428699999999997, 6.327849999999999, 4.876733333333333, 4.2145, 3.8385600000000006, 3.6098166666666662, 3.3983857142857152, 3.2382874999999993, 3.1671000000000005, 3.1546199999999995, 3.1414727272727276, 3.0049666666666663, 2.968484615384616, 2.9391785714285716, 2.871186666666667])
    arm_mpb_t_std = np.array([1.144133239618533, 0.4706367999848716, 0.14700383970797806, 0.11682392734367397, 0.1553377944996001, 0.13968908233008837, 0.1601386790585342, 0.13881945502612375, 0.09439224112809899, 0.13325357931402818, 0.12566505526876778, 0.11650441574845508, 0.09220699860485634, 0.15198450495000457, 0.10387002048500599])
    arm_orig_times = np.array([43.710000/100.*24.]*15)
    arm_orig_t_std = np.array([0.02]*15)
    """

    pns = list(range(1,max_pointings+1))
    
    #Benchmarks for 60s of beamforming

    # Galaxy ------------------------------------------------------
    galaxy_orig_times = np.array([80.024356]*20)*pns
    galaxy_orig_t_std = np.array([17.407244805764243]*20)
    galaxy_ipfb_times = np.array([176.90150595833333]*20)*pns
    galaxy_ipfb_t_std = np.array([4.81427011600363]*20)
    galaxy_mpb_times = np.array([87.92744508333332, 130.752930625, 142.68917829166665, 151.09457491666666, 175.9199586666667,
                        190.59483129166665, 211.239435, 220.92454999999998, 237.09312275, 241.56801558333336,
                        263.16444916666666, 278.83963575, 295.05736541666664, 324.809873875, 338.5593945,
                        361.0044414583333, 381.0780106666666, 394.8712732916667, 415.33012729166666, 422.7796177916666])
    galaxy_mpb_t_std = np.array([15.547161555917718, 17.83388703707663, 14.421354619971813, 14.253923430757352, 14.94715700944577,
                        16.52515973362016, 11.542991266680383, 8.241497860268732, 7.9050627160796205, 8.75072582114132,
                        9.984636448279785, 8.207385441325538, 6.7822706699207576, 13.355794591824555, 9.976878315324516,
                        8.90584301787109, 16.0897684666793, 7.04310014752728, 13.628790278733309, 7.023970736723734])
    galaxy_mpb_gpu_times = np.array([0.46495833333333336, 0.7200000000000001, 0.976625, 1.232, 1.49,
                            1.7459999999999998, 2.002, 2.256, 2.511, 2.770000000000001,
                            3.0250000000000004, 3.2760000000000002, 3.5425000000000004, 3.7940000000000005, 4.050000000000001,
                            4.304, 4.556, 4.824, 5.073, 5.328333333333333])
    galaxy_mpb_gpu_t_std = np.array([0.01119584390223042, 0.1784189825476351, 0.25025506010602105, 0.31014109692202996, 0.3648629879831606,
                            0.4159326868617084, 0.46489246068311324, 0.511547217327546, 0.5568886782831917, 0.6008529816482194,
                            0.6438652809400426, 0.6856262344261145, 0.727397257922611, 0.7669700486691483, 0.806613290245084,
                            0.8453303101943838, 0.8839070463949624, 0.9214119599831552, 0.9587272741041174, 0.995414713383192])

    # Garrawarla --------------------------------------------------
    garrawarla_orig_times = np.array([36.61813077083333]*20)*pns
    garrawarla_orig_t_std = np.array([7.0984338765348305]*20)
    garrawarla_ipfb_times = np.array([45.22205279166667]*20)*pns
    garrawarla_ipfb_t_std = np.array([4.3707582012192505]*20)
    garrawarla_mpb_times = np.array([36.603427875, 31.402163833333333, 38.723946041666665, 42.19923445833333, 44.64464158333334,
                            48.18021158333334, 47.220577291666665, 55.08318729166667, 54.272442250000005, 63.23153141666668,
                            65.77896220833334, 70.14031154166668, 72.28991529166667, 76.61966149999999, 75.985575875,
                            80.09525933333333, 85.08665087499999, 90.50430125000001, 92.93317920833333, 95.44081741666666])
    garrawarla_mpb_t_std = np.array([14.936836223789541, 9.197380107737061, 8.004670682815602, 6.719730478378461, 8.626076415601606,
                            9.571553176642036, 9.41658963492266, 7.565667642952214, 9.876546570038121, 12.460640966743567,
                            10.892631988997312, 10.085567742822462, 8.259131706006734, 8.231193363069687, 8.785111966270026,
                            9.823887533315723, 9.279506162073455, 8.04683163106655, 9.213397394332917, 9.857900903074347])
    garrawarla_mpb_gpu_times =  np.array([0.12008333333333332, 0.15225, 0.18600000000000003, 0.22033333333333335, 0.25333333333333335,
                            0.28775, 0.32025, 0.35300000000000004, 0.38587499999999997, 0.42125000000000007,
                            0.4519166666666667, 0.484, 0.5189166666666667, 0.5489166666666666, 0.585,
                            0.6173333333333333, 0.6467083333333333, 0.6832499999999998, 0.7156666666666666, 0.7433333333333333])
    garrawarla_mpb_gpu_t_std = np.array([0.01595655001147989, 0.08219217570012903, 0.10941368363539668, 0.13163703214498143, 0.15082157540845217,
                           0.169256949214595, 0.18621798330757708, 0.20307807471354336, 0.21838650938881146, 0.23479197053048526,
                           0.2493708759036287, 0.2641602813839541, 0.278770219153938, 0.2930138104385162, 0.3070779624134562,
                           0.32091418923343, 0.3345482989319772, 0.34741005249171475, 0.3601623782279182, 0.3720426292590419])

    # Calculate improvement uncertainties and plot them
    """
    ozstar_improvement = ozstar_mpb_times/ozstar_orig_times
    ozstar_per_unc = ozstar_mpb_t_std/ozstar_mpb_times + ozstar_orig_t_std/ozstar_orig_times
    plt.errorbar(pns, ozstar_improvement, yerr=ozstar_improvement*ozstar_per_unc,
                 color='green', label='OzSTAR super computer')
    """
    galaxy_improvement = galaxy_orig_times/galaxy_mpb_times
    galaxy_per_unc = galaxy_mpb_t_std/galaxy_mpb_times + galaxy_orig_t_std/galaxy_orig_times
    plt.errorbar(pns, galaxy_improvement, yerr=galaxy_improvement*galaxy_per_unc,
                 color='blue', label='Galaxy super computer')
    garrawarla_improvement = garrawarla_orig_times/garrawarla_mpb_times
    garrawarla_per_unc = garrawarla_mpb_t_std/garrawarla_mpb_times + garrawarla_orig_t_std/garrawarla_orig_times
    plt.errorbar(pns, garrawarla_improvement, yerr=garrawarla_improvement*garrawarla_per_unc,
                 color='purple', label='Garrawarla  super computer')
    """
    plt.errorbar(pns, arm_orig_times/arm_mpb_times, yerr=arm_mpb_t_std/5,
                 color='red', label='CSRC prototype')
    
    # Galaxy garra comparison
    gg_improvement = galaxy_mpb_times/garrawarla_mpb_times
    gg_per_unc = galaxy_mpb_t_std/galaxy_mpb_times + garrawarla_mpb_t_std/garrawarla_mpb_times
    plt.errorbar(pns, gg_improvement, yerr=gg_improvement*gg_per_unc,
                 color='blue', label='Total processing improvement')
    gg_gpu_improvement = galaxy_mpb_gpu_times/garrawarla_mpb_gpu_times
    gg_gpu_per_unc = galaxy_mpb_gpu_t_std/galaxy_mpb_gpu_times + garrawarla_mpb_gpu_t_std/garrawarla_mpb_gpu_times
    plt.errorbar(pns, gg_gpu_improvement, yerr=garrawarla_mpb_gpu_t_std,
                 color='green', label='GPU processing improvement')
    """

    plt.ylabel("Factor of improved processing efficiency")
    plt.xlabel("Number of simultaneous tied-array beams")
    #plt.legend(loc='upper right', bbox_to_anchor=(0.95, 0.85))
    plt.legend(loc='upper left', bbox_to_anchor=(0.005, 0.995))
    plt.savefig("Beamformer_benchmark.eps")
    plt.savefig("Beamformer_benchmark.png", bbox_inches='tight', dpi=1000)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Plot beamformer benchmarks""")
    parser.add_argument('--max_pointing_num', type=int, default=20,
            help="Max number of pointings for multipixel beamformer")
    args = parser.parse_args()

    plot_benchmarks(args.max_pointing_num)