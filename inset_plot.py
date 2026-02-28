import matplotlib.pyplot as plt
from matplotlib.image import imread
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

def create_inset_plot(main_plot_file,inset_plot_file,inset_position,filename):
    # Load images
    img_main = imread(main_plot_file)
    img_inset = imread(inset_plot_file)

    fig = plt.figure(figsize=(15,15))

    ax_main = fig.add_axes([0, 0, 1, 1])
    ax_main.imshow(img_main)
    ax_main.axis("off")

    # Inset image (x, y, width, height in figure coords)
    ax_inset = fig.add_axes(inset_position)
    ax_inset.imshow(img_inset)
    ax_inset.axis("off")

    plt.savefig(filename, transparent=True)


# main_plot_file ="graphs/stress_tensor_manuscript/shear_stress_vs_radial_bins_p0.80.png"
# inset_plot_file ="graphs/stress_tensor_manuscript/shear_stress_vs_radial_bins_p0.00.png"
# filename ="graphs/stress_tensor_manuscript/shear_stress_vs_radial_bins_p0.80_inset.png"
# inset_pos =[0.16, 0.454, 0.45, 0.45]
# create_inset_plot(main_plot_file,inset_plot_file,inset_pos,filename)
if __name__ == "__main__":
    for l in [4,5,6]:
        main_plot_file = "Panel_2/SD_s0_scatter_l_{}.png".format(l)
        inset_plot_file = "Panel_2/s0_histogram_periodic_l_{}.png".format(l)
        filename ="Panel_2/SD_s0_with_inset_l_{}.png".format(l)
        width = 0.5
        inset_pos =[0.4, 0.15,width,width]
        create_inset_plot(main_plot_file,inset_plot_file,inset_pos,filename)
    # for l in [5,6]:
    #     main_plot_file = "Panel_5/SD_s0_scatter_l_{}.png".format(l)
    #     inset_plot_file = "Panel_5/s0_histogram_spheroid_l_{}.png".format(l)
    #     filename ="Panel_5/SD_s0_with_inset_l_{}.png".format(l)
    #     width = 0.55
    #     inset_pos =[0.38, 0.36,width,width]
    #     create_inset_plot(main_plot_file,inset_plot_file,inset_pos,filename)