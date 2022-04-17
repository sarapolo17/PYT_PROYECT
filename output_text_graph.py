###################################################
############## TEXT AND PLOT OUTPUT ###############
###################################################
# Libraries
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from statistics import mean
import numpy as np

def create_output_text (prefix, complete_flex_list):

    '''Creates parseable text file with the flexibility scores'''

    fd = open(prefix + ".txt", "w")
    fd.write("Position\t\t\tResidue\t\t\tFlexibility\t\tConfidence\n")

    for flex_tuple in complete_flex_list:
        fd.write(str(flex_tuple[0]) + "\t\t\t\t" + flex_tuple[1] + "\t\t\t" + "%+.3f\t\t\t" % (flex_tuple[2]) + str(flex_tuple[3]) + "\n")

    fd.close()


def draw_flex_line(prefix, flexibility_results):

    '''Creates plot of the flexibility results'''

    pos = []
    res = []
    flex = []
    conf = []

    for flex_tuple in flexibility_results:

        pos.append(flex_tuple[0])
        res.append(flex_tuple[1])
        flex.append(flex_tuple[2])
        conf.append(flex_tuple[3])

    plt.plot(pos,flex)
    plt.axhline(y=0, color='r', linestyle='-')
    plt.xlabel("Residue position")
    plt.ylabel("Flexibility score")

    plt.savefig(prefix + '.png')
    plt.show()

def draw_flex_line_col(prefix, flexibility_results):

    '''Creates plot of the flexibility results colored by confidence
    score'''

    pos = []
    res = []
    flex = []
    conf = []

    for flex_tuple in flexibility_results:

        pos.append(flex_tuple[0])
        res.append(flex_tuple[1])
        flex.append(flex_tuple[2])
        conf.append(flex_tuple[3])

    # Arrays

    pos = np.asarray(pos)
    flex = np.asarray(flex)
    conf = np.asarray(conf)


    # Color map

    cmap = ListedColormap(['r', 'y', 'g'])
    norm = BoundaryNorm([0, max(conf)/2, 0.75 * max(conf), max(conf)], cmap.N)

    # Line segments

    points = np.array([pos, flex]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(conf)

    # Plot

    plt.gca().add_collection(lc)
    plt.ylim(min(flex), max(flex))

    plt.xlim(0,len(pos))
    plt.axhline(y=0, color='k', linestyle='-')
    plt.xlabel("Residue position")
    plt.ylabel("Flexibility score")

    plt.savefig(prefix + '.png')
    plt.show()
