import os
import shutil
import stat
import subprocess
import numpy as np

from sklearn.metrics import euclidean_distances
from sklearn.isotonic import IsotonicRegression

from . import write
from .config import parse
from shelljob import job


max_iter = 5

CMD_PM = ('%s -o %s -w 8 '
          '-r %d '
          '-k %s '
          '-i %s '
          '-q %s '
          '-d %f '
          '-c %s -y 1 -a %f -b %f > %s')


CMD_MDS = ('%s -o %s -w 8 '
           '-r %d '
           '-k %s '
           '-i %s '
           '-q %s '
           '-d %f '
           '-c %s -y 1 > %s')


def centerInteractions(Outpath):
    if os.path.exists(Outpath):
       spatialModel = np.loadtxt(Outpath, skiprows=1)
       print "model has dim of %s %s" % (spatialModel.shape[0], spatialModel.shape[1])
       center  = spatialModel[:,2:5].sum(axis=0)/spatialModel.shape[0]
       spatialModel[:, 2:5] = spatialModel[:, 2:5] - center

       rows = np.array(['chrom', 'locus', '3D_x', '3D_y','3D_z'], dtype='|S20')[:, np.newaxis]
       with open(Outpath, 'w') as f:
            np.savetxt(f, np.hstack((rows, spatialModel)), delimiter=', ', fmt='%s')


def run_mds(directory):
    print "Running MDS"

    if os.path.exists(os.path.join(directory,
                                   "config.ini")):
        config_file = os.path.join(directory, "config.ini")
    else:
        config_file = None

    options = parse(config_file)

    random_state = np.random.RandomState(seed=options["seed"])

    # First, compute MDS
    counts = np.load(os.path.join(directory,
                                  options["counts"]))
    wish_distances = compute_wish_distances(counts, alpha=options["alpha"],
                                            beta=options["beta"])
    wish_distances *= 1000
    wish_distances[np.isinf(wish_distances) | np.isnan(wish_distances)] = 0
    lengths = np.loadtxt(
        os.path.join(directory, options["organism_structure"]))

    try:
        len(lengths)
    except TypeError:
        lengths = np.array([lengths])

    write(wish_distances,
          os.path.join(directory,
                       'wish_distances.txt'),
          lengths=lengths,
          resolution=options["resolution"])

    # Write initial point
    X = 1. - 2. * random_state.rand(len(counts) * 3)
    np.savetxt(os.path.join(directory,
                            "MDS." + options["output_name"] + ".temp.txt"),
               X)

    locus_coord = options["output_name"].replace(".pdb",".bed")

    if not os.path.exists(os.path.join(directory,"model")) :
        os.mkdir(os.path.join(directory,"model"))

    locus_coord = os.path.join("model", locus_coord)

    cmd = CMD_MDS % (options["binary_mds"],
                     os.path.join(directory,
                                  "MDS." + options["output_name"]),
                     options["resolution"],
                     os.path.join(directory,
                                  options["organism_structure"]),
                     os.path.join(directory,
                                  "wish_distances.txt"),
                     os.path.join(directory,
                                  locus_coord),
                     options["adjacent_beads"],
                     options["chromosomes"],
                     os.path.join(directory,
                                  'MDS.log'))

    filename = os.path.join(directory, 'MDS.sh')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    fileptr = open(filename, 'wb')
    fileptr.write(cmd)
    fileptr.close()
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IXUSR)
    p = subprocess.Popen(filename.split(), shell='True')
    p.wait()
    ## TODO: center interactions
    #print "centering %s" % locus_coord
    #centerInteractions(locus_coord)
    print "Done!"
    return True


def mds(wish_distances, X=None, random_state=None):
    """
    Computes the MDS solution

    Parameters
    ----------
    wish_distances : ndarray
        Array of wish distances
    """
    if len(wish_distances.shape) != 2:
        raise ValueError("The wish distances should be a 2D ndarray.")
    if wish_distances.shape[0] == wish_distances.shape[1]:
        raise ValueError(
            "The wish distances should be a squared array (of "
            "shape (n by n).")
    n = wish_distances.shape[0]
    if X is not None and len(X) != n:
        raise ValueError(
            "The length of ")


def compute_wish_distances(counts, alpha=-3, beta=1):
    """
    Computes wish distances from a counts matrix

    Parameters
    ----------
    counts : ndarray
        Interaction counts matrix

    alpha : float, optional, default: -3
        Coefficient of the power law

    beta : float, optional, default: 1
        Scaling factor

    Returns
    -------
    wish_distances
    """
    if beta == 0:
        raise ValueError('beta cannot be equal to 0.')
    wish_distances = counts.copy()/ beta
    wish_distances[wish_distances != 0] **= 1. / alpha
    return wish_distances


def run_nmds(directory):
    print directory

    if os.path.exists(os.path.join(directory,
                                   "config.ini")):
        config_file = os.path.join(directory, "config.ini")
    else:
        config_file = None

    options = parse(config_file)
    run_mds(directory)

    for i in range(0, max_iter):
        if i == 0:
            try:
                X = np.loadtxt(
                    os.path.join(directory,
                                 "MDS." + options["output_name"] + ".txt"))
            except IOError:
                return
        else:
            X = np.loadtxt(
                os.path.join(directory,
                             '%d.NMDS.' % (i) + options["output_name"] +
                             ".txt"))

        X = X.reshape((len(X) / 3, 3))

        dis = euclidean_distances(X) * 1000
        counts = np.load(
            os.path.join(directory, options["counts"]))
        counts[np.isnan(counts)] = 0

        wish_distances = np.zeros(counts.shape)

        print "Fitting isotonic regression..."
        ir = IsotonicRegression()
        wish_distances[counts != 0] = ir.fit_transform(
            1. / counts[counts != 0],
            dis[counts != 0])
        print "writing wish distances"

        lengths = np.loadtxt(
            os.path.join(directory, options["organism_structure"]))

        try:
            len(lengths)
        except TypeError:
            lengths = np.array([lengths])

        write(wish_distances,
              os.path.join(directory,
                           '%d.NMDS.wish_distances.txt' % i),
              lengths=lengths, resolution=options["resolution"])

        if i == 0:
            shutil.copy(
                os.path.join(directory,
                             "MDS." + options["output_name"] + ".txt"),
                os.path.join(directory,
                             '%d.NMDS.' % (i + 1) + options["output_name"] +
                             ".temp.txt"))
        else:
            shutil.copy(
                os.path.join(directory,
                             '%d.NMDS.' % i + options["output_name"] + ".txt"),
                os.path.join(directory,
                             '%d.NMDS.' % (i + 1) + options["output_name"] +
                             ".temp.txt"))

        locus_coord = options["output_name"].replace(".pdb",".bed")

        cmd = CMD_MDS % (options["binary_mds"],
                         os.path.join(directory,
                                      "%d.NMDS." % (i + 1) +
                                      options["output_name"]),
                         options["resolution"],
                         os.path.join(directory,
                                      options["organism_structure"]),
                         os.path.join(directory,
                                      "%d.NMDS.wish_distances.txt" % (i)),
                         os.path.join(directory,
                                      locus_coord),
                         options["adjacent_beads"],
                         options["chromosomes"],
                         os.path.join(directory,
                                      str(i + 1) + '.NMDS.log'))

        filename = os.path.join(directory, str(i + 1) + '.NMDS.sh')
        fileptr = open(filename, 'wb')
        fileptr.write(cmd)
        fileptr.close()
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IXUSR)
        p =subprocess.Popen(filename.split(), shell='True')
        p.wait()
def run_pm1(directory):
    if os.path.exists(os.path.join(directory,
                                   "config.ini")):
        config_file = os.path.join(directory, "config.ini")
    else:
        config_file = None

    options = parse(config_file)
    run_mds(directory)
    shutil.copy(
        os.path.join(directory,
                     "MDS." + options["output_name"] + ".txt"),
        os.path.join(directory,
                     'PM1.' + options["output_name"] +
                     ".temp.txt"))

    counts = np.load(os.path.join(directory, options["counts"]))
    lengths = np.loadtxt(
        os.path.join(directory,
                     options["organism_structure"]))

    try:
        len(lengths)
    except TypeError:
        lengths = np.array([lengths])

    write(counts,
          os.path.join(directory,
                       'counts.txt'),
          lengths=lengths,
          resolution=options["resolution"])

    locus_coord = options["output_name"].replace(".pdb",".bed")
    cmd = CMD_PM % (options["binary_pm"],
                    os.path.join(directory,
                                 "PM1." + options["output_name"]),
                    options["resolution"],
                    os.path.join(directory, options["organism_structure"]),
                    os.path.join(directory,
                                 'counts.txt'),
                    os.path.join(directory,
                                 locus_coord),
                    options["adjacent_beads"],
                    options["chromosomes"],
                    options["alpha"],
                    options["beta"],
                    os.path.join(directory,
                                 'PM1.log'))

    filename = os.path.join(directory, 'PM1.sh')
    fileptr = open(filename, 'wb')
    fileptr.write(cmd)
    fileptr.close()
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IXUSR)
    p =subprocess.Popen(filename.split(), shell='True')
    p.wait()
    return True


def run_pm2(directory):
    from .poisson_model_power_law import estimate_alpha_beta
    if os.path.exists(os.path.join(directory,
                                   "config.ini")):
        config_file = os.path.join(directory, "config.ini")
    else:
        config_file = None

    options = parse(config_file)
    alpha, beta = options["alpha"], options["beta"]
    run_mds(directory)

    # Save counts in a format the C++ code can use
    counts = np.load(
        os.path.join(directory, options["counts"]))
    lengths = np.loadtxt(
        os.path.join(directory,
                     options["organism_structure"]))

    try:
        len(lengths)
    except TypeError:
        lengths = np.array([lengths])

    write(counts,
          os.path.join(directory,
                       'counts.txt'),
          lengths=lengths,
          resolution=options["resolution"])

    for i in range(0, max_iter):
        if i == 0:
            try:
                X = np.loadtxt(
                    os.path.join(directory,
                                 "MDS." + options["output_name"] + ".txt"))
            except IOError:
                return
        else:
            X = np.loadtxt(
                os.path.join(directory,
                             '%d.PM2.' % (i) + options["output_name"] +
                             ".txt"))

        X = X.reshape((len(X) / 3, 3))

        # Fit alpha and beta
        alpha, beta = estimate_alpha_beta(counts, X, ini=np.array([alpha]))

        if i == 0:
            shutil.copy(
                os.path.join(directory,
                             "MDS." + options["output_name"] + ".txt"),
                os.path.join(directory,
                             '%d.PM2.' % (i + 1) + options["output_name"] +
                             ".temp.txt"))
        else:
            shutil.copy(
                os.path.join(directory,
                             '%d.PM2.' % i + options["output_name"] + ".txt"),
                os.path.join(directory,
                             '%d.PM2.' % (i + 1) + options["output_name"] +
                             ".temp.txt"))

        locus_coord = options["output_name"].replace(".pdb",".bed")

        cmd = CMD_PM % (options["binary_pm"],
                        os.path.join(directory,
                                     "%d.PM2." % (i + 1) +
                                     options["output_name"]),
                        options["resolution"],
                        os.path.join(directory,
                                     options["organism_structure"]),
                        os.path.join(directory,
                                     'counts.txt'),
                        os.path.join(directory,
                                     locus_coord),
                        options["adjacent_beads"],
                        options["chromosomes"],
                        alpha,
                        beta,
                        os.path.join(directory,
                                     str(i + 1) + '.PM2.log'))

        filename = os.path.join(directory, str(i + 1) + '.PM2.sh')
        fileptr = open(filename, 'wb')
        fileptr.write(cmd)
        fileptr.close()
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IXUSR)
        p =subprocess.Popen(filename.split())
        p.wait()
