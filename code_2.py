import gc, csv, re, pathlib, sys
import numpy as np
from osgeo import gdal

AOH_QUANTI = True
WARP_FACTOR = 5
BOOTS = "pred_boots_bobyqa.csv"             # Outputs of logistic traits model
PREDS = "pred_full_data.csv"                # Outputs of logistic traits model
SPCONVERT = "conversion_pious_avonet.csv"   # Species name conversion
AOH = "./wn_piou_aoh_warped"                # Path to AOH 500m resolution files
OUTPUT_DIR = "./wn_piou_boot_aoh"           # Path to bootstrapped BRI files
STEP = 10                                   # Size of batches
NBOOTS = 1000                               # Number of bootstrap replicates

msktif = gdal.Open("masque_warped.tif")
mskbnd = msktif.GetRasterBand(1)
mskdata = mskbnd.ReadAsArray()

predrisk = {}
with open(PREDS) as csvfile:
    reader = csv.reader(csvfile, delimiter=";")
    header = next(reader)
    for pred in reader:
        predrisk[pred[0]] = [float(pred[1])]
with open(BOOTS) as csvfile:
    reader = csv.reader(csvfile, delimiter=";")
    header = next(reader)
    for pred in reader:
        predrisk[pred[0]] = predrisk[pred[0]]+[float(x) for x in pred[1:]]

aohdir = [str(x) for x in pathlib.Path(AOH).iterdir() if x.is_file]
spconvert = {}
with open(SPCONVERT) as csvfile:
    for sp in csv.DictReader(csvfile):
        spconvert[sp["pious"]] = sp["avonet"]

index = -1
while index < NBOOTS:
    if index == -1:
        ixrange = [-1]
    elif index == 0:
        ixrange = [0]
    else:
        ixrange = [i for i in range(index, index+STEP) if i <= NBOOTS]
    piouvec = [np.zeros(mskdata.shape, dtype=np.float32) for ix in ixrange]
    print("Computing", [ixrange[0], ixrange[-1]])
    sys.stdout.flush()
    for piou, pred in sorted(predrisk.items()):
        species = spconvert[piou]
        prefix = species.replace(" ", "_")
        nbdata = 0
        spdata = None
        for filename in sorted(aohdir):
            if re.search("^%s/%s.tif$" % (AOH, prefix), filename) or \
               re.search("^%s/%s_B.tif$" % (AOH, prefix), filename) or \
               re.search("^%s/%s_R.tif$" % (AOH, prefix), filename):
                tif = gdal.Open(filename)
                bnd = tif.GetRasterBand(1)
                data = bnd.ReadAsArray()
                data = np.where(data==255, 0, data)
                nbdata = nbdata + 1
                if spdata is None:
                    spdata = data
                else:
                    spdata = spdata + data
        if AOH_QUANTI:
            spdata = spdata / (nbdata * float(WARP_FACTOR**2))
        else:
            spdata = np.where(spdata > 0, 1, 0)
        for i in range(len(ixrange)):
            ix = ixrange[i]
            if ix == -1:
                pr = np.float32(pred[ix])
            else:
                pr = np.float32(1.)
            piouvec[i] = piouvec[i] + pr*spdata
    for i in range(len(ixrange)):
        ix = ixrange[i]
        piouarr = piouvec[i] / np.max(piouvec[i])
        piouarr = np.where(mskdata==200, np.float32(-1), piouarr)
        driver = gdal.GetDriverByName('GTiff')
        if ix == -1:
            tifname =  OUTPUT_DIR + "/richness.tif"
        elif ix == 0:
            tifname =  OUTPUT_DIR + "/predicted.tif"
        else:
            tifname = OUTPUT_DIR + "/b%d.tif" % ix
        print("saving ", tifname)
        sys.stdout.flush()
        raster = driver.Create(tifname, piouarr.shape[1], piouarr.shape[0], 
                               1, gdal.GDT_Float32)
        raster.SetProjection(msktif.GetProjection())
        raster.SetGeoTransform(msktif.GetGeoTransform())
        band = raster.GetRasterBand(1)
        band.WriteArray(piouarr)
        band.SetNoDataValue(-1)
        band.FlushCache()
        raster = None
    index = ixrange[-1]+1

