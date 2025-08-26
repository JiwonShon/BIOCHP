import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
from skimage.io import *
from skimage.color import *
from skimage.morphology import *
from skimage.filters.rank import *
from scipy.ndimage import label, generate_binary_structure, binary_fill_holes

NUM_CLASSES = 4
EC_SIZE_THRESHOLD = 20

def mask_postprocess(img):
    #if ecDNA is touching chromosome/nuclei, mark that whole
    #component as that class
    def merge_comp(img, class_id=2):  
        I = img
        mask_id = 1
        if(class_id == 1): mask_id = 2
        temp = I == mask_id
        I[temp] = 0 #remove everything with mask_id
        O = I
        s = generate_binary_structure(2,2)
        labeled_array, num_features = label(I,  structure=s)
        for i in range(1, num_features):
            ind = (labeled_array == i)
            if(np.any(I[ind]==class_id)):
                O[ind] = class_id
        img[opening(O, diamond(1)) == class_id] = class_id #reset nuclei and chromosomes in main image
        img[temp] = mask_id 
        return img

    #fill holes in connected components
    def fill_holes(img, class_id):
        temp = binary_fill_holes(img == class_id)
        img[temp == 1] = class_id
        return img

    def size_thresh(img):
        nuc_regs = measure.regionprops(measure.label(img == 1)) #Set small nucleic regions as background
        chrom_regs = measure.regionprops(measure.label(img == 2))

        c_area = [c.area for c in chrom_regs]
        if len(c_area):
            avg_chrom_size = np.mean(c_area)
        else:
            avg_chrom_size = 0
            
        for r in nuc_regs:
            if(r.area < avg_chrom_size):
                img[tuple(r.coords.T)] = 0 

        chrom_regs = measure.regionprops(measure.label(img == 2)) #Set small chromosomal regions as ecDNA
        ec_regs = measure.regionprops(measure.label(img == 3))

        c_area = [c.area for c in ec_regs]
        if len(c_area):
            avg_ec_size = np.mean(c_area)
        else:
            avg_ec_size = 0

        for r in chrom_regs:
            if(r.area < avg_ec_size):
                img[tuple(r.coords.T)] = 3

        for r in ec_regs: # set small ecs to background
            if(r.area < EC_SIZE_THRESHOLD):
                img[tuple(r.coords.T)] = 0
        return img

    img = fill_holes(fill_holes(img, 1), 2) #fill holes
    img = size_thresh(img)
    #dilation XOR erosion of ecDNA. Smoothens ecDNA borders by enlarging and then evenly eroding
    img[binary_dilation(img == 3, diamond(1)) ^ binary_erosion(img == 3, diamond(1))] = 0

    chrom_regs = measure.regionprops(measure.label(img == 2))
    nuc_regs = measure.regionprops(measure.label(img == 1))
    c_y = [c.centroid[0] for c in chrom_regs]; c_x= [c.centroid[1] for c in chrom_regs]
    n_cent = [n.centroid for n in nuc_regs]

    ######check if there is a nuclei in the middle of the metaphase###############
    min_chrom_count = 5; v=70
    for idx, n in enumerate(n_cent):
        count = 0
        ######## check if chroms are close to metaphase and if chroms are surrounding the nuclei###
        left = (len(np.where((c_x > n[1]) & (c_x < n[1]+v))[0])> min_chrom_count)
        right = (len(np.where((c_x < n[1]) & (c_x > n[1]-v))[0]) > min_chrom_count)
        bottom = (len(np.where((c_y < n[0]) & (c_y > n[0]-v))[0])> min_chrom_count)
        top = (len(np.where((c_y > n[0]) & (c_y < n[0]+v))[0]) > min_chrom_count)
        if((left*bottom & right*top) or (bottom*right & top*left)): #check if two opposing quadrants contain chrs
            img[tuple(nuc_regs[idx].coords.T)] = 0
    img = merge_comp(merge_comp(img, 1), 2)
    img[binary_dilation(img == 3, diamond(1))] = 3
    return img


def count_features(image: np.ndarray, threshold: int):
    bin = np.where(image > 32, 1, 0)
    bin = measure.label(bin) # Connected component

    intensity = []
    blob_size = []
    component_num = 0
    mean_intensity = 0
    mean_blob_size = 0

    for idx in np.unique(bin)[1:]:
        component = np.where(bin==idx, image, 0)
        pixel_num = np.count_nonzero(component)
        if pixel_num > threshold:
            component_num += 1
            intensity.append(np.sum(component) / pixel_num)
            blob_size.append(pixel_num)

    if component_num:
        mean_intensity = np.mean(np.array(intensity))
        mean_blob_size = np.mean(np.array(blob_size))

    return component_num, mean_intensity, mean_blob_size

def chrom_count(mask):
    return len(np.unique(measure.label(mask))) - 1

def chrom_inten(chrom_image):
    return np.sum(chrom_image)

def nuclei_count(mask):
    return len(np.unique(measure.label(mask))) - 1

def nuclei_inten(nuclei_image):
    return np.sum(nuclei_image)

def interphase_count(image, mask):
    mask = np.where(mask==1, 1, 0)
    mask = binary_fill_holes(mask).astype(int)
    mask = mask[...,np.newaxis]
    image = image * mask
    r = image[...,0]
    g = image[...,1]

    r_comp, r_inten, r_size = count_features(r, EC_SIZE_THRESHOLD)
    g_comp, g_inten, g_size = count_features(g, EC_SIZE_THRESHOLD)

    nuclei_num = nuclei_count(mask)
    if nuclei_num:
        mean_nuclei_inten = nuclei_inten(image) / nuclei_num
    else: mean_nuclei_inten = 0

    return r_comp, g_comp, r_inten, g_inten, r_size, g_size, nuclei_num, mean_nuclei_inten


def metaphase_count(image, mask):
    chrom_mask = np.where(mask==2, 1, 0)
    chrom_mask = binary_fill_holes(chrom_mask).astype(int)
    chrom_mask = chrom_mask[..., np.newaxis]
    chrom_image = image * chrom_mask

    nuclei_mask = np.where(mask==1, 1, 0)
    nuclei_mask = nuclei_mask[..., np.newaxis]

    ec_mask = np.logical_not(np.logical_or(chrom_mask, nuclei_mask))
    ec_image = image * ec_mask

    r = ec_image[...,0]
    g = chrom_image[...,1]

    r_comp, r_inten, r_size = count_features(r, EC_SIZE_THRESHOLD)
    g_comp, g_inten, g_size = count_features(g, EC_SIZE_THRESHOLD)

    chrom_num = chrom_count(chrom_mask)

    if chrom_num:
        mean_chrom_inten = chrom_count(chrom_image) / chrom_num
    else: mean_chrom_inten = 0

    return r_comp, g_comp, r_inten, g_inten, r_size, g_size, chrom_num, mean_chrom_inten



def feature_extraction(image, mask):
    # mask = mask_postprocess(mask)
    return np.array([interphase_count(image, mask), metaphase_count(image, mask)]).reshape(-1)
