from copy import deepcopy
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
from skimage import morphology
import scipy
import random
import copy
import time


def isolateBody(ctScan):
    img = nib.load(ctScan)
    img_data = img.get_data()
    low_values_flags = img_data < -500
    high_values_flags = img_data > 2000
    img_data[:] = 1
    img_data[low_values_flags] = 0
    img_data[high_values_flags] = 0
    label, connectivity = measure.label(img_data, return_num=True)
    properties = measure.regionprops(label)

    biggest_CC = 0
    print(biggest_CC)

    for prop in properties:
        if (prop.area > biggest_CC):
            biggest_CC = prop.area
        print("bgst is ", prop.area)
    label = morphology.remove_small_objects(label, biggest_CC - 2)
    label, connectivity = measure.label(label, return_num=True)
    print("connecticity ", connectivity)

    name_format = img.get_filename().split('.')[0] + "bodysegmentation.nii.gz"

    low_values_flags = label == False
    high_values_flags = label == True

    img_data[low_values_flags] = 0
    img_data[high_values_flags] = 1

    print(name_format)
    nib.save(img, name_format)

    return img_data


def find_ROI(CT_scan, aorta_seg):
    print("looking for roi")

    aorta = nib.load(aorta_seg)
    aorta_data = aorta.get_data()
    img = nib.load(CT_scan)
    img_data = img.get_data()
    tuple_of_aorta_seg = np.nonzero(aorta_data)
    aorta_x = tuple_of_aorta_seg[0]
    aorta_y = tuple_of_aorta_seg[1]
    aorta_z = tuple_of_aorta_seg[2]
    minimal_z = np.amin(aorta_z)
    maximal_z = np.amax(aorta_z)
    minimal_x = np.amin(aorta_x)
    maximal_x = np.amax(aorta_x)
    minimal_y = np.amin(aorta_y)
    maximal_y = np.amax(aorta_y)
    half_z = int(np.floor((minimal_z + maximal_z) / 2))
    quarter_z = int(np.floor((minimal_z + maximal_z) / 4))
    eighth_z = int(np.floor((minimal_z + maximal_z) / 8))
    third_z = int(np.floor((minimal_z + maximal_z) / 3))
    tenth_z = int(np.floor((minimal_z + maximal_z) / 10))
    tent_to_2 = int(tenth_z / 2)
    seventh_z = int((minimal_z + maximal_z) / 7)
    sixth_z = int((minimal_z + maximal_z) / 6)
    fifth_z = int((minimal_z + maximal_z) / 5)
    print(third_z)
    third_z_2 = 2 * third_z

    img_data[:] = 0
    print(minimal_z, eighth_z, quarter_z)
    img_data[int(1.2 * maximal_x):int(1.6 * maximal_x), int(1.2 * minimal_y):maximal_y,
    minimal_z + fifth_z:half_z + tent_to_2] = 1
    name_format = img.get_filename().split('.')[0] + "ROI.nii.gz"

    nib.save(img, name_format)

    return img_data


def findSeeds(CT_scan, ROI):
    wanted_points = 200

    img = nib.load(CT_scan)
    img_data = img.get_data()
    seeds = copy.deepcopy(img_data)
    seeds[:] = 0
    img_roi = nib.load(ROI)
    roi_data = img_roi.get_data()
    tuple_of_roi_seg = np.nonzero(roi_data)
    roi_x = tuple_of_roi_seg[0]
    roi_y = tuple_of_roi_seg[1]
    roi_z = tuple_of_roi_seg[2]
    minimal_z = np.amin(roi_z)
    maximal_z = np.amax(roi_z)
    minimal_x = np.amin(roi_x)
    maximal_x = np.amax(roi_x)
    minimal_y = np.amin(roi_y)
    maximal_y = np.amax(roi_y)
    z_axis = list(range(minimal_z, maximal_z))
    print(z_axis)
    x_axis = list(range(minimal_x, maximal_x))
    y_axis = list(range(minimal_y, maximal_y))

    cur_points = 0
    while cur_points <= wanted_points:
        # print("cur_points " , cur_points)
        z_point = random.choice(z_axis)
        x_point = random.choice(x_axis)
        y_point = random.choice(y_axis)

        # print ("x y z " , x_point , y_point , z_point )
        gray_level = img_data[x_point, y_point, z_point]
        print(gray_level)
        if (gray_level > 0) and (gray_level < 150):
            seeds[x_point, y_point, z_point] = 1
            cur_points += 1
        else:
            print("point has not been choosefd")
            continue

    img_data[::] = seeds
    name_format = img.get_filename().split('.')[0] + "seeds.nii.gz"

    nib.save(img, name_format)
    return 0


def get_neighboors():
    pass


def find_similarity_func(voxel, ct_data, seeds):
    mone = ct_data[voxel] - np.mean(ct_data[seeds])
    return np.abs(mone / np.std(ct_data[seeds]))


def find_std(img_data, labeled_array):
    pass


def multipleSeedsRG(CT_scan, Aorta):
    img = nib.load(CT_scan)
    ROI = img.get_filename().split('.')[0] + "ROI.nii.gz"

    find_ROI(CT_scan, Aorta)

    img_data = img.get_data().astype(np.uint8)
    xdim, ydim, zdim = img_data.shape
    findSeeds(CT_scan, ROI)

    checked_voxels = np.zeros((xdim, ydim, zdim))
    checked_voxels.astype(np.uint8)
    new_neighboors_voxels = np.zeros((xdim, ydim, zdim)).astype(np.uint8)
    new_neighboors_voxels.astype(np.uint8)
    seeds_name_format = img.get_filename().split('.')[0] + "seeds.nii.gz"

    seeds = nib.load(seeds_name_format)
    seeds_data = seeds.get_data()
    checked_voxels[seeds_data == 1] = 1
    print("checked voxels after insert seeds", np.sum(checked_voxels), " while seeds is", np.sum(seeds_data))
    next_segment = morphology.dilation(seeds_data, morphology.cube(3, np.uint8))
    new_neighboors_voxels = np.subtract(next_segment, seeds_data)
    new_neighboors_voxels[checked_voxels == 1] = 0
    checked_voxels[new_neighboors_voxels == 1] = 1
    print("checked voxels after insert first neighboors", np.sum(checked_voxels), " while seeds is", np.sum(seeds_data))
    iteration = 0
    iterations = []
    new_voxels_added_to_list = []
    new_voxels_added_to_segmentation = []
    while (np.sum(new_neighboors_voxels) > 0):
        iteration += 1
        print("sum of seeds ", np.sum(seeds_data == 1), "iteration", iteration)
        mean = np.mean(img_data[seeds_data == 1])

        print("mean is", mean)
        std = np.std(img_data[seeds_data == 1])
        neighboors_grays = np.zeros((xdim, ydim, zdim))

        # print ("ng shape " , new_neighboors_voxels.shape)

        neighboors_grays[new_neighboors_voxels == 1] = img_data[new_neighboors_voxels == 1]
        seeds_before = np.sum(seeds_data)
        seeds_data[np.absolute(neighboors_grays - mean) < 20] = 1
        seeds_after = np.sum(seeds_data)
        seeds_just_added = seeds_after - seeds_before
        print("seeds added", seeds_just_added)
        # print ("seeds sum aftet insertion via mean", np.sum(seeds_data))

        next_segment = morphology.dilation(seeds_data, morphology.cube(3, np.uint8))
        # print("seeds sum aftet dialition", np.sum(seeds_data))
        new_neighboors_voxels = np.subtract(next_segment, seeds_data)
        # print("new neibors after sumstract dialition from seeds", np.sum(new_neighboors_voxels))
        new_neighboors_voxels[checked_voxels == 1] = 0
        # print ("new neibors after substract who checked " , np.sum(new_neighboors_voxels))

        checked_voxels[new_neighboors_voxels == 1] = 1
        iterations.append(iteration)
        new_voxels_added_to_list.append(np.sum(new_neighboors_voxels))
        new_voxels_added_to_segmentation.append(seeds_just_added)

        if (iteration == 100):
            break
        # if(iteration%30==0):
        #
        #     name_format = img.get_filename().split('.')[0] + "after adding" + str(np.sum(seeds_data==1))+"in iteration" + str(iteration) +".nii.gz"
        #     nib.save(seeds, name_format)

    name_format = img.get_filename().split('.')[0] + "AfterRegionGrwing.nii.gz"
    nib.save(seeds, name_format)
    # lines = plt.plot(iterations,new_voxels_added_to_list , iterations, new_voxels_added_to_segmentation,'o')
    # plt.setp(lines[0],linewidth = 4)
    # plt.setp(lines[1],linewidth = 2)
    # plt.legend(("new_voxels_added_to_list","new_voxels_added_to_segmentation"),loc='upper right')
    # plt.title("data")
    # plt.show()
    return 0


def oneConnectorComponent(segmentation):
    """
    opens a segmentation and make 1 connector component, as well as other morphological operations :
    erosion and dialition
    :param segmentation: already been saved
    :param aorta:
    :param after:
    :return:
    """

    seg = nib.load(segmentation)
    seg_data = seg.get_data()

    label, connectivity = measure.label(seg_data, return_num=True)
    print ("first connectivity" , connectivity)
    label = morphology.binary_erosion(label)
    label, connectivity = measure.label(label, return_num=True)
    print ("connectivity after 1 erusion connectivity" , connectivity)
    # label = morphology.binary_erosion(label)
    # label, connectivity = measure.label(label, return_num=True)
    # print ("connectivity after 2 erusions connectivity" , connectivity)
    properties = measure.regionprops(label)

    biggest_CC = 0
    print(biggest_CC)

    for prop in properties:
        if (prop.area > biggest_CC):
            biggest_CC = prop.area
        print("bgst is ", biggest_CC)
    label = morphology.remove_small_objects(label, biggest_CC -2) #todo delete!
    label, connectivity = measure.label(label, return_num=True)
    print("connecticity ", connectivity)
    for i in range (4):

        label = morphology.binary_dilation(label)
    # label = morphology.binary_dilation(label)

    name_format = seg.get_filename().split('_')[0] + "after_1_component.nii.gz"
    # low_values_flags = label == False
    # high_values_flags = label == True

    seg_data[label == False] = 0
    seg_data[label == True] = 2

    print(name_format)
    nib.save(seg, name_format)

    return seg_data


def segment_liver(ctFileName, AortaFileName, outputFileName):
    img = nib.load(ctFileName)
    img_data = img.get_data()
    aorta = nib.load(AortaFileName)
    multipleSeedsRG(ctFileName, AortaFileName)
    name_format = name_format = img.get_filename().split('.')[0] + "AfterRegionGrwing.nii.gz"
    print ( "name after region groeing" , name_format)
    time.sleep(5)
    oneConnectorComponent(name_format)

    name_format = img.get_filename().split('_')[0] + "after_1_component.nii.gz"
    print("name to load" ,name_format)
    img = nib.load(name_format)
    nib.save(img, outputFileName)

    return 0


def evaluateSegmentation(groundTruthSeg, estimatedSeg):
    groundTruth = nib.load(groundTruthSeg)
    groundTruthData = groundTruth.get_data()
    seg = nib.load(estimatedSeg)
    seg_data = seg.get_data()
    union = np.sum(seg_data) + np.sum(groundTruthData)
    intersection = np.logical_and(seg_data, groundTruthData)
    intersection = np.sum(intersection)

    divis = np.float(np.float(intersection) / np.float(union))

    #
    DC = 2 * (divis)
    VOD = (1 - divis)
    print("d and v ", DC, VOD)

    return VOD, DC




if __name__ == "__main__":
    # isolateBody('Case1_CT.nii.gz')
    # oneConnectorComponent('Case1_CTAftersegmentation50.nii.gz')
    # findSeeds('Case1_CT.nii.gz','Case1_CTROI.nii.gz')
    # print("hello world")
    # find_ROI('Case1_CT.nii.gz','Case1_Aorta.nii.gz' )
    # segment_liver('Case1_CT.nii.gz', 'Case1_Aorta.nii.gz', "first.nii.gz")
    # segment_liver('Case2_CT.nii.gz', 'Case2_Aorta.nii.gz', "sec.nii.gz")
    # segment_liver('Case3_CT.nii.gz', 'Case3_Aorta.nii.gz', "third.nii.gz")
    # segment_liver('Case4_CT.nii.gz', 'Case4_Aorta.nii.gz', "fourth.nii.gz")
    evaluateSegmentation('Case1_liver_segmentation.nii.gz', 'first.nii.gz')

