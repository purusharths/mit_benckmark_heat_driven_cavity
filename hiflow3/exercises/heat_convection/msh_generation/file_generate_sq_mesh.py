import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta

mat_num = {"bottom_b": 11, "top_b": 12,
           "right_b":13, "left_b":14 }

def denseboundspace(size, start, end, alpha):
    x = np.linspace(0, 1, size)
    return start + beta.cdf(x, 2.-alpha, 2.-alpha) * (end-start)

def create_point_kernel(size=10, start=0, end=1, alpha = 0.5):
    return denseboundspace(size, start, end, alpha)

def cartesian_coord(*arrays):
    grid = np.meshgrid(*arrays)
    coord_list = [entry.ravel() for entry in grid]
    points = np.vstack(coord_list).T
    return points

def print_point_coords(fobj, cc, xsize, ysize):
    # cartesian_coord(n_x, n_y)
    offset = 0
    count = 0
    for y in range(0,ysize):
        # print("\n*** {}\n  ".format(y))
        for i in range(0,xsize):
            tmp = cc[i+offset]
            fobj.write("{} {} {} 0\n".format(count, tmp[0], tmp[1]))
            count+=1
        offset+=xsize    

def print_quad_box(fobj, x_max, y_max):
    count = 0
    offset = 0
    n = x_max
    for i in range(0,y_max-1):
        for j in range(0,n-1):
            fobj.write("{} -1 quad {} {} {} {}\n".format(count, j+offset, (j+offset)+1, ((j+offset)+1)+n, (j+offset)+n))
            count+=1
        offset+=n

def print_boundary_lines(fobj, total_sq_count, xsize, ysize):
    counter = total_sq_count
    # bottom boundary
    for i in range(xsize-1):
        fobj.write("{} {} line {} {}\n".format(counter, mat_num["bottom_b"], i, i+1))
        counter+=1
    # top boundary
    tb_offset = (ysize-1)*xsize
    for i in range(xsize-1):
        fobj.write("{} {} line {} {}\n".format(counter, mat_num["top_b"], i+tb_offset, (i+1)+tb_offset))
        counter+=1
    # left boundary
    left_b = list(range(0, tb_offset+1, xsize))
    for i in range(0, ysize-1):
        fobj.write("{} {} line {} {}\n".format(counter, mat_num["left_b"], left_b[i], left_b[i+1]))
        counter+=1
    # right boundary
    right_b = [x+(xsize-1) for x in left_b]
    for i in range(0, ysize-1):
        fobj.write("{} {} line {} {}\n".format(counter, mat_num["right_b"], right_b[i], right_b[i+1]))    
        counter+=1
    
def main(fobj, xsize, ysize, alpha):
    # create spaced out points
    n_x = create_point_kernel(size=xsize, alpha=alpha);
    n_y = create_point_kernel(size=ysize,start = 0, end = 8, alpha=alpha);
    # create cartesian coord with the points
    cartesian_product =  cartesian_coord(n_x, n_y)
    # claculating number of quad/line objects in inp
    num_sq_x = math.ceil((((xsize-1)*2)+1)/2) - 1
    total_sq = num_sq_x * (ysize-1)
    total_lines = ((xsize-1) * 2) + ((ysize-1) * 2) 
    # print first line
    fobj.write("{} {} 0 0 0\n".format(xsize*ysize, total_sq+total_lines, 4))
    print_point_coords(fobj, cartesian_product, xsize, ysize) # print coordinates
    print_quad_box(fobj, xsize, ysize) # print quad boxes
    print_boundary_lines(fobj, total_sq, xsize, ysize) # print all boundary lines

if __name__ == '__main__':
    xsize = 50
    ysize = 400
    alpha = 0.005
    with open("channel.inp", "w") as fobj:
        main(fobj,xsize, ysize, alpha)

