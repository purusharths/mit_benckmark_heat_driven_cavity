import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta
import math

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


def print_point_coords(cc, xsize, ysize):
    # cartesian_coord(n_x, n_y)
    offset = 0
    count = 0
    for y in range(0,ysize):
        # print("\n*** {}\n  ".format(y))
        for i in range(0,xsize):
            tmp = cc[i+offset]
            print("{} {} {} 0".format(count, tmp[0], tmp[1]))
            count+=1
        offset+=xsize    


def print_quad_box(x_max, y_max):
    count = 0
    offset = 0
    n = x_max
    for i in range(0,y_max-1):
        for j in range(0,n-1):
            print("{} -1 quad {} {} {} {} ".format(count, j+offset, (j+offset)+1, ((j+offset)+1)+n, (j+offset)+n))
            count+=1
        offset+=n

def print_boundary_lines(total_sq_count, xsize, ysize, bb_num = 11, tb_num = 12, rb_num = 13, lb_num = 14):
    counter = total_sq_count
    # bottom boundary
    for i in range(xsize-1):
        print("{} {} line {} {}".format(counter, bb_num, i, i+1))
        counter+=1

    # top boundary
    tb_offset = (ysize-1)*xsize
    for i in range(xsize-1):
        print("{} {} line {} {}".format(counter, tb_num, i+tb_offset, (i+1)+tb_offset))
        counter+=1

    # left boundary
    left_b = list(range(0, tb_offset+1, xsize))
    for i in range(0, ysize-1):
        print("{} {} line {} {}".format(counter, lb_num, left_b[i], left_b[i+1]))
        counter+=1
    # print("Right")
    # right boundary
    right_b = [x+(xsize-1) for x in left_b]
    for i in range(0, ysize-1):
        print("{} {} line {} {}".format(counter, rb_num, right_b[i], right_b[i+1]))    
        counter+=1
    

if __name__ == '__main__':
    xsize = 500
    ysize = 800
    x_range = [0,1]
    y_range = [0,8]
    n_x = create_point_kernel(size=xsize, alpha=0.01);
    n_y = create_point_kernel(size=ysize,start = 0, end = 8, alpha=0.02);
    # cartesian_product = list(itertools.product(n_y,n_x))
    # plt.scatter(*zip(*cartesian_product))
    cartesian_product =  cartesian_coord(n_x, n_y)
    num_sq_x = math.ceil((((xsize-1)*2)+1)/2) - 1
    total_sq = num_sq_x * (ysize-1)
    total_lines = ((xsize-1) * 2) + ((ysize-1) * 2) 
    print("{} {} 0 0 0".format(xsize*ysize, total_sq+total_lines, 4))
    print_point_coords(cartesian_product, xsize, ysize)
    print_quad_box(xsize, ysize)
    print_boundary_lines(total_sq, xsize, ysize)

