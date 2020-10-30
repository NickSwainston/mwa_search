from math import cos, sin
import numpy as np

# grid movements all in rad ----------------------
def left(ra_in, dec_in, fwhm):
    dec_out = dec_in
    ra_out = ra_in - fwhm/cos(dec_in)
    return [ra_out,dec_out]

def right(ra_in, dec_in, fwhm):
    dec_out = dec_in
    ra_out = ra_in + fwhm/cos(dec_in)
    return [ra_out,dec_out]

def up(ra_in, dec_in, fwhm):
    dec_out = dec_in + fwhm / cos(dec_in + np.radians(26.7))**2
    ra_out = ra_in
    return [ra_out,dec_out]

def down(ra_in, dec_in, fwhm):
    dec_out = dec_in - fwhm / cos(dec_in + np.radians(26.7))**2
    ra_out = ra_in
    return [ra_out,dec_out]

def up_left(ra_in, dec_in, fwhm):
    half_fwhm_approx = fwhm/2.#/cos(dec_in)
    dec_out = dec_in + sin(np.radians(60.))*sin(half_fwhm_approx) / sin(np.radians(30.)) \
              / cos(dec_in + np.radians(26.7))**2
    ra_out = ra_in - fwhm/2./cos(dec_out)
    return [ra_out,dec_out]

def up_right(ra_in, dec_in, fwhm):
    half_fwhm_approx = fwhm/2.#/cos(dec_in)
    dec_out = dec_in + sin(np.radians(60.))*sin(half_fwhm_approx) / sin(np.radians(30.)) \
              / cos(dec_in + np.radians(26.7))**2
    ra_out = ra_in + fwhm/2./cos(dec_out)
    return [ra_out,dec_out]

def down_left(ra_in, dec_in, fwhm):
    half_fwhm_approx = fwhm/2.#/cos(dec_in)
    dec_out = dec_in - sin(np.radians(60.))*sin(half_fwhm_approx) / sin(np.radians(30.)) \
              / cos(dec_in + np.radians(26.7))**2
    ra_out = ra_in - fwhm/2./cos(dec_out)
    return [ra_out,dec_out]

def down_right(ra_in, dec_in, fwhm):
    half_fwhm_approx = fwhm/2.#/cos(dec_in)
    dec_out = dec_in - sin(np.radians(60.))*sin(half_fwhm_approx) / sin(np.radians(30.)) \
              / cos(dec_in + np.radians(26.7))**2
    ra_out = ra_in + fwhm/2./cos(dec_out)
    return [ra_out,dec_out]

# ------------------------------------------------

def cross_grid(ra0,dec0,centre_fwhm, loop):
    #start location list [loop number][shape corner (6 for hexagon 4 for square)][number from corner]
    #each item has [ra,dec,fwhm] in radians
    pointing_list = [[[[ra0,dec0,centre_fwhm]]]]
    print("Calculating the tile positions")
    for l in range(loop):
        loop_temp = []
        for c in range(4):
            if c == 0:
                ra,dec =left(pointing_list[l][c][0][0],
                             pointing_list[l][c][0][1],centre_fwhm)
            elif c == 1:
                if l == 0:
                    c = 0
                ra,dec =up(pointing_list[l][c][0][0],
                           pointing_list[l][c][0][1],centre_fwhm)
            elif c == 2:
                if l == 0:
                    c = 0
                ra,dec =right(pointing_list[l][c][0][0],
                              pointing_list[l][c][0][1],centre_fwhm)
            elif c == 3:
                if l == 0:
                    c = 0
                ra,dec =down(pointing_list[l][c][0][0],
                             pointing_list[l][c][0][1],centre_fwhm)
            loop_temp.append([[ra, dec]])
        pointing_list.append(loop_temp)
    return pointing_list


def hex_grid(ra0,dec0,centre_fwhm, loop):
    #start location list [loop number][shape corner (6 for hexagon 4 for square)][number from corner]
    #each item has [ra,dec,fwhm] in radians
    pointing_list = [[[[ra0,dec0,centre_fwhm]]]]
    print("Calculating the tile positions")

    for l in range(loop):
        #different step for each corner
        loop_temp = []
        for c in range(6):
            corner_temp = []
            for n in range(l + 1):
                if l == 0:
                    #First loop so all c = 0
                    if c == 0:
                        ra,dec =left(pointing_list[l][0][n][0],
                                     pointing_list[l][0][n][1],centre_fwhm)
                    elif c == 1:

                        ra,dec =up_left(pointing_list[l][0][n][0],
                                        pointing_list[l][0][n][1],centre_fwhm)
                    elif c == 2:
                        ra,dec =up_right(pointing_list[l][0][n][0],
                                         pointing_list[l][0][n][1],centre_fwhm)
                    elif c == 3:
                        ra,dec =right(pointing_list[l][0][n][0],
                                      pointing_list[l][0][n][1],centre_fwhm)
                    elif c == 4:
                        ra,dec =down_right(pointing_list[l][0][n][0],
                                           pointing_list[l][0][n][1],centre_fwhm)
                    elif c == 5:
                        ra,dec =down_left(pointing_list[l][0][n][0],
                                          pointing_list[l][0][n][1],centre_fwhm)

                elif l == n:
                    #change the 2 for each loop
                    #uses next corner
                    if c == 0:
                        ra,dec =left(pointing_list[l][c+1][0][0],
                                     pointing_list[l][c+1][0][1],centre_fwhm)
                    elif c == 1:
                        ra,dec =up_left(pointing_list[l][c+1][0][0],
                                     pointing_list[l][c+1][0][1],centre_fwhm)
                    elif c == 2:
                        ra,dec =up_right(pointing_list[l][c+1][0][0],
                                     pointing_list[l][c+1][0][1],centre_fwhm)
                    elif c == 3:
                        ra,dec =right(pointing_list[l][c+1][0][0],
                                     pointing_list[l][c+1][0][1],centre_fwhm)
                    elif c == 4:
                        ra,dec =down_right(pointing_list[l][c+1][0][0],
                                     pointing_list[l][c+1][0][1],centre_fwhm)
                    elif c == 5:
                        ra,dec =down_left(pointing_list[l][0][0][0],
                                 pointing_list[l][0][0][1],centre_fwhm)

                else:
                    if c == 0:
                        ra,dec =left(pointing_list[l][c][n][0],
                                     pointing_list[l][c][n][1],centre_fwhm)
                    elif c == 1:

                        ra,dec =up_left(pointing_list[l][c][n][0],
                                        pointing_list[l][c][n][1],centre_fwhm)
                    elif c == 2:
                        ra,dec =up_right(pointing_list[l][c][n][0],
                                         pointing_list[l][c][n][1],centre_fwhm)
                    elif c == 3:
                        ra,dec =right(pointing_list[l][c][n][0],
                                      pointing_list[l][c][n][1],centre_fwhm)
                    elif c == 4:
                        ra,dec =down_right(pointing_list[l][c][n][0],
                                           pointing_list[l][c][n][1],centre_fwhm)
                    elif c == 5:
                        ra,dec =down_left(pointing_list[l][c][n][0],
                                          pointing_list[l][c][n][1],centre_fwhm)
                corner_temp.append([ra,dec])
            loop_temp.append(corner_temp)
        pointing_list.append(loop_temp)
    return pointing_list


def square_grid(ra0,dec0,centre_fwhm, loop):
    #start location list [loop number][shape corner (4 for square)][number from corner]
    #each item has [ra,dec,fwhm] in radians
    pointing_list = [[[[ra0,dec0,centre_fwhm]]]]
    print("Calculating the tile positions")

    for l in range(loop):
        #different step for each corner
        loop_temp = []
        for c in range(4):
            corner_temp = []
            for n in range((l + 1) * 2):
                if n == 0:
                    #grab from previous corner
                    cfrom = (c + 3)%4
                    if l == 0:
                        #if first loop set all corners to zero
                        cfrom = 0
                    #and last number from corner
                    nfrom = l * 2 - 1
                    #First loop so all c = 0
                    if c == 0:
                        ra,dec = left(pointing_list[l][cfrom][nfrom][0],
                                      pointing_list[l][cfrom][nfrom][1],centre_fwhm)
                    elif c == 1:
                        ra,dec = up(pointing_list[l][cfrom][nfrom][0],
                                    pointing_list[l][cfrom][nfrom][1],centre_fwhm)
                    elif c == 2:
                        ra,dec = right(pointing_list[l][cfrom][nfrom][0],
                                       pointing_list[l][cfrom][nfrom][1],centre_fwhm)
                    elif c == 3:
                        ra,dec = down(pointing_list[l][cfrom][nfrom][0],
                                      pointing_list[l][cfrom][nfrom][1],centre_fwhm)

                else:
                    #moves to the edges
                    if c == 0:
                        ra,dec = up(corner_temp[n-1][0],
                                    corner_temp[n-1][1],centre_fwhm)
                    elif c == 1:
                        ra,dec = right(corner_temp[n-1][0],
                                       corner_temp[n-1][1],centre_fwhm)
                    elif c == 2:
                        ra,dec = down(corner_temp[n-1][0],
                                      corner_temp[n-1][1],centre_fwhm)
                    elif c == 3:
                        ra,dec = left(corner_temp[n-1][0],
                                      corner_temp[n-1][1],centre_fwhm)


                corner_temp.append([ra,dec])
            loop_temp.append(corner_temp)
        pointing_list.append(loop_temp)
    return pointing_list


def get_grid(ra, dec, grid_sep, loop, grid_type='hex'):
    """
    ra: Right Acension in radians
    dec: Declination in radians
    grid_sep: seperation between grid pointings in radians
    loop: number of pointing loops
    grid_type: Possible grid types from ['hex', 'cross', 'squaire']

    return [rads, decds]
    RAs and Decs in degrees
    """
    #calc grid positions
    if grid_type == 'hex':
        pointing_list = hex_grid(   ra, dec, grid_sep, loop)
    elif grid_type == 'cross':
        pointing_list = cross_grid( ra, dec, grid_sep, loop)
    elif grid_type == 'square':
        pointing_list = square_grid(ra, dec, grid_sep, loop)
    else:
        print("Unrecognised grid type. Exiting.")
        quit()
    #TODO add square

    rads = []; decds = []

    print("Converting ra dec to degrees")
    for loop in pointing_list:
        for corner in loop:
            for num in corner:
                #format grid pointings
                rad = np.degrees(num[0])
                decd = np.degrees(num[1])

                if (decd < 90.) and (decd > -90.):
                    # Only include ra and dec within the real decs
                    rads.append(rad)
                    decds.append(decd)
    return rads, decds
