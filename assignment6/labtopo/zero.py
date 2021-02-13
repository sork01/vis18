#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import sys


class vec2:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y
    
    def length(self):
        return math.sqrt(self.x*self.x + self.y*self.y)
    
    def __getitem__(self, key):
        if key == 0:
            return self.x
        else:
            return self.y
    
    def __setitem__(self, key, value):
        if key == 0:
            self.x = value
        else:
            self.y = value
    
    def __add__(self, other):
        if type(other) == vec2:
            return vec2(self.x + other.x, self.y + other.y)
        else:
            return vec2(self.x + other, self.y + other)
    
    def __radd__(self, other):
        if type(other) == vec2:
            return vec2(self.x + other.x, self.y + other.y)
        else:
            return vec2(self.x + other, self.y + other)
    
    def __sub__(self, other):
        if type(other) == vec2:
            return vec2(self.x - other.x, self.y - other.y)
        else:
            return vec2(self.x - other, self.y - other)
    
    def __mul__(self, other):
        if type(other) == vec2:
            return vec2(self.x * other.x, self.y * other.y)
        else:
            return vec2(self.x * other, self.y * other)
    
    def __rmul__(self, other):
        if type(other) == vec2:
            return vec2(self.x * other.x, self.y * other.y)
        else:
            return vec2(self.x * other, self.y * other)
    
    def __truediv__(self, other):
        return vec2(self.x / other, self.y / other)
    
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
    
    def __str__(self):
        return "vec2({:.8f}, {:.8f})".format(self.x, self.y)
    
    def __repr__(self):
        return str(self)



class bounds:
    def __init__(self, bottom_left, width, height):
        self.bottom_left = bottom_left
        self.width = width
        self.height = height
    
    def getCornerPoints(self):
        return (
            vec2(self.bottom_left.x, self.bottom_left.y),
            vec2(self.bottom_left.x + self.width, self.bottom_left.y),
            vec2(self.bottom_left.x, self.bottom_left.y + self.height),
            vec2(self.bottom_left.x + self.width, self.bottom_left.y + self.height))
    
    def getCenter(self):
        return self.bottom_left + vec2(self.width / 2.0, self.height / 2.0)
    
    def __str__(self):
        return "bounds({}, {}, {})".format(self.bottom_left, self.width, self.height)
    
    def __repr__(self):
        return str(self)



class clsvol:
    def __init__(self, basis, voldata):
        self.basis = basis
        self.voldata = voldata
    
    def getDimensions(self):
        return vec2(len(self.voldata), len(self.voldata[0]))
    
    def getBasis(self):
        return self.basis



def sampleFromField(vol, position):
    # auto dims = vr->getDimensions();
    dims = vol.getDimensions()
    
    # auto base = vol->getBasis();
    base = vol.getBasis()
    
    # vec2 cellSize(base[0][0] / dims[0], base[1][1] / dims[1]);
    cellSize = vec2(base[0][0] / dims[0], base[1][1] / dims[1]);
    
    # // Sampled outside the domain!
    # if (position[0] < 0 || position[0] > dims[0] - 1 || position[1] < 0 ||
        # position[1] > dims[1] - 1)
    # {
        # return vec2(0, 0);
    # }
    if position[0] < 0 or position[0] > dims[0] - 1 or position[1] < 0 or position[1] > dims[1] - 1:
        return vec2(0,0)
    
    # size2_t p((size_t)position[0], (size_t)position[1]);
    p = [int(position[0]), int(position[1])]
    
    # // Leads to accessing only inside the volume
    # // Coefficients computation takes care of using the correct values
    # for (int d = 0; d < 2; ++d)
        # p[d] = std::min(p[d], dims[d] - 2);
    for d in range(2):
        p[d] = min(p[d], dims[d] - 2)
    
    # const auto f00 = vr->getAsDVec2(size3_t(p[0], p[1], 0));
    # const auto f10 = vr->getAsDVec2(size3_t(p[0] + 1, p[1], 0));
    # const auto f01 = vr->getAsDVec2(size3_t(p[0], p[1] + 1, 0));
    # const auto f11 = vr->getAsDVec2(size3_t(p[0] + 1, p[1] + 1, 0));
    
    f00 = vol.voldata[p[0]    ][p[1]]
    f10 = vol.voldata[p[0] + 1][p[1]]
    f01 = vol.voldata[p[0]    ][p[1] + 1]
    f11 = vol.voldata[p[0] + 1][p[1] + 1]
    
    # const float x = position[0] - p[0];
    # const float y = position[1] - p[1];
    x = position[0] - p[0]
    y = position[1] - p[1]
    
    # vec2 f;
    f = vec2()
    
    # for (int i = 0; i < 2; i++)
    # {
    for i in range(2):
        f[i] = f00[i] * (1 - x) * (1 - y) + f01[i] * (1 - x) * y + f10[i] * x * (1 - y) + f11[i] * x * y

        # // Bring vector back to grid space.
        f[i] /= cellSize[i]
    # }

    # return f;
    return f



def ANoise2CT4_dat():
    voldata = list()
    for i in range(5):
        voldata.append([0] * 5)
    
    voldata[0][0] = vec2(0.965104,-0.976884)
    voldata[0][1] = vec2(0.999949,-1)
    voldata[0][2] = vec2(-0.196601,0.992637)
    voldata[0][3] = vec2(0.999882,-1)
    voldata[0][4] = vec2(0.886238,-0.78125)
    voldata[1][0] = vec2(0.833346,-0.57363)
    voldata[1][1] = vec2(0.767033,-0.273978)
    voldata[1][2] = vec2(-0.47048,0.767801)
    voldata[1][3] = vec2(-0.951299,-0.955765)
    voldata[1][4] = vec2(0.629356,0.320677)
    voldata[2][0] = vec2(-0.888844,-0.790239)
    voldata[2][1] = vec2(-0.158863,0.996859)
    voldata[2][2] = vec2(0.979501,-0.991886)
    voldata[2][3] = vec2(-0.938846,-0.931424)
    voldata[2][4] = vec2(-0.17767,0.995087)
    voldata[3][0] = vec2(-0.793876,-0.397837)
    voldata[3][1] = vec2(0.401217,0.874827)
    voldata[3][2] = vec2(0.517895,0.665508)
    voldata[3][3] = vec2(0.86423,-0.700173)
    voldata[3][4] = vec2(0.377817,0.90112)
    voldata[4][0] = vec2(0.783264,-0.349071)
    voldata[4][1] = vec2(0.763329,-0.256811)
    voldata[4][2] = vec2(0.227892,0.986719)
    voldata[4][3] = vec2(0.817933,-0.506373)
    voldata[4][4] = vec2(0.852376,-0.65313)
    
    return clsvol([(2,0), (0,2)], voldata)



def rk4(vol, position, stepsize):
    v1 = sampleFromField(vol, position)
    v2 = sampleFromField(vol, position + (stepsize/2.0)*v1)
    v3 = sampleFromField(vol, position + (stepsize/2.0)*v2)
    v4 = sampleFromField(vol, position + stepsize*v3)
    
    return position + stepsize * ((v1/6.0) + (v2/3.0) + (v3/3.0) + (v4/6.0))



def possiblezero(vol, thebounds): # vec2[4] bounds p00 p10 p01 p11
    corners = thebounds.getCornerPoints()
    f00 = sampleFromField(vol, corners[0])
    f10 = sampleFromField(vol, corners[1])
    f01 = sampleFromField(vol, corners[2])
    f11 = sampleFromField(vol, corners[3])
    
    s00 = vec2(1 if f00[0] > 0 else 0, 1 if f00[1] > 0 else 0)
    s10 = vec2(1 if f10[0] > 0 else 0, 1 if f10[1] > 0 else 0)
    s01 = vec2(1 if f01[0] > 0 else 0, 1 if f01[1] > 0 else 0)
    s11 = vec2(1 if f11[0] > 0 else 0, 1 if f11[1] > 0 else 0)
    
    # print(s00)
    # print(s10)
    # print(s01)
    # print(s11)
    
    vz = vec2(1,1)
    
    if s00 + s10 == vz or s00 + s01 == vz or s00 + s11 == vz or s10 + s01 == vz or s10 + s11 == vz or s01 + s11 == vz:
        return True
    else:
        return False



def decomposition_start(vol, initialsubs, threshold):
    subbounds = list()
    
    dims = vol.getDimensions()
    
    subdim = math.sqrt(initialsubs)
    
    if subdim - int(subdim) != 0:
        raise ValueError("initialsubs must be an integer square")
    
    subsize = vec2((dims.x - 1)  / subdim, (dims.y - 1) / subdim)
    subbounds = list()
    
    xoffset = 0.0
    yoffset = 0.0
    
    for ystep in range(int(subdim)):
        xoffset = 0.0
        
        for xstep in range(int(subdim)):
            # print(vec2(xoffset, yoffset))
            # subbounds.append( (vec2(xoffset, yoffset), vec2(xoffset + subsize.x, yoffset), vec2(xoffset, yoffset + subsize.y), vec2(xoffset + subsize.x, yoffset + subsize.y)) )
            subbounds.append(bounds(vec2(xoffset, yoffset), subsize.x, subsize.y))
            xoffset += subsize.x
        
        yoffset += subsize.y
    
    results = list()
    
    for thebounds in subbounds:
        result = decomposition(vol, thebounds, threshold)
        
        if result is not None:
            results.append(result)
    
    return results



def decomposition(vol, thebounds, threshold):
    halfwidth = thebounds.width * 0.5
    halfheight = thebounds.height * 0.5
    
    if possiblezero(vol, thebounds):
        if thebounds.width <= threshold:
            return thebounds.getCenter()
        else:
            test = decomposition(vol, bounds(thebounds.bottom_left, halfwidth, halfheight), threshold)
            
            if test is not None:
                return test
                
            test = decomposition(vol, bounds(thebounds.bottom_left + vec2(halfwidth, 0), halfwidth, halfheight), threshold)
            
            if test is not None:
                return test
                
            test = decomposition(vol, bounds(thebounds.bottom_left + vec2(0, halfheight), halfwidth, halfheight), threshold)
            
            if test is not None:
                return test
                
            test = decomposition(vol, bounds(thebounds.bottom_left + vec2(halfwidth, halfheight), halfwidth, halfheight), threshold)
            
            if test is not None:
                return test
            
            return None
    else:
        return None





vol = ANoise2CT4_dat()
print(decomposition_start(vol, 64, 0.00001))


### ALL POSSIBLE ZEROS?
# (0,0)+(1,1) = (1,1)
# (0,1)+(1,0) = (1,1)
# (1,0)+(0,1) = (1,1)
###



### FIND THE SINK IN ANOISE2CT4 NEAR 0.5,2.5
# vol = ANoise2CT4_dat()

# oldpos = vec2(0.5, 2.5)

# for i in range(10000):
  # newpos = rk4(vol, oldpos, 0.01)
  # sample = sampleFromField(vol, newpos)
  # print(newpos, sample.length())
  # oldpos = newpos
###

