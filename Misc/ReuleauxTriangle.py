#!/usr/bin/env python
"""ReuleauxTriangle.py."""
from __future__ import division
import math as math
import matplotlib.pyplot as plt


def equalateral(center, sideLength, angle=0):
    """Defining relations for the tringle can be found at http://mathworld.wolfram.com/EquilateralTriangle.html
    Function takes up to 3 parameters
    center: Center point of the triangle.
        - needs to be a tuple (x, y) coordinate pair
    sideLenght: Side length of triangle
    angle: angle anticlockwise from horizontal to one of the corners of the triangle
        - needsto be in degrees

    """
    R = math.sqrt(3) / 3 * sideLength  # Circumference radius
    # r = math.sqrt(3)/6 * sideLength # inner radius
    # h = math.sqrt(3)/2 * sideLength  # Height  h=R + r
    # Coordinates of each apex
    A = (center[0] + R * math.cos(math.radians(angle)), center[1] + R * math.sin(math.radians(angle)))
    B = (center[0] + R * math.cos(math.radians(angle + 120)), center[1] + R * math.sin(math.radians(angle + 120)))
    C = (center[0] + R * math.cos(math.radians(angle + 240)), center[1] + R * math.sin(math.radians(angle + 240)))

    return A, B, C


def reuleauxTriangle(center, sideLength, angle=0):
    A, B, C = equalateral(center, sideLength, angle=0)  # points

    # The equations for the 3 Circles
    # 1:(x-A[0])**2 + (y-A[1])**2 = sideLength**2
    # 2:(x-B[0])**2 + (y-B[1])**2 = sideLength**2
    # 3:(x-C[0])**2 + (y-C[1])**2 = sideLength**2

    # Conditions
    # The edges of the curve will satisfy one of these
    # A = (Ax, Ay) -> A[0] = Ax etc.
    # (x-A[0])**2 + (y-A[1])**2 = sideLength**2 and
    # (x-B[0])**2 + (y-B[1])**2 < sideLength**2 and
    # (x-C[0])**2 + (y-C[1])**2 < sideLength**2
    #                     or
    # (x-A[0])**2 + (y-A[1])**2 < sideLength**2 and
    # (x-B[0])**2 + (y-B[1])**2 = sideLength**2 and
    # (x-C[0])**2 + (y-C[1])**2 < sideLength**2
    #                     or
    # (x-A[0])**2 + (y-A[1])**2 < sideLength**2 and
    # (x-B[0])**2 + (y-B[1])**2 < sideLength**2 and
    # (x-C[0])**2 + (y-C[1])**2 = sideLength**2

    return


# Testing a Triangle
center = (5, -11)
Length = 5
theta = 57

A, B, C = equalateral(center, Length, angle=theta)

print(A, B, C)
plt.figure()
plt.plot(center[0], center[1], "*", label="Center")
plt.plot(A[0], A[1], "h", label="A")
plt.plot(B[0], B[1], "o", label="B")
plt.plot(C[0], C[1], ">", label="C")
plt.plot([A[0], B[0], C[0], A[0]], [A[1], B[1], C[1], A[1]], label="line")
plt.xlim(min([A[0], B[0], C[0]]) - Length, max([A[0], B[0], C[0]]) + Length)
plt.ylim(min([A[1], B[1], C[1]]) - Length, max([A[1], B[1], C[1]]) + Length)
plt.legend()
# plt.show()


# Drawing circles on plot
circle1 = plt.Circle(A, Length, color='r', alpha=0.2)
circle2 = plt.Circle(B, Length, color='b', alpha=0.2)
circle3 = plt.Circle(C, Length, color='g', alpha=0.2, clip_on=False)
fig = plt.gcf()
fig.gca().add_artist(circle1)
fig.gca().add_artist(circle2)
fig.gca().add_artist(circle3)
plt.show()
