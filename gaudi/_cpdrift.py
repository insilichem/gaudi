#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
# 
# https://github.com/insilichem/gaudi
#
# Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############

"""
Coherent Point Drift (affine and rigid) Python2/3 implementation,
adapted from `kwohlfahrt's <https://github.com/kwohlfahrt/coherent-point-drift>`_.

Only 3D points are supported in this version.

Depends on:

- Python 2.7, 3.4+
- Numpy
- Matplotlib (plotting only)
"""

from __future__ import division
import numpy as np
import math
from itertools import product, repeat, islice
from functools import reduce


def coherent_point_drift(X, Y, w=0.0, B=None, guess_steps=5, max_iterations=20, 
                         method='affine'):
    register, transform = {'affine': (affine_cpd, affine_xform),
                           'rigid':  (rigid_cpd, rigid_xform)}[method]

    best_xforms, best_transformed_Y, best_rmsd = None, None, 1000
    for rotation in spaced_rotations(guess_steps):
        RB = rotation_matrix(*rotation)
        xforms = last(islice(register(X, Y, w, RB), max_iterations))
        transformed_Y = transform(Y, *xforms)
        rmsd = RMSD(X, transformed_Y)
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_xforms = xforms
            best_transformed_Y = transformed_Y
    return best_transformed_Y, best_xforms, best_rmsd


def rigid_cpd(X, Y, w=0.0, R=None):
    D = X.shape[1]
    N, M = len(X), len(Y)

    R = np.eye(D) if R is None else R
    s = std(X) / std(rigid_xform(Y, R=R))
    t = X.mean(axis=0) - rigid_xform(Y, R=R, s=s).mean(axis=0)
    sigma_sq = pairwise_sqdist(rigid_xform(Y, R, t, s), X).sum() / (D * M * N)

    old_exceptions = np.seterr(divide='ignore', over='ignore', under='ignore', 
                               invalid='raise')
    while True:
        try:
            transformed_Y = rigid_xform(Y, R, t, s)
            P, N_p, mu_x, mu_y, X_hat, Y_hat = common_steps(X, Y, transformed_Y, 
                                                            w, sigma_sq)
        except FloatingPointError:
            np.seterr(**old_exceptions)
            break

        A = X_hat.T.dot(P.T).dot(Y_hat)
        U, _, VT = np.linalg.svd(A)
        C = np.eye(D)
        C[-1, -1] = np.linalg.det(U.dot(VT))
        R = U.dot(C).dot(VT)
        s = np.trace(A.T.dot(R)) / np.trace(Y_hat.T.dot(np.diag(P.sum(axis=1))).dot(Y_hat))
        t = mu_x - s * R.dot(mu_y)
        sigma_sq = (np.trace(X_hat.T.dot(np.diag(P.T.sum(axis=1))).dot(X_hat))
                    - s * np.trace(A.T.dot(R))) / (N_p * D)

        yield R, t, s


def affine_cpd(X, Y, w=0.0, B=None):
    D = X.shape[1]
    N = len(X)
    M = len(Y)

    B = np.eye(D) if B is None else B
    s = (std(X) / std(affine_xform(Y, B=B)))  # scale
    B = s * B
    t = X.mean(axis=0) - affine_xform(Y, B=B).mean(axis=0)

    sigma_sq = pairwise_sqdist(affine_xform(Y, B, t), X).sum() / (D * M * N)
    old_exceptions = np.seterr(divide='ignore', over='ignore', under='ignore', 
                               invalid='raise')
    while True:
        try:
            transformed_Y = affine_xform(Y, B, t)
            P, N_p, mu_x, mu_y, X_hat, Y_hat = common_steps(X, Y, transformed_Y, w, 
                                                            sigma_sq)
        except FloatingPointError:
            np.seterr(**old_exceptions)
            break

        B = (X_hat.T.dot(P.T).dot(Y_hat)
             .dot(np.linalg.inv(Y_hat.T.dot(np.diag(P.sum(axis=1))).dot(Y_hat))))
        t = mu_x - B.dot(mu_y)
        sigma_sq = (np.trace(X_hat.T.dot(np.diag(P.T.sum(axis=1))).dot(X_hat))
                    - np.trace(X_hat.T.dot(P.T).dot(Y_hat).dot(B.T))) / (N_p * D)

        yield B, t


def common_steps(X, Y, Y_, w, sigma_sq):
    # E step
    D = X.shape[1]
    N = len(X)
    M = len(Y_)

    dist = pairwise_sqdist(Y_, X)
    overlap = np.exp(-dist / (2 * sigma_sq))

    # The original algorithm expects unit variance,
    # so normalize (2πς**2)**D/2 to compensate
    P = overlap / (overlap.sum(axis=0)
                   + (2 * math.pi * sigma_sq) ** (D / 2)
                   * w / (1-w) * M / N / std(X) ** D)

    # M-step
    N_p = P.sum()
    mu_x = 1 / N_p * X.T.dot(P.T.sum(axis=1))
    mu_y = 1 / N_p * Y.T.dot(P.sum(axis=1))
    X_hat = X - mu_x.T
    Y_hat = Y - mu_y.T

    return P, N_p, mu_x, mu_y, X_hat, Y_hat


def last(sequence):
    return reduce(lambda acc, x: x, sequence)


def std(x):
    return np.sqrt(np.var(x, axis=0).mean())


def rigid_xform(X, R=np.array(1), t=0.0, s=1.0):
    return s * R.dot(X.T).T + t


def affine_xform(X, B=np.array(1), t=0):
    return B.dot(X.T).T + t


def spaced_rotations(N):
    # Ken Shoemake
    # Graphics Gems III, pp 124-132
    for x_theta in product(frange(0, 1, 1/N), 
                           *repeat(frange(0, 2*math.pi, 2*math.pi/N), 2)):
        X, theta = x_theta[0], x_theta[1:]
        R0, R1 = math.sqrt(1-X), math.sqrt(X)
        yield Quaternion(math.sin(theta[0]) * R0, math.cos(theta[0]) * R0,
                         math.sin(theta[1]) * R1, math.cos(theta[1]) * R1).axis_angle


def rotation_matrix(*angles):
    theta, (x, y, z) = angles
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c+x*x*(1-c),   x*y*(1-c)-z*s, (1-c)*x*z+y*s],
                     [y*x*(1-c)+z*s, c+y*y*(1-c),   y*z*(1-c)-x*s],
                     [z*x*(1-c)-y*s, z*y*(1-c)+x*s, c+z*z*(1-c)]])


def pairwise_sqdist(X, Y):
    # R[i, j] = distance(X[i], Y[j]) ** 2
    return ((X[:, None, :] - Y[None, :, :]) ** 2).sum(axis=2)


def RMSD(X, Y):
    dist = pairwise_sqdist(X, Y)
    # Minimum RMSD for each point in X
    min_rmsd = np.sqrt(dist.min(axis=1).mean())
    # Normalize for scale
    return min_rmsd / std(X)


def plot(x, y, t):
    """
    Plot the initial datasets and registration results.

    Parameters
    ----------
    x : ndarray
        The static shape that y will be registered to. Expected array shape is 
        [n_points_x, n_dims]
    y : ndarray
        The moving shape. Expected array shape is [n_points_y, n_dims]. 
        Note that n_dims should be equal for x and y, but n_points does 
        not need to match.
    t : ndarray
        The transformed version of y. Output shape is [n_points_y, n_dims].
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    ax1 = Axes3D(plt.figure(1))
    ax1.plot(x[:, 0], x[:, 1], x[:, 2], 'yo')
    ax1.plot(y[:, 0], y[:, 1], y[:, 2], 'r+')
    ax1.set_title("Before registration", fontdict=None, loc='center')
    ax2 = Axes3D(plt.figure(2))
    ax2.plot(x[:, 0], x[:, 1], x[:, 2], 'yo')
    ax2.plot(t[:, 0], t[:, 1], t[:, 2], 'r+')
    ax2.set_title("After registration", fontdict=None, loc='center')
    plt.show()


class Quaternion(object):

    def __init__(self, s, i, j, k):
        self.s = s
        self.i = i
        self.j = j
        self.k = k

    @classmethod
    def fromAxisAngle(cls, v, theta):
        if len(v) != 3:
            raise ValueError("v must be a 3D vector.")

        theta = theta % (2 * math.pi) / 2
        try:
            # Numpy-compatible array
            v = v / v.norm() * math.sin(theta)
        except AttributeError:
            norm = math.sqrt(sum(map(lambda x: x ** 2, v)))
            v = map(lambda x: x / norm * math.sin(theta), v)
        return cls(math.cos(theta), *v)

    def __iter__(self):
        for i in (self.s, self.i, self.j, self.k):
            yield i

    def __repr__(self):
        return ("Quaternion({}, {}i, {}j, {}k)"
                .format(self.s, self.i, self.j, self.k))

    def __eq__(self, other):
        return (self.s == other.s
                and self.i == other.i
                and self.j == other.j
                and self.k == other.k)

    def __hash__(self):
        return hash(self.s, self.i, self.j, self.k)

    def __add__(self, other):
        return Quaternion(self.s + other.s,
                          self.i + other.i,
                          self.j + other.j,
                          self.k + other.k)

    def __sub__(self, other):
        return Quaternion(self.s - other.s,
                          self.i - other.i,
                          self.j - other.j,
                          self.k - other.k)

    def __mul__(self, other):
        try:
            return Quaternion(
                self.s * other.s - self.i * other.i - self.j * other.j - self.k * other.k,
                self.s * other.i + self.i * other.s + self.j * other.k - self.k * other.j,
                self.s * other.j - self.i * other.k + self.j * other.s + self.k * other.i,
                self.s * other.k + self.i * other.j - self.j * other.i + self.k * other.s)
        except AttributeError:
            return Quaternion(other * self.s, other * self.i, 
                              other * self.j, other * self.k)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return Quaternion(self.s / other, self.i / other, self.j / other, self.k / other)

    def __abs__(self):
        return self.norm()

    def __round__(self, ndigits=0):
        return Quaternion(*(round(i, ndigits) for i in (self.s, self.i, self.j, self.k)))

    def norm(self):
        return math.sqrt(sum(map(lambda x: x ** 2), (self.s, self.i, self.j, self.k)))

    def reciprocal(self):
        return self.conjugate() / self.norm() ** 2

    def conjugate(self, other=None):
        if not other:
            return Quaternion(self.s, -self.i, -self.j, -self.k)
        else:
            return self * other * self.reciprocal()

    def unit(self):
        return self / self.norm()

    def matrix(self):
        r, i, j, k = self
        return [[1 - 2*j*j - 2*k*k, 2*(i*j - k*r), 2*(i*k + j*r)],
                [2*(i*j + k*r), 1 - 2*i*i - 2*k*k, 2*(j*k - i*r)],
                [2*(i*k - j*r), 2*(j*k + i*r), 1 - 2*i*i - 2*j*j]]

    @property
    def vector(self):
        return self.i, self.j, self.k

    @vector.setter
    def vector(self, *v):
        self.i, self.j, self.k = v

    @property
    def axis_angle(self):
        theta = 2 * math.acos(self.s)
        try:
            v = tuple(map(lambda x: x / math.sin(theta / 2), self.vector))
        except ZeroDivisionError:
            # No rotation, so v doesn't matter
            v = (1, 0, 0)
        return theta, v


class frange(object):

    def __init__(self, start, stop, step):
        self.start = start
        self.stop = stop
        self.step = step

    def __len__(self):
        return int((self.stop - self.start) / self.step)

    def __getitem__(self, idx):
        if not 0 <= idx < (self.stop - self.start) / (self.step):
            raise IndexError("Index {} out of range".format(idx))
        return self.start + self.step * idx

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]
