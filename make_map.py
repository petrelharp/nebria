#!/usr/bin/python3
import numpy as np
import PIL.Image as Image


H, W = (10, 15)
res = 0.1
p = 0.005
h, w = (int(H / res), int(W / res))
habitat = np.zeros(h * w, dtype=np.uint8)
habitat[np.random.uniform(size=h * w) < p] = 255
habitat = habitat.reshape((h, w))
im = Image.fromarray(habitat)
im.save("habitat.png")


