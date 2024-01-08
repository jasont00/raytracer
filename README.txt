Input Syntax:

NEAR n
LEFT l
RIGHT r
BOTTOM b 
TOP t 
RES x y
SPHERE name posx posy posz sclx scly sclz r g b ka kd ks kr n 
LIGHT name posx posy posz ir ig ib
BACK r g b
AMBIENT ir ig ib 
OUTPUT name

Sphere syntax: name (string), position (posx, posy, posz), scale (sclx, scly, sclz),
color (r, g, b), ambient (ka), diffuse (kb), specular (ks), reflectance (kr), specular exponent (n) (Phong illumination value).

Light syntax: name (string), position (posx, posy, posz), intensity (rgb).

How to Run:
- type the following in console
- ./raytracer.exe ./file-directory/filename.txt